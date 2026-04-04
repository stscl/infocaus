/******************************************************************************
 * File: surd.hpp
 *
 * Synergistic–Unique–Redundant Decomposition (SURD) of mutual
 * information between a target variable and multiple source variables.
 *
 * The SURD framework decomposes the mutual information structure into:
 *
 *   0 = Redundant information
 *       Information shared among multiple variables.
 *
 *   1 = Unique information
 *       Information uniquely provided by a single variable.
 *
 *   2 = Synergistic information
 *       Information that only emerges from combinations of variables.
 *
 *   3 = Information loss
 *       Remaining uncertainty in the target after conditioning on all sources.
 *
 * The algorithm operates as follows:
 *
 *   1. Generate all combinations of source variables up to max_order.
 *   2. Compute mutual information I(Y ; X_set) for each combination.
 *   3. Group combinations by order (number of variables).
 *   4. Sort each group by mutual information.
 *   5. Perform ladder-style decomposition:
 *
 *        Order 1:
 *           Incremental differences determine redundant vs unique
 *           contributions.
 *
 *        Higher orders:
 *           Information exceeding the maximum lower-order mutual
 *           information is attributed to synergy.
 *
 *   6. Optionally normalize contributions so they sum to one.
 *   7. Append information loss as the final component.
 *
 * Two data types are supported:
 *
 *   Discrete data
 *      Mutual information computed using joint entropy estimators.
 *
 *   Continuous data
 *      Mutual information estimated via the KSG k-nearest neighbor method.
 *
 * Input data format:
 *
 *   Row 0  : target variable
 *   Row 1+ : source variables
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ******************************************************************************/

#ifndef INFOXTR_SURD_HPP
#define INFOXTR_SURD_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infoxtr/combn.hpp"
#include "infoxtr/numericutils.hpp"
#include "infoxtr/infotheo.hpp"
#include "infoxtr/ksginfo.hpp"
#include <RcppThread.h>

namespace infoxtr 
{

namespace surd 
{

    using DiscMat = std::vector<std::vector<uint64_t>>;
    using ContMat = std::vector<std::vector<double>>;

    /***********************************************************
     * Result structure
     ***********************************************************/
    struct SURDRes
    {
        std::vector<double> values;
        std::vector<uint8_t> types;
        std::vector<std::vector<size_t>> var_indices;

        size_t size() const noexcept { return values.size(); }
    };

    /***********************************************************
     * Synergistic-Unique-Redundant Decomposition Utilities
     ***********************************************************/
    inline SURDRes mutualinfo_decomposition(
        const std::vector<std::vector<size_t>>& combs,
        const std::vector<double>& mi,
        bool normalize = false)
    {
        struct Entry
        {
            size_t idx;
            size_t order;
            double mi;
        };

        const size_t n = combs.size();

        std::vector<Entry> entries(n);

        size_t max_order = 0;

        for (size_t i = 0; i < n; ++i)
        {
            entries[i] = {i, combs[i].size(), mi[i]};
            max_order = std::max(max_order, combs[i].size());
        }

        /***********************************************************
         * Group by order
         ***********************************************************/
        std::vector<std::vector<Entry*>> groups(max_order + 1);

        for (auto& e : entries)
            groups[e.order].push_back(&e);

        for (auto& g : groups)
        {
            std::sort(
                g.begin(),
                g.end(),
                [](Entry* a, Entry* b)
                {   
                    if (!infoxtr::numericutils::doubleNearlyEqual(a->mi, b->mi))
                        return a->mi < b->mi;

                    return a->idx < b->idx;
                });
        }

        auto get_max = [&](size_t m)
        {
            if (m >= groups.size() || groups[m].empty())
                return 0.0;

            return groups[m].back()->mi;
        };

        SURDRes result;
        result.values.reserve(combs.size());
        result.types.reserve(combs.size());
        result.var_indices.reserve(combs.size());

        /***********************************************************
         * Order 1 decomposition
         ***********************************************************/
        if (!groups[1].empty())
        {
            double prev = 0.0;

            for (size_t i = 0; i < groups[1].size(); ++i)
            {
                auto* e = groups[1][i];

                double delta = e->mi - prev;

                if (!infoxtr::numericutils::doubleNearlyEqual(delta,0.0) 
                    && delta < 0.0 )
                {
                    delta = 0.0;
                }

                result.values.push_back(delta);

                if (i == groups[1].size() - 1)
                    result.types.push_back(1);  // unique
                else
                    result.types.push_back(0);  // redundant

                result.var_indices.push_back(combs[e->idx]);

                prev = e->mi;
            }
        }

        /***********************************************************
         * Higher order synergy
         ***********************************************************/
        for (size_t m = 2; m <= max_order; ++m)
        {
            if (groups[m].empty())
                continue;

            double max_prev = get_max(m - 1);

            for (size_t i = 0; i < groups[m].size(); ++i)
            {
                auto* e = groups[m][i];

                double prev =
                    (i > 0) ? groups[m][i - 1]->mi : 0.0;

                double delta = 0.0;

                if (!infoxtr::numericutils::doubleNearlyEqual(e->mi, max_prev)\
                     && e->mi > max_prev)
                {
                    if (prev >= max_prev)
                        delta = e->mi - prev;
                    else
                        delta = e->mi - max_prev;
                }
                
                result.values.push_back(delta);
                result.types.push_back(2);  // synergistic
                result.var_indices.push_back(combs[e->idx]);
            }
        }

        /***********************************************************
        * Optional normalization
        ***********************************************************/
        if (normalize)
        {
            double sum = 0.0;

            for (size_t i = 0; i < result.values.size(); ++i)
                if (result.types[i] != 3)
                    sum += result.values[i];

            if (!infoxtr::numericutils::doubleNearlyEqual(sum, 0.0))
            {
                for (size_t i = 0; i < result.values.size(); ++i)
                    if (result.types[i] != 3)
                        result.values[i] /= sum;
            }
        }

        return result;
    }

    /***************************************************************
     * Synergistic-Unique-Redundant Decomposition for Discrete Data
     ***************************************************************/
    inline SURDRes surd(
        const DiscMat& mat,
        size_t max_order = std::numeric_limits<size_t>::max(),
        size_t threads = 1,
        double base = 2.0,
        bool normalize = false)
    {
        if (threads == 0) threads = 1;
        size_t hw = std::thread::hardware_concurrency();
        if (hw > 0) threads = std::min(threads, hw);

        if (mat.size() < 2)
            throw std::invalid_argument("[SURD] SURD needs >2 variables");

        const size_t n_sources = mat.size() - 1;
        max_order = std::min(max_order, n_sources);

        // Construct variable combination vector
        std::vector<size_t> ag_idx(n_sources);
        std::iota(ag_idx.begin(), ag_idx.end(), 1);
        const std::vector<std::vector<size_t>> combs =
            infoxtr::combn::genSubsets(ag_idx, max_order);

        // Compute joint entropies
        std::vector<double> H_sources(combs.size(), 
                                      std::numeric_limits<double>::quiet_NaN());
        std::vector<double> H_joints(combs.size(), 
                                     std::numeric_limits<double>::quiet_NaN());
        double H_target = infoxtr::infotheo::je(mat, {0}, base, true);

        if (threads <= 1) 
        {   
            for (size_t i = 0; i < combs.size(); ++i)
            {   
                std::vector<size_t> joint_idx = {0};
                joint_idx.insert(joint_idx.end(), combs[i].begin(), combs[i].end());
                H_sources[i] = infoxtr::infotheo::je(mat, combs[i], base, true);
                H_joints[i] = infoxtr::infotheo::je(mat, joint_idx, base, true);
            }
        } 
        else 
        {
            RcppThread::parallelFor(0, combs.size()), [&](size_t i) {
                std::vector<size_t> joint_idx = {0};
                joint_idx.insert(joint_idx.end(), combs[i].begin(), combs[i].end());
                H_sources[i] = infoxtr::infotheo::je(mat, combs[i], base, true);
                H_joints[i] = infoxtr::infotheo::je(mat, joint_idx, base, true);
            }
        }

        // Compute mutual information
        std::vector<double> mi_combs(combs.size(), 
                                     std::numeric_limits<double>::quiet_NaN());
        
        for (size_t i = 0; i < combs.size(); ++i)
        {   
            mi_combs[i] = H_target + H_sources[i] - H_joints[i];
        }
        
        // Decompose mutual information
        SURDRes result = mutualinfo_decomposition(combs, mi_combs, normalize);

        // Information loss
        double leak = 0.0;
        if (!infoxtr::numericutils::doubleNearlyEqual(H_target, 0.0))
        {
            if (max_order < n_sources) 
            {
                std::vector<size_t> joint_idx = {0};
                joint_idx.insert(joint_idx.end(), ag_idx.begin(), ag_idx.end());
                leak = (infoxtr::infotheo::je(mat, joint_idx, base, true) - 
                        infoxtr::infotheo::je(mat, ag_idx, base, true)) / H_target;
            }
            else 
            {   
                leak = (H_joints.back() - H_sources.back()) / H_target;
            }
        }
        leak = std::max(0.0, std::min(1.0, leak));

        result.values.push_back(leak);
        result.types.push_back(3);
        result.var_indices.push_back(ag_idx);

        return result;
    }

    /*****************************************************************
     * Synergistic-Unique-Redundant Decomposition for Continuous Data
     *****************************************************************/
    inline SURDRes surd(
        const ContMat& mat,
        size_t max_order = std::numeric_limits<size_t>::max(),
        size_t k = 3,
        size_t alg = 0,
        size_t threads = 1,
        double base = 2.0,
        bool normalize = false)
    {
        if (threads == 0) threads = 1;
        size_t hw = std::thread::hardware_concurrency();
        if (hw > 0) threads = std::min(threads, hw);

        if (mat.size() < 2)
            throw std::invalid_argument("[SURD] SURD needs >2 variables");
        if (mat[0].size() < k + 1)
            throw std::invalid_argument("[SURD] SURD needs >k observations");
        
        const size_t n_sources = mat.size() - 1;
        max_order = std::min(max_order, n_sources);

        // Construct variable combination vector
        std::vector<size_t> ag_idx(n_sources);
        std::iota(ag_idx.begin(), ag_idx.end(), 1);
        const std::vector<std::vector<size_t>> combs =
            infoxtr::combn::genSubsets(ag_idx, max_order);

        // Compute mutual information
        std::vector<double> mi_combs(combs.size(), 
                                     std::numeric_limits<double>::quiet_NaN());

        if (threads <= 1) 
        {   
            for (size_t i = 0; i < combs.size(); ++i)
            {   
                mi_combs[i] = infoxtr::ksginfo::mi(
                    mat, {0}, combs[i], k, alg, base, false);
            }
        } 
        else 
        {
            RcppThread::parallelFor(0, combs.size()), [&](size_t i) {
                mi_combs[i] = infoxtr::ksginfo::mi(
                    mat, {0}, combs[i], k, alg, base, false);
            }
        }
        
        // Decompose mutual information
        SURDRes result = mutualinfo_decomposition(combs, mi_combs, normalize);

        // Information loss
        double leak = infoxtr::ksginfo::ce(mat, {0}, ag_idx, k, alg, base);

        leak /= infoxtr::ksginfo::je(mat, {0}, k, alg, base);

        leak = std::max(0.0, std::min(1.0, leak));

        result.values.push_back(leak);
        result.types.push_back(3);
        result.var_indices.push_back(ag_idx);

        return result;
    }
    
} // namespace surd

}

#endif // INFOXTR_SURD_HPP
