/******************************************************************************
 * File: surd.hpp
 *
 * Synergistic–Unique–Redundant Decomposition (SURD)
 * -------------------------------------------------
 *
 * This module implements the SURD framework for decomposing the mutual
 * information between a target variable and multiple source variables into
 * interpretable components:
 *
 *   Redundant information
 *       Information simultaneously provided by multiple sources.
 *
 *   Unique information
 *       Information provided exclusively by a single source.
 *
 *   Synergistic information
 *       Information that only emerges when variables are considered jointly.
 *
 *   Information leak
 *       Remaining uncertainty in the target after conditioning on all sources.
 *
 * ---------------------------------------------------------------------------
 * Algorithm overview
 * ---------------------------------------------------------------------------
 *
 * Let Y denote the target variable and X = {X1, X2, ..., Xn} the source
 * variables. The algorithm proceeds in the following stages:
 *
 * 1. Subset enumeration
 *
 *      Generate all subsets of source variables up to order max_order.
 *      Each subset represents a candidate information channel.
 *
 * 2. Joint distribution construction
 *
 *      A joint state table is constructed for the full variable set using
 *      the joint entropy utilities in infotheo.hpp.
 *
 * 3. Conditional pointwise mutual information
 *
 *      For each target state s, compute pointwise mutual information
 *
 *          I_s(X_set ; Y)
 *
 *      for every subset X_set using grouped projections of the joint
 *      state table.
 *
 *      Mutual information is accumulated as:
 *
 *          I(X_set ; Y) = sum_s p(s) * I_s(X_set ; Y)
 *
 * 4. Monotonic SURD filtering
 *
 *      Subsets are sorted by their pointwise mutual information values.
 *      A monotonic constraint is enforced so that higher-order subsets
 *      cannot contain less information than the maximum of lower-order
 *      subsets. Violations are clipped to zero.
 *
 * 5. Information layer decomposition
 *
 *      A ladder-style decomposition is applied to the sorted values.
 *      Incremental differences between successive layers determine
 *      the information contributions.
 *
 *          delta_i = max(I_i - I_{i-1}, 0)
 *
 *      Contributions are classified as:
 *
 *          |subset| = 1   → redundant / unique layer
 *          |subset| > 1   → synergistic layer
 *
 * 6. Aggregation across target states
 *
 *      Contributions are weighted by p(s) and accumulated across
 *      target states.
 *
 * 7. Information leak
 *
 *      Remaining uncertainty in the target is measured as
 *
 *          H(Y | X_all) / H(Y)
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
#include <RcppThread.h>

namespace infoxtr 
{

namespace surd 
{

    using DiscMat = std::vector<std::vector<uint64_t>>;

    /***********************************************************
     * Result structure
     ***********************************************************/
    struct SURDRes
    {
        std::vector<std::vector<size_t>> unique_vars;
        std::vector<double>              unique_vals;

        std::vector<std::vector<size_t>> redundant_vars;
        std::vector<double>              redundant_vals;

        std::vector<std::vector<size_t>> synergy_vars;
        std::vector<double>              synergy_vals;

        std::vector<std::vector<size_t>> mi_vars;
        std::vector<double>              mi_vals;

        double info_leak = 0.0;
    };

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
            throw std::invalid_argument("SURD requires >=2 variables");

        const size_t n_vars    = mat.size();
        const size_t n_sources = n_vars - 1;

        max_order = std::min(max_order , n_sources);

        const double log_base = std::log(base);

        /***********************************************************
         * Construct variable combination vector
         ***********************************************************/
        std::vector<size_t> ag_idx(n_sources);
        std::iota(ag_idx.begin(), ag_idx.end(), 1);

        const std::vector<std::vector<size_t>> combs =
            infoxtr::combn::genSubsets(ag_idx, max_order);

        const size_t n_combs = combs.size();
        
        // Initialize containers for result
        std::vector<double> info(n_combs , 0.0);
        std::vector<double> I_unique(n_sources , 0.0);
        std::vector<std::vector<size_t>> unique_vars(n_sources);

        for (size_t i = 0; i < n_sources; i++)
            unique_vars[i] = {i + 1};
        
        std::vector<double> I_R(n_combs , 0.0);
        std::vector<double> I_S(n_combs , 0.0);

        /***********************************************************
         * Joint table
         ***********************************************************/
        infoxtr::infotheo::JointTable jt =
            infoxtr::infotheo::joint_table(mat);

        const size_t k        = jt.k;
        const size_t n_states = jt.counts.size();

        size_t total_n = 0;
        for (auto c : jt.counts)
            total_n += c;

        /***********************************************************
         * Target states
         ***********************************************************/
        std::vector<uint64_t> s_vals;
        std::vector<size_t>   s_counts;

        for (size_t i = 0; i < n_states; i++)
        {
            uint64_t s = jt.states[i * k];

            auto it = std::find(s_vals.begin(), s_vals.end(), s);

            if (it == s_vals.end())
            {
                s_vals.push_back(s);
                s_counts.push_back(jt.counts[i]);
            }
            else
            {
                size_t idx = it - s_vals.begin();
                s_counts[idx] += jt.counts[i];
            }
        }

        const size_t n_s = s_vals.size();

        /***********************************************************
         * state -> s index
         ***********************************************************/
        std::vector<size_t> state_s_index(n_states);

        for (size_t i = 0; i < n_states; i++)
        {
            uint64_t s = jt.states[i * k];

            for (size_t j = 0; j < n_s; j++)
                if (s_vals[j] == s)
                    state_s_index[i] = j;
        }

        /***********************************************************
         * Precompute projection + px groups
         ***********************************************************/
        struct PxGroup
        {
            std::vector<uint64_t> state;
            size_t count = 0;
            std::vector<size_t> rows;
        };

        std::vector<std::vector<PxGroup>> px_groups(n_combs);

        if (threads <= 1) 
        {   
            for (size_t ci = 0; ci < n_combs; ci++)
            {
                auto & subset = combs[ci];
                const size_t kk = subset.size();

                std::vector<std::vector<uint64_t>> proj(n_states);

                for (size_t i = 0; i < n_states; i++)
                {
                    size_t base = i * k;

                    proj[i].resize(kk);

                    for (size_t j = 0; j < kk; j++)
                        proj[i][j] = jt.states[base + subset[j]];
                }

                std::vector<size_t> order(n_states);

                for (size_t i = 0; i < n_states; i++)
                    order[i] = i;

                std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b)
                {
                    for (size_t j = 0; j < kk; j++)
                    {
                        if (proj[a][j] < proj[b][j]) return true;
                        if (proj[a][j] > proj[b][j]) return false;
                    }
                    return false;
                });

                PxGroup g;
                g.state = proj[order[0]];
                g.count = jt.counts[order[0]];
                g.rows.push_back(order[0]);

                for (size_t i = 1; i < n_states; i++)
                {
                    bool same = true;

                    for (size_t j = 0; j < kk; j++)
                        if (proj[order[i]][j] != g.state[j])
                            same = false;

                    if (same)
                    {
                        g.count += jt.counts[order[i]];
                        g.rows.push_back(order[i]);
                    }
                    else
                    {
                        px_groups[ci].push_back(g);

                        g.state = proj[order[i]];
                        g.count = jt.counts[order[i]];
                        g.rows.clear();
                        g.rows.push_back(order[i]);
                    }
                }

                px_groups[ci].push_back(g);
            }
        } 
        else  
        {
            RcppThread::parallelFor(0, n_combs, [&](size_t ci) {
                auto & subset = combs[ci];
                const size_t kk = subset.size();

                std::vector<std::vector<uint64_t>> proj(n_states);

                for (size_t i = 0; i < n_states; i++)
                {
                    size_t base = i * k;

                    proj[i].resize(kk);

                    for (size_t j = 0; j < kk; j++)
                        proj[i][j] = jt.states[base + subset[j]];
                }

                std::vector<size_t> order(n_states);

                for (size_t i = 0; i < n_states; i++)
                    order[i] = i;

                std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b)
                {
                    for (size_t j = 0; j < kk; j++)
                    {
                        if (proj[a][j] < proj[b][j]) return true;
                        if (proj[a][j] > proj[b][j]) return false;
                    }
                    return false;
                });

                PxGroup g;
                g.state = proj[order[0]];
                g.count = jt.counts[order[0]];
                g.rows.push_back(order[0]);

                for (size_t i = 1; i < n_states; i++)
                {
                    bool same = true;

                    for (size_t j = 0; j < kk; j++)
                        if (proj[order[i]][j] != g.state[j])
                            same = false;

                    if (same)
                    {
                        g.count += jt.counts[order[i]];
                        g.rows.push_back(order[i]);
                    }
                    else
                    {
                        px_groups[ci].push_back(g);

                        g.state = proj[order[i]];
                        g.count = jt.counts[order[i]];
                        g.rows.clear();
                        g.rows.push_back(order[i]);
                    }
                }

                px_groups[ci].push_back(g);
            }, threads);
        } 

        /***********************************************************
         * Loop target states
         ***********************************************************/
        for (size_t si = 0; si < n_s; si++)
        {
            double p_s =
                static_cast<double>(s_counts[si]) / total_n;

            std::vector<double> I_s(n_combs , 0.0);
            
            if (threads <= 1)
            {
                for (size_t ci = 0; ci < n_combs; ci++)
                {
                    double sum = 0.0;

                    for (auto & g : px_groups[ci])
                    {
                        size_t psx_count = 0;

                        for (auto r : g.rows)
                            if (state_s_index[r] == si)
                                psx_count += jt.counts[r];

                        if (psx_count == 0)
                            continue;

                        double psx =
                            static_cast<double>(psx_count) / total_n;

                        double px =
                            static_cast<double>(g.count) / total_n;

                        sum += psx * std::log(psx / (p_s * px));
                    }

                    double pointwise = sum / p_s / log_base;

                    I_s[ci] = pointwise;

                    // accumulate MI 
                    info[ci] += p_s * pointwise;
                }
            }
            else 
            {
                RcppThread::parallelFor(0, n_combs, [&](size_t ci) {
                    double sum = 0.0;

                    for (auto & g : px_groups[ci])
                    {
                        size_t psx_count = 0;

                        for (auto r : g.rows)
                            if (state_s_index[r] == si)
                                psx_count += jt.counts[r];

                        if (psx_count == 0)
                            continue;

                        double psx =
                            static_cast<double>(psx_count) / total_n;

                        double px =
                            static_cast<double>(g.count) / total_n;

                        sum += psx * std::log(psx / (p_s * px));
                    }

                    double pointwise = sum / p_s / log_base;

                    I_s[ci] = pointwise;

                    // accumulate MI 
                    info[ci] += p_s * pointwise;
                }, threads);
            }

            /**************************************************
             * SURD decomposition
             **************************************************/

            struct Node
            {
                size_t idx;
                size_t len;
                double val;
            };

            std::vector<Node> nodes;
            nodes.reserve(n_combs);

            for (size_t i = 0; i < n_combs; i++)
            {
                Node n;
                n.idx = i;
                n.len = combs[i].size();
                n.val = I_s[i];
                nodes.push_back(n);
            }

            // sort once 
            std::sort(nodes.begin(), nodes.end(),
            [](const Node & a, const Node & b)
            {
                if (!infoxtr::numericutils::doubleNearlyEqual(a.val, b.val))
                        return a.val < b.val;

                    return a.idx < b.idx;
            });

            // find max subset length 
            size_t max_len = 0;

            for (auto & n : nodes)
                if (n.len > max_len)
                    max_len = n.len;

            /**************************************************
             * SURD monotonic filter
             **************************************************/

            for (size_t l = 1; l < max_len; l++)
            {
                double Il1max = -std::numeric_limits<double>::infinity();

                for (auto & n : nodes)
                    if (n.len == l)
                        if (n.val > Il1max)
                            Il1max = n.val;

                for (auto & n : nodes)
                    if (n.len == l + 1)
                        if (n.val < Il1max)
                            n.val = 0.0;
            }

            /**************************************************
             * SURD information layers
             **************************************************/

            std::vector<size_t> red_vars(n_sources);
            for(size_t i = 0; i < n_sources; i++)
                red_vars[i] = i+1;

            double prev = 0.0;

            for (auto & n : nodes)
            {
                double delta = n.val - prev;

                if (delta < 0)
                    delta = 0;

                double info_add = delta * p_s;

                const auto & subset = combs[n.idx];

                if (subset.size() == 1)
                {
                    // redundant information
                    for (size_t ri = 0; ri < n_combs; ri++)
                    {
                        if (combs[ri] == red_vars)
                        {
                            I_R[ri] += info_add;
                            break;
                        }
                    }

                    auto it = std::find(
                        red_vars.begin(),
                        red_vars.end(),
                        subset[0]);

                    // NOTE (GCC13/CRAN):
                    // Avoid std::vector::erase() here. GCC13 may emit
                    // -Wstringop-overflow warnings due to internal
                    // memmove static analysis. We remove elements
                    // using swap + pop_back instead, which avoids
                    // the memmove path in STL.

                    // if (it != red_vars.end())
                    //     red_vars.erase(it);

                    if (it != red_vars.end())
                    {
                        std::swap(*it, red_vars.back());
                        red_vars.pop_back();
                    }
                }
                else
                {
                    // synergy information
                    I_S[n.idx] += info_add;
                }

                if (n.val > prev)
                    prev = n.val;
            }
        }

        SURDRes result;

        result.unique_vars.reserve(n_sources);
        result.redundant_vars.reserve(n_combs);
        result.synergy_vars.reserve(n_combs);
        result.mi_vars.reserve(n_combs);

        /*
        NOTE

        We use `insert(end(), ...)` instead of push_back/emplace_back
        for nested containers (vector<vector<size_t>>).

        Some GCC versions (11–14) may emit false-positive warnings
        such as -Wstringop-overflow or -Warray-bounds when copying
        nested vectors via push_back/emplace_back due to STL
        memmove analysis.

        Using insert(end(), ...) avoids those static analyzer paths
        and keeps the code portable across CRAN build systems.
        */

        for (size_t i = 0; i < n_combs; ++i)
        {
            const auto & c = combs[i];

            if (I_R[i] > 0)
            {
                if (c.size() == 1)
                {
                    result.unique_vars.insert(result.unique_vars.end(), c);
                    result.unique_vals.emplace_back(I_R[i]);
                }
                else
                {
                    result.redundant_vars.insert(result.redundant_vars.end(), c);
                    result.redundant_vals.emplace_back(I_R[i]);
                }
            }

            if (c.size() > 1 && I_S[i] > 0)
            {
                result.synergy_vars.insert(result.synergy_vars.end(), c);
                result.synergy_vals.emplace_back(I_S[i]);
            }

            result.mi_vars.insert(result.mi_vars.end(), c);
            result.mi_vals.emplace_back(info[i]);
        }

        // Information leak
        std::vector<size_t> target_idx = {0};
        double leak = infoxtr::infotheo::ce(mat, target_idx, ag_idx, base, false) 
                    / infoxtr::infotheo::je(mat, target_idx, base, false);
        result.info_leak = std::max(0.0, std::min(1.0, leak));

        /***********************************************************
         * Normalize (optional)
         ***********************************************************/
        if (normalize)
        {
            double max_mi = 0.0;

            for (double v : result.mi_vals)
                if (v > max_mi)
                    max_mi = v;

            if (max_mi <= 0.0)
                max_mi = 1.0;

            for (auto & v : result.unique_vals)
                v /= max_mi;

            for (auto & v : result.redundant_vals)
                v /= max_mi;

            for (auto & v : result.synergy_vals)
                v /= max_mi;

            // for (auto & v : result.mi_vals)
            //     v /= max_mi;
        }

        return result;
    }
    
} // namespace surd

}

#endif // INFOXTR_SURD_HPP
