/**********************************************************************
 *  File: infotheo.hpp
 *
 *  High performance information theoretic measurements
 *  for discrete state data (uint64_t encoded).
 *
 *  Data layout:
 *      Series = std::vector<uint64_t>
 *      Matrix = std::vector<std::vector<uint64_t>>   // mat[var][obs]
 *
 *  Missing value convention:
 *      Value 0 is treated as NA when na_rm=true
 *      Ensure your discrete encoding reserves 0 for missing data
 *
 *  Functions:
 *      entropy  Shannon Entropy
 *      je       Joint Entropy
 *      ce       Conditional Entropy
 *      mi       Mutual Information
 *      cmi      Conditional Mutual Information
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 **********************************************************************/

#ifndef INFOXTR_INFOTHEO_HPP
#define INFOXTR_INFOTHEO_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace infoxtr
{

namespace infotheo
{

    using Series = std::vector<uint64_t>;
    using Matrix = std::vector<std::vector<uint64_t>>;

    /***********************************************************
     * Utilities
     ***********************************************************/
    inline double convert_log_base(double x, double base)
    {
        if (x <= 0.0) return 0.0;

        if (!(base > 0.0) || std::abs(base - 1.0) < 1e-12)
            throw std::invalid_argument("Log base must be positive and not equal to 1.");

        return x / std::log(base);
    }

    /***********************************************************
     * Entropy
     ***********************************************************/
    inline double entropy(
        const Series& series,
        double base = 2.0,
        bool na_rm = true)
    {
        if (series.empty())
            return std::numeric_limits<double>::quiet_NaN();

        std::unordered_map<uint64_t, size_t> freq;
        freq.reserve(series.size() * 1.3);

        size_t n_valid = 0;

        for (uint64_t x : series)
        {
            if (na_rm && x == 0)
                continue;

            ++freq[x];
            ++n_valid;
        }

        if (n_valid == 0)
            return std::numeric_limits<double>::quiet_NaN();

        double H = 0.0;
        for (const auto& kv : freq)
        {
            double p = static_cast<double>(kv.second) / n_valid;
            H -= p * std::log(p);
        }

        return convert_log_base(H, base);
    }

    /***********************************************************
     * Joint Entropy
     ***********************************************************/
    inline double je(
        const Matrix& mat,
        const std::vector<size_t>& vars,
        double base = 2.0,
        bool na_rm = true)
    {
        if (mat.empty() || vars.empty())
            return std::numeric_limits<double>::quiet_NaN();

        const size_t n_obs  = mat[0].size();
        const size_t n_cols = mat.size();

        /*------------------------------------------------------
         * Validate, sort, and deduplicate variable indices
         *-----------------------------------------------------*/
        std::vector<size_t> clean_vars;
        clean_vars.reserve(vars.size());

        for (size_t v : vars)
            if (v < n_cols)
                clean_vars.push_back(v);

        // Sort + unique: remove duplicates and ensure deterministic order
        std::sort(clean_vars.begin(), clean_vars.end());
        clean_vars.erase(
            std::unique(clean_vars.begin(), clean_vars.end()),
            clean_vars.end()
        );

        if (clean_vars.empty())
            return std::numeric_limits<double>::quiet_NaN();

        const size_t k = clean_vars.size();

        /*------------------------------------------------------
         * Construct joint state table
         * Layout:
         *   states[i*k + j]
         *-----------------------------------------------------*/
        std::vector<uint64_t> states;
        states.reserve(n_obs * k);

        size_t n_valid = 0;
        for (size_t i = 0; i < n_obs; ++i)
        {
            // Phase 1: Check for NA first (read-only)
            bool skip = false;
            for (size_t j = 0; j < k; ++j)
            {
                if (na_rm && mat[clean_vars[j]][i] == 0)
                {
                    skip = true;
                    break;
                }
            }
            if (skip) continue;
            
            // Phase 2: Push valid values (guaranteed complete)
            for (size_t j = 0; j < k; ++j)
                states.push_back(mat[clean_vars[j]][i]);
            
            ++n_valid;
        }

        if (n_valid == 0)
            return std::numeric_limits<double>::quiet_NaN();

        /*------------------------------------------------------
         * Create index array for sorting
         *-----------------------------------------------------*/
        std::vector<size_t> order(n_valid);
        for (size_t i = 0; i < n_valid; ++i)
            order[i] = i;

        /*------------------------------------------------------
         * Sort joint states lexicographically
         *-----------------------------------------------------*/
        std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b)
            {
                size_t ia = a * k;
                size_t ib = b * k;

                for (size_t j = 0; j < k; ++j)
                {
                    uint64_t va = states[ia + j];
                    uint64_t vb = states[ib + j];

                    if (va < vb) return true;
                    if (va > vb) return false;
                }
                return false;
            });

        /*------------------------------------------------------
         * Run-length frequency counting
         *-----------------------------------------------------*/
        double H = 0.0;
        size_t run = 1;

        auto equal_state = [&](size_t a, size_t b)
        {
            size_t ia = a * k;
            size_t ib = b * k;

            for (size_t j = 0; j < k; ++j)
                if (states[ia + j] != states[ib + j])
                    return false;

            return true;
        };

        for (size_t i = 1; i < n_valid; ++i)
        {
            if (equal_state(order[i], order[i-1]))
            {
                ++run;
            }
            else
            {
                double p = static_cast<double>(run) / n_valid;
                H -= p * std::log(p);
                run = 1;
            }
        }

        double p = static_cast<double>(run) / n_valid;
        H -= p * std::log(p);

        return convert_log_base(H, base);
    }

    /***********************************************************
     * Conditional Entropy
     ***********************************************************/
    inline double ce(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& conds,
        double base = 2.0,
        bool na_rm = true)
    {
        if (mat.empty() || target.empty())
            return std::numeric_limits<double>::quiet_NaN();

        if (conds.empty())
            throw std::invalid_argument("conds cannot be empty");

        std::vector<size_t> tc = conds;
        tc.insert(tc.end(), target.begin(), target.end());

        return je(mat, tc, base, na_rm) - je(mat, conds, base, na_rm);
    }

    /***********************************************************
     * Mutual Information
     ***********************************************************/
    inline double mi(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& interact,
        double base = 2.0,
        bool na_rm = true,
        bool normalize = false)
    {
        if (mat.empty() || target.empty())
            return std::numeric_limits<double>::quiet_NaN();

        if (interact.empty())
            throw std::invalid_argument("interact cannot be empty");

        std::vector<size_t> ti = interact;
        ti.insert(ti.end(), target.begin(), target.end());

        double ht = je(mat, target, base, na_rm); 
        double hi = je(mat, interact, base, na_rm); 
        double hti = je(mat, ti, base, na_rm); 
        
        double mival = ht + hi - hti; 
        
        if (!normalize) return mival; 
        if (hti <= 0) return mival; 

        return mival / hti;
    }

    /***********************************************************
     * Conditional Mutual Information
     ***********************************************************/
    inline double cmi(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& interact,
        const std::vector<size_t>& conds,
        double base = 2.0,
        bool na_rm = true,
        bool normalize = false)
    {
        if (mat.empty() || target.empty())
            return std::numeric_limits<double>::quiet_NaN();

        if (interact.empty() || conds.empty())
            throw std::invalid_argument("interact / conds cannot be empty");

        std::vector<size_t> ct  = conds;
        std::vector<size_t> ci  = conds;
        std::vector<size_t> cti = conds;

        ct.insert(ct.end(), target.begin(), target.end());
        ci.insert(ci.end(), interact.begin(), interact.end());
        cti.insert(cti.end(), target.begin(), target.end());
        cti.insert(cti.end(), interact.begin(), interact.end());

        double h_ct  = je(mat, ct, base, na_rm);
        double h_ci  = je(mat, ci, base, na_rm);
        double h_c   = je(mat, conds, base, na_rm);
        double h_cti = je(mat, cti, base, na_rm);

        double cmival = h_ct + h_ci - h_c - h_cti;

        if (!normalize) return cmival;

        std::vector<size_t> ti = interact;
        ti.insert(ti.end(), target.begin(), target.end());

        double h_ti_c = ce(mat, ti, conds, base, na_rm);
        if (h_ti_c <= 0) return cmival;

        return cmival / h_ti_c;
    }

} // namespace infotheo

}

#endif // INFOXTR_INFOTHEO_HPP
