

#ifndef TRANSFERENTROPY_HPP
#define TRANSFERENTROPY_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infotheo.hpp"
#include "ksginfo.hpp"

namespace TE
{   
    using DiscMat = std::vector<std::vector<uint64_t>>;
    using ContMat = std::vector<std::vector<double>>;

    /***********************************************************
     * Transfer Entropy for Discrete Data
     ***********************************************************/
    inline double TE4Disc(
        const DiscMat& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& agent,
        size_t lag = 3,
        double base = 2.0,
        bool na_rm = true,
        bool normalize = false)
    {
        if (mat.empty())
            return std::numeric_limits<double>::quiet_NaN();

        const size_t n_obs  = mat[0].size();
        const size_t n_cols = mat.size();

        // Validate, sort, and deduplicate variable indices
        std::vector<size_t> tg;
        tg.reserve(target.size());
        for (size_t v : target)
            if (v < n_cols)
                tg.push_back(v);
        std::sort(tg.begin(), tg.end());
        tg.erase(
            std::unique(tg.begin(), tg.end()),
            tg.end()
        );

        std::vector<size_t> ag;
        ag.reserve(agent.size());
        for (size_t v : agent)
            if (v < n_cols)
                ag.push_back(v);
        std::sort(ag.begin(), ag.end());
        ag.erase(
            std::unique(ag.begin(), ag.end()),
            ag.end()
        );

        if (tg.empty() || ag.empty())
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
    inline double CE(
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

        return JE(mat, tc, base, na_rm)
            - JE(mat, conds, base, na_rm);
    }

    /***********************************************************
     * Mutual Information
     ***********************************************************/
    inline double MI(
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

        double ht = JE(mat, target, base, na_rm); 
        double hi = JE(mat, interact, base, na_rm); 
        double hti = JE(mat, ti, base, na_rm); 
        
        double mi = ht + hi - hti; 
        
        if (!normalize) return mi; 
        if (hti <= 0) return mi; 

        return mi / hti;
    }

    /***********************************************************
     * Conditional Mutual Information
     ***********************************************************/
    inline double CMI(
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

        double h_ct  = JE(mat, ct, base, na_rm);
        double h_ci  = JE(mat, ci, base, na_rm);
        double h_c   = JE(mat, conds, base, na_rm);
        double h_cti = JE(mat, cti, base, na_rm);

        double cmi = h_ct + h_ci - h_c - h_cti;

        if (!normalize) return cmi;

        std::vector<size_t> ti = interact;
        ti.insert(ti.end(), target.begin(), target.end());

        double h_ti_c = CE(mat, ti, conds, base, na_rm);
        if (h_ti_c <= 0) return cmi;

        return cmi / h_ti_c;
    }

} // namespace InfoTheo

#endif // TRANSFERENTROPY_HPP
