

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
        
        // Construct joint state matrix
        DiscMat pm(tg.size()*2 + ag.size());

        for (size_t i = 0; i < tg.size(); ++i)
        {   
            pm[i] = mat[tg[i]];
        }
        for (size_t i = 0; i < ag.size(); ++i)
        {   
            pm[i + tg.size()] = mat[ag[i]];
        }

        
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
