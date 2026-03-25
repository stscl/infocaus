/**********************************************************************
 *  File: transferentropy.hpp
 *
 *  Transfer Entropy estimation for discrete and continuous time series.
 *
 *  Algorithms:
 *
 *      Transfer Entropy via Conditional Mutual Information
 *
 *      TE(X → Y) = I(Y_t ; X_{t-lag} | Y_{t-lag})
 *
 *  Data layout:
 *      DiscMat = std::vector<std::vector<uint64_t>>
 *                 // mat[var][obs]  (discrete series)
 *
 *      ContMat = std::vector<std::vector<double>>
 *                 // mat[var][obs]  (continuous observations)
 *
 *  Backend estimators:
 *
 *      Discrete Transfer Entropy
 *      -------------------------
 *      Uses discrete information theoretic estimators
 *      implemented in:
 *
 *          InfoTheo::CMI
 *
 *      Continuous Transfer Entropy
 *      ---------------------------
 *      Uses k-nearest neighbour estimators implemented in:
 *
 *          KSGInfo::CMI
 *
 *      Based on:
 *
 *          Kraskov–Stögbauer–Grassberger (KSG) estimator
 *
 *  Parameters:
 *
 *      lag_p, lag_q
 *          Time lag for target and agent variables.
 *
 *      k
 *          k-nearest neighbours for KSG estimator.
 *
 *      alg
 *          Estimator variant for continuous TE:
 *
 *              0  KSG estimator I (KSG1)
 *              1  KSG estimator II (KSG2)
 *
 *      normalize
 *          If true, normalize TE by conditional entropy.
 *
 *  Dependencies:
 *
 *      infotheo.hpp
 *      ksginfo.hpp
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 **********************************************************************/

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
        size_t lag_p = 3,
        size_t lag_q = 3,
        double base = 2.0,
        bool na_rm = true,
        bool normalize = false)
    {
        const size_t n_obs  = mat[0].size();
        const size_t n_cols = mat.size();

        if (mat.empty() || n_obs <= lag_p + 1 || n_obs <= lag_q + 1)
            return std::numeric_limits<double>::quiet_NaN();
        if (lag_p == 0 || lag_q == 0)
            return 0.0;

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
        size_t t0 = std::max(lag_p, lag_q);
        size_t N  = n_obs - t0;
        DiscMat pm(tg.size()*2 + ag.size(),
                   std::vector<uint64_t>(N, 0));
        
        // Y_t
        for (size_t i = 0; i < tg.size(); ++i)
        {   
            for (size_t t = t0; t < n_obs; ++t)
            {
                pm[i][t - t0] = mat[tg[i]][t];
            }
        }

        // X_{t-lag}
        for (size_t i = 0; i < ag.size(); ++i)
        {   
            for (size_t t = 0; t < N; ++t)
            {
                uint64_t v = mat[ag[i]][t - lag_q];
                if (v != 0)
                    pm[i + tg.size()][t - t0] = v;
            }
        }

        // Y_{t-lag}
        for (size_t i = 0; i < tg.size(); ++i)
        {   
            for (size_t t = 0; t < N; ++t)
            {
                uint64_t v = mat[tg[i]][t];
                if (v != 0)
                    pm[i + tg.size() + ag.size()][t] = v;
            }
        }

        // Construct variable index vector for CMI
        std::vector<size_t> tg_idx(tg.size());
        std::iota(tg_idx.begin(), tg_idx.end(), 0);

        std::vector<size_t> ag_idx(ag.size());
        std::iota(ag_idx.begin(), ag_idx.end(), tg.size());

        std::vector<size_t> tgl_idx(tg.size());
        std::iota(tgl_idx.begin(), tgl_idx.end(), tg.size() + ag.size());

        // Compute conditional mutual information
        return InfoTheo::CMI(pm, tg_idx, ag_idx, tgl_idx, base, na_rm, normalize);
    }

    /***********************************************************
     * Transfer Entropy for Continuous Data
     ***********************************************************/
    inline double TE4Cont(
        const ContMat& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& agent,
        size_t lag_p = 3,
        size_t lag_q = 3,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0,
        bool normalize = false)
    {
        const size_t n_obs  = mat[0].size();
        const size_t n_cols = mat.size();

        if (mat.empty() || n_obs <= lag_p + k + 1 || n_obs <= lag_q + k + 1)
            return std::numeric_limits<double>::quiet_NaN();
        if (lag_p == 0 || lag_q == 0)
            return 0.0;

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
        size_t t0 = std::max(lag_p, lag_q);
        size_t N  = n_obs - t0;
        ContMat pm(tg.size()*2 + ag.size(),
                std::vector<double>(N,
                std::numeric_limits<double>::quiet_NaN()));
        
        // Y_t
        for (size_t i = 0; i < tg.size(); ++i)
        {   
            for (size_t t = t0; t < n_obs; ++t)
            {
                pm[i][t - t0] = mat[tg[i]][t];
            }
        }

        // X_{t-lag}
        for (size_t i = 0; i < ag.size(); ++i)
        {   
            for (size_t t = t0; t < n_obs; ++t)
            {
                double v = mat[ag[i]][t - lag_q];
                if (!std::isnan(v))
                    pm[i + tg.size()][t - t0] = v;
            }
        }

        // Y_{t-lag}
        for (size_t i = 0; i < tg.size(); ++i)
        {   
            for (size_t t = t0; t < n_obs; ++t)
            {
                double v = mat[tg[i]][t - lag_p];
                if (!std::isnan(v))
                    pm[i + tg.size() + ag.size()][t - t0] = v;
            }
        }

        // Construct variable index vector for CMI
        std::vector<size_t> tg_idx(tg.size());
        std::iota(tg_idx.begin(), tg_idx.end(), 0);

        std::vector<size_t> ag_idx(ag.size());
        std::iota(ag_idx.begin(), ag_idx.end(), tg.size());

        std::vector<size_t> tgl_idx(tg.size());
        std::iota(tgl_idx.begin(), tgl_idx.end(), tg.size() + ag.size());

        // Compute conditional mutual information
        return KSGInfo::CMI(pm, tg_idx, ag_idx, tgl_idx, k, alg, base, normalize);
    }    

} // namespace TE

#endif // TRANSFERENTROPY_HPP
