/**********************************************************************
 *  File: transferentropy.hpp
 *
 *  Transfer Entropy estimation for discrete and continuous time series.
 *
 *  Algorithms:
 *
 *      Transfer Entropy via Conditional Mutual Information
 *
 *      TE(X → Y) = I(Y_t ; X_{t-q}^{(q)} | Y_{t-p}^{(p)})
 *
 *      where
 *
 *          X_{t-q}^{(q)} = (X_{t-1}, X_{t-2}, ..., X_{t-q})
 *          Y_{t-p}^{(p)} = (Y_{t-1}, Y_{t-2}, ..., Y_{t-p})
 *
 *  Data layout:
 *
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
 *      lag_p
 *          Target history length.
 *
 *      lag_q
 *          Agent history length.
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
 *      lag_single
 *
 *          true  → use only single lag
 *
 *              X_{t-lag_q}
 *              Y_{t-lag_p}
 *
 *          false → use full lag embedding
 *
 *              (X_{t-1} ... X_{t-lag_q})
 *              (Y_{t-1} ... Y_{t-lag_p})
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
        bool normalize = false,
        bool lag_single = false)
    {   
        const size_t n_obs = mat[0].size();
        const size_t n_cols = mat.size();

        if (mat.empty())
            return std::numeric_limits<double>::quiet_NaN();

        size_t t0 = std::max(lag_p, lag_q);
        if (n_obs <= t0 + 1)
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
        size_t tg_lag = lag_single ? tg.size() : tg.size() * lag_p;
        size_t ag_lag = lag_single ? ag.size() : ag.size() * lag_q;
        size_t N = n_obs - t0;
        DiscMat pm(tg.size() + ag_lag + tg_lag,
                   std::vector<uint64_t>(N,0));
        
        // Y_present
        for (size_t i = 0; i < tg.size(); ++i)
        {   
            for (size_t t = t0; t < n_obs; ++t)
            {
                pm[i][t - t0] = mat[tg[i]][t];
            }
        }

        size_t col = tg.size();

        // X_past
        for (size_t i = 0; i < ag.size(); ++i)
        {
            if (lag_single)
            {
                for (size_t t = t0; t < n_obs; ++t)
                    pm[col + i][t - t0] = mat[ag[i]][t - lag_q];
            }
            else
            {
                for (size_t l = 1; l <= lag_q; ++l)
                    for (size_t t = t0; t < n_obs; ++t)
                        pm[col + i * lag_q + (l - 1)][t - t0] = mat[ag[i]][t - l];
            }
        }

        col += ag_lag;

        // Y_past
        for (size_t i = 0; i < tg.size(); ++i)
        {
            if (lag_single)
            {
                for (size_t t = t0; t < n_obs; ++t)
                    pm[col + i][t - t0] = mat[tg[i]][t - lag_p];
            }
            else
            {
                for (size_t l = 1; l <= lag_p; ++l)
                    for (size_t t = t0; t < n_obs; ++t)
                        pm[col + i * lag_p + (l - 1)][t - t0] = mat[tg[i]][t - l];
            }
        }

        // Construct variable indices vector for CMI
        std::vector<size_t> tg_idx(tg.size());
        std::iota(tg_idx.begin(), tg_idx.end(), 0);

        std::vector<size_t> ag_idx(ag_lag);
        std::iota(ag_idx.begin(), ag_idx.end(), tg.size());

        std::vector<size_t> tgl_idx(tg_lag);
        std::iota(tgl_idx.begin(), tgl_idx.end(), tg.size() + ag_lag);

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
        bool normalize = false,
        bool lag_single = false)
    {   
        const size_t n_obs = mat[0].size();
        const size_t n_cols = mat.size();

        if (mat.empty())
            return std::numeric_limits<double>::quiet_NaN();

        size_t t0 = std::max(lag_p, lag_q);
        if (n_obs <= t0 + k + 1)
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
        size_t tg_lag = lag_single ? tg.size() : tg.size() * lag_p;
        size_t ag_lag = lag_single ? ag.size() : ag.size() * lag_q;
        size_t N = n_obs - t0;
        ContMat pm(tg.size() + ag_lag + tg_lag,
                   std::vector<double>(N, std::numeric_limits<double>::quiet_NaN()));
        
        // Y_present
        for (size_t i = 0; i < tg.size(); ++i)
        {   
            for (size_t t = t0; t < n_obs; ++t)
            {
                pm[i][t - t0] = mat[tg[i]][t];
            }
        }

        size_t col = tg.size();

        // X_past
        for (size_t i = 0; i < ag.size(); ++i)
        {
            if (lag_single)
            {
                for (size_t t = t0; t < n_obs; ++t)
                    pm[col + i][t - t0] = mat[ag[i]][t - lag_q];
            }
            else
            {
                for (size_t l = 1; l <= lag_q; ++l)
                    for (size_t t = t0; t < n_obs; ++t)
                        pm[col + i * lag_q + (l - 1)][t - t0] = mat[ag[i]][t - l];
            }
        }

        col += ag_lag;

        // Y_past
        for (size_t i = 0; i < tg.size(); ++i)
        {
            if (lag_single)
            {
                for (size_t t = t0; t < n_obs; ++t)
                    pm[col + i][t - t0] = mat[tg[i]][t - lag_p];
            }
            else
            {
                for (size_t l = 1; l <= lag_p; ++l)
                    for (size_t t = t0; t < n_obs; ++t)
                        pm[col + i * lag_p + (l - 1)][t - t0] = mat[tg[i]][t - l];
            }
        }

        // Construct variable indices vector for CMI
        std::vector<size_t> tg_idx(tg.size());
        std::iota(tg_idx.begin(), tg_idx.end(), 0);

        std::vector<size_t> ag_idx(ag_lag);
        std::iota(ag_idx.begin(), ag_idx.end(), tg.size());

        std::vector<size_t> tgl_idx(tg_lag);
        std::iota(tgl_idx.begin(), tgl_idx.end(), tg.size() + ag_lag);

        // Compute conditional mutual information
        return KSGInfo::CMI(pm, tg_idx, ag_idx, tgl_idx, k, alg, base, normalize);
    }    

} // namespace TE

#endif // TRANSFERENTROPY_HPP
