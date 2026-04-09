#ifndef INFOXTR_KOCMI_HPP
#define INFOXTR_KOCMI_HPP

#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <thread>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infoxtr/numericutils.hpp"
#include "infoxtr/distance.hpp"
#include "infoxtr/infotheo.hpp"
#include <RcppThread.h>

namespace infoxtr 
{

namespace kocmi 
{
    using DiscVec = std::vector<uint64_t>;
    using DiscMat = std::vector<std::vector<uint64_t>>;
    using ContVec = std::vector<double>;
    using ContMat = std::vector<std::vector<double>>;

    /******************************************************************
     * Result structure
     ******************************************************************/
    struct KOCMIRes
    {
        double t_stat = std::numeric_limits<double>::quiet_NaN();
        double p_value = std::numeric_limits<double>::quiet_NaN();
    };

    /*****************************************************************
     * Utilities
     *****************************************************************/
    inline KOCMIRes permutation_test_mean(
        const std::vector<double>& diffs,
        size_t nboots = 10000,
        size_t threads = 1,
        uint64_t seed = 123456789) 
    {
        std::vector<double> vec;
        vec.reserve(diffs.size());
        
        double sum = 0.0;
        for (double v : diffs) 
        {
            if (!std::isnan(v)) 
            {
                vec.push_back(v);
                sum += v;
            }
        }

        const size_t n = vec.size();
        
        KOCMIRes result;
        // If all values are NaN or only have one non-NaN value
        if (n < 2) {
            return result;
        }

        if (threads == 0) threads = 1;
        size_t hw = std::thread::hardware_concurrency();
        if (hw > 0) threads = std::min(threads, hw);
        
        // Compute observed statistic
        const double observed_mean = sum / static_cast<double>(n);
        const double observed_stat = std::abs(observed_mean);

        double observed_sdval = 0.0;
        for (double v : vec) 
        {
            observed_sdval += (v - observed_mean) * (v - observed_mean);
        }
        observed_sdval = std::sqrt(observed_sdval / static_cast<double>(n - 1));
        
        if (infoxtr::numericutils::doubleNearlyEqual(observed_sdval, 0.0))
        {
            result.t_stat = observed_mean / observed_sdval;
        }

        /*
         * The single RNG version implementation is retained here for reference.
         * To facilitate parallelization of the permutation process without compromising randomness,
         * we adopted a strategy of pre-constructing an RNG pool to compute permutation statistics.
         * The pre-constructed RNG pool is initialized independently using the user-provided seed.
        */
        // std::mt19937_64 rng(seed);
        // std::uniform_int_distribution<int> sign_dist(0, 1);

        // size_t count = 0;

        // for (size_t b = 0; b < nboots; ++b) {

        //     double perm_sum = 0.0;

        //     for (size_t i = 0; i < n; ++i) {
        //         int sign = sign_dist(rng) ? 1 : -1;
        //         perm_sum += vec[i] * sign;
        //     }

        //     double perm_stat = std::abs(perm_sum / static_cast<double>(n));

        //     if (perm_stat >= observed_stat) {
        //         ++count;
        //     }
        // }
        
        // Prebuild 64-bit RNG pool for reproducibility
        std::vector<std::mt19937_64> rng_pool(nboots);
        for (size_t b = 0; b < nboots; ++b) 
        {
            std::seed_seq seq{static_cast<uint64_t>(seed), static_cast<uint64_t>(b)};
            rng_pool[b] = std::mt19937_64(seq);
        }

        std::uniform_int_distribution<int> sign_dist(0, 1);

        std::vector<size_t> perm_flags(nboots, 0);

        // Perform permutation
        if (threads <= 1)
        {
            for (size_t b = 0; b < nboots; ++b) 
            {
                double perm_sum = 0.0;
                std::mt19937_64& rng = rng_pool[b];

                for (size_t i = 0; i < n; ++i) 
                {
                    int sign = sign_dist(rng) ? 1 : -1;
                    perm_sum += vec[i] * sign;
                }

                double perm_stat = std::abs(perm_sum / static_cast<double>(n));

                perm_flags[b] = (perm_stat >= observed_stat) ? 1 : 0;
            }
        }
        else 
        {
            RcppThread::parallelFor(0, nboots, [&](size_t b) {

                double perm_sum = 0.0;
                std::mt19937_64& rng = rng_pool[b];

                for (size_t i = 0; i < n; ++i) 
                {
                    int sign = sign_dist(rng) ? 1 : -1;
                    perm_sum += vec[i] * sign;
                }

                double perm_stat = std::abs(perm_sum / static_cast<double>(n));

                perm_flags[b] = (perm_stat >= observed_stat) ? 1 : 0;

            }, threads);
        }
        
        // Compute p-value
        size_t perm_count = 0;
        for (size_t f : perm_flags) {
            perm_count += f;
        }
        
        result.p_value = static_cast<double>(perm_count) / static_cast<double>(nboots);

        return result;
    }

    inline double cmi(
        const ContVec& target,
        const ContVec& interact,
        const ContMat& conds,
        size_t k = 3,
        size_t alg = 0)
    {   
        ContMat xyz = conds;
        xyz.push_back(target);  
        xyz.push_back(interact);  

        ContMat xz = conds;
        xz.push_back(target);  

        std::vector<size_t> yz = conds;
        yz.push_back(interact);

        auto d_xyz = infoxtr::distance::distance(xyz,"maximum",true,false);
        auto d_xz  = infoxtr::distance::distance(xz,"maximum",true,false);
        auto d_yz  = infoxtr::distance::distance(yz,"maximum",true,false);
        auto d_z   = infoxtr::distance::distance(conds,"maximum",true,false);

        const size_t n = d_xyz.size();

        double sum = 0.0;

        for (size_t i = 0; i < n; ++i)
        {   
            auto& row = d_xyz[i];
            row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbor count");

            std::nth_element(row.begin(),row.begin()+k-1,row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            size_t nxz = 0, nyz = 0, nz = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i == j) continue;

                if (alg == 0)
                {
                    if (!std::isnan(d_xz[i][j]) && d_xz[i][j] < eps) nxz++;
                    if (!std::isnan(d_yz[i][j]) && d_yz[i][j] < eps) nyz++;
                    if (!std::isnan(d_z[i][j])  && d_z[i][j]  < eps) nz++;
                }
                else
                {
                    if (!std::isnan(d_xz[i][j]) && d_xz[i][j] <= eps) nxz++;
                    if (!std::isnan(d_yz[i][j]) && d_yz[i][j] <= eps) nyz++;
                    if (!std::isnan(d_z[i][j])  && d_z[i][j]  <= eps) nz++;    
                } 
            }

            if (alg == 0)
                sum += infoxtr::numericutils::digamma(nxz+1)
                     + infoxtr::numericutils::digamma(nyz+1)
                     - infoxtr::numericutils::digamma(nz+1);
            else
                sum += infoxtr::numericutils::digamma(nxz)
                     + infoxtr::numericutils::digamma(nyz)
                     - infoxtr::numericutils::digamma(nz);
        }

        double cmival = infoxtr::numericutils::digamma(k) - sum / n;
        if (alg == 1) cmival -= 1.0 / k;

        return cmival;
    }

    

    /*****************************************************************
     * Knockoff Conditional Mutual Information for Continuous Data
     *****************************************************************/
    inline KOCMIRes kocmi(
        const ContVec& target,
        const ContVec& agent,
        const ContMat& conds,
        const ContMat& all_knockoff,
        const ContMat& target_knockoff,
        size_t k = 3,
        size_t alg = 0,
        size_t threads = 1,
        double base = 2.0)
    {

    }

     /*****************************************************************
     * Knockoff Conditional Mutual Information for Discrete Data
     *****************************************************************/



} // namespace kocmi

}

#endif // INFOXTR_KOCMI_HPP
