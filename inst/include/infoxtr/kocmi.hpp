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
            return KOCMIRes;
        }
        
        // Compute observed statistic
        const double observed_mean = sum / static_cast<double>(n)
        const double observed_stat = std::abs(observed_mean);

        double observed_sdval = 0.0;
        for (double v : vec) 
        {
            observed_sdval += (v - observed_mean) * (v - observed_mean);
        }
        result.t_stat = observed_mean / std::sqrt(observed_sdval / static_cast<double>(n - 1));

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

        // Compute p-value
        size_t perm_count = 0;
        for (size_t f : perm_flags) {
            perm_sum += f;
        }
    }

    

    /*****************************************************************
     * Knockoff Conditional Mutual Information for Continuous Data
     *****************************************************************/

     /*****************************************************************
     * Knockoff Conditional Mutual Information for Discrete Data
     *****************************************************************/



} // namespace kocmi

}

#endif // INFOXTR_KOCMI_HPP
