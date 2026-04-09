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
    inline double permutation_test_mean(const std::vector<double>& diffs,
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

        // If all values are NaN
        if (vec.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const size_t n = vec.size();

        // Compute observed statistic
        const double observed_stat = std::abs(sum / static_cast<double>(n));

        /*
         * The single RNG version implementation is retained here for reference.
         * To facilitate parallelization of the permutation process without compromising randomness,
         * we adopted a strategy of pre-constructing an RNG pool to compute permutation statistics.
         * The pre-constructed RNG pool is initialized independently using the user-provided seed.
        */
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<int> sign_dist(0, 1);

        size_t count = 0;

        for (size_t b = 0; b < nboots; ++b) {

            double perm_sum = 0.0;

            for (size_t i = 0; i < n; ++i) {
                int sign = sign_dist(rng) ? 1 : -1;
                perm_sum += vec[i] * sign;
            }

            double perm_stat = std::abs(perm_sum / static_cast<double>(n));

            if (perm_stat >= observed_stat) {
                ++count;
            }
        }

        return static_cast<double>(count) / static_cast<double>(nboots);
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
