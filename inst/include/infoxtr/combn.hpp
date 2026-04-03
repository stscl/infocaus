/**************************************************************************
 * File: combn.hpp
 *
 * Combinatorial utilities for generating combinations and subsets.
 *
 * Provides helper template functions for:
 *   - Generating all m-combinations from a given vector.
 *   - Generating all non-empty subsets of a vector.
 *
 * Implemented using recursive backtracking with minimal dependencies.
 * Suitable for general-purpose combinatorial enumeration tasks.
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 *************************************************************************/

#ifndef INFOXTR_COMBN_HPP
#define INFOXTR_COMBN_HPP

#include <vector>
#include <limits>
#include <algorithm> 
#include <functional>

namespace infoxtr
{

namespace combn
{
    
    // ==============================
    // Combination generation
    // ==============================

    /**
     * @brief Generate all combinations of m elements from a given vector.
     *
     * Uses recursive backtracking to enumerate all size-m subsets
     * while preserving the original element order.
     *
     * @note: If m > vec.size() or m == 0, an empty vector is returned.
     *        This follows standard combinatorial mathematics (C(n, m) = 0).
     *
     * @tparam T Element type
     * @param vec Input vector
     * @param m Number of elements per combination
     * @return std::vector<std::vector<T>> All possible m-combinations
     */
    template <typename T>
    inline std::vector<std::vector<T>> combn(const std::vector<T>& vec, size_t m)
    {
        std::vector<std::vector<T>> result;
        std::vector<T> current;

        const size_t vec_size = vec.size();

        if (m == 0 || m > vec_size)
            return result;

        std::function<void(size_t)> helper = [&](size_t start)
        {
            if (current.size() == m)
            {
                result.push_back(current);
                return;
            }

            size_t remaining = m - current.size();

            for (size_t i = start; i + remaining <= vec_size; ++i)
            {
                current.push_back(vec[i]);
                helper(i + 1);
                current.pop_back();
            }
        };

        helper(0);
        return result;
    }

    // ==============================
    // Subset generation
    // ==============================

    /**
    * @brief Generate all non-empty subsets of a given vector up to a maximum order.
    *
    * Iteratively calls combn for sizes 1 to min(n, max_order),
    * where n is the size of the input vector.
    *
    * @tparam T Element type
    * @param vec Input vector
    * @param max_order Maximum subset size (default: no limit)
    * @return std::vector<std::vector<T>> All non-empty subsets up to max_order
    */
    template <typename T>
    inline std::vector<std::vector<T>> genSubsets(
        const std::vector<T>& vec,
        size_t max_order = std::numeric_limits<size_t>::max())
    {
        std::vector<std::vector<T>> allSubsets;
        const size_t n = vec.size();
        const size_t limit = std::min(n, max_order);

        for (size_t m = 1; m <= limit; ++m) {
            std::vector<std::vector<T>> combs = combn(vec, m);
            allSubsets.insert(allSubsets.end(), combs.begin(), combs.end());
        }

        return allSubsets;
    }

} // namespace combn

}

#endif // INFOXTR_COMBN_HPP
