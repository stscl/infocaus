/********************************************************************************
 * File: numericutils.hpp
 *
 * Utility functions for safe and consistent floating-point operations.
 *
 * Provides helper functions for:
 *   - Floating-point comparison with combined relative and absolute tolerance.
 *   - Portable numeric constants (epsilon and tolerance).
 *   - Special mathematical functions such as the digamma function ψ(x).
 *
 * Intended for scientific computation where double precision stability matters.
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ********************************************************************************/

#ifndef INFOXTR_NUMERICUTILS_HPP
#define INFOXTR_NUMERICUTILS_HPP

#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <initializer_list>

namespace infoxtr
{

namespace numericutils
{

    // ==============================
    // Common numeric constants
    // ==============================

    constexpr double DOUBLE_EPS     = std::numeric_limits<double>::epsilon(); // ≈ 2.22e-16
    constexpr double DOUBLE_TOL_ABS = 1.5e-16;  // Absolute tolerance
    constexpr double DOUBLE_TOL_REL = 1.5e-8;   // Relative tolerance

    constexpr double PI = 3.141592653589793238462643383279502884; // π

    // ==============================
    // Floating-point comparison
    // ==============================

    /**
     * @brief Compare two double values with combined relative and absolute tolerance.
     *
     * Implements a numerically stable test for near equality:
     *
     * |x - y| <= max(rel_tol * max(|x|, |y|, 1.0), abs_tol)
     *
     * @param x First value
     * @param y Second value
     * @param rel_tol Relative tolerance (default DOUBLE_TOL_REL)
     * @param abs_tol Absolute tolerance (default DOUBLE_TOL_ABS)
     * @return true if x and y are considered equal within tolerance
     */
    inline bool doubleNearlyEqual(double x,
                                  double y,
                                  double rel_tol = DOUBLE_TOL_REL,
                                  double abs_tol = DOUBLE_TOL_ABS) noexcept
    {
        double diff  = std::fabs(x - y);
        double scale = std::max({1.0, std::fabs(x), std::fabs(y)});
        return diff <= std::max(rel_tol * scale, abs_tol);
    }

    // ==============================
    // Special mathematical functions
    // ==============================

    /**
     * @brief Computes the digamma function ψ(x) with high numerical accuracy.
     *
     * The digamma function is defined as
     *
     *      ψ(x) = d/dx log Γ(x)
     *
     * Implementation strategy
     *
     * 1. Reflection formula for negative arguments
     *
     *      ψ(1 - x) = ψ(x) + π cot(πx)
     *
     *    This allows evaluation for x < 0 except at poles
     *    (0, -1, -2, ...).
     *
     * 2. Recurrence relation
     *
     *      ψ(x) = ψ(x + 1) − 1/x
     *
     *    Used to shift small x values upward into a stable region.
     *
     * 3. Asymptotic Bernoulli expansion for large x
     *
     *      ψ(x) ≈ log(x) − 1/(2x)
     *             − 1/(12x²)
     *             + 1/(120x⁴)
     *             − 1/(252x⁶)
     *             + ...
     *
     * Accuracy is typically around 1e-15 for double precision inputs.
     *
     * @param x Input value
     * @return digamma value ψ(x)
     */
    inline double digamma(double x) noexcept
    {
        // Handle NaN
        if (std::isnan(x))
            return std::numeric_limits<double>::quiet_NaN();

        // Pole at non-positive integer
        if (std::floor(x) == x && x <= 0)
            return std::numeric_limits<double>::infinity();

        // Reflection formula for negative arguments
        if (x < 0.0)
        {
            return digamma(1.0 - x) - PI / std::tan(PI * x);
        }

        double result = 0.0;

        // Recurrence to shift x into asymptotic region
        while (x < 8.0)
        {
            result -= 1.0 / x;
            x += 1.0;
        }

        double inv_x  = 1.0 / x;
        double inv_x2 = inv_x * inv_x;

        // Bernoulli asymptotic expansion
        double series =
            inv_x2 * (-1.0 / 12.0 +
            inv_x2 * ( 1.0 / 120.0 +
            inv_x2 * (-1.0 / 252.0 +
            inv_x2 * ( 1.0 / 240.0 +
            inv_x2 * (-5.0 / 660.0 +
            inv_x2 * ( 691.0 / 32760.0 +
            inv_x2 * (-1.0 / 12.0 +
                      3617.0 / 8160.0)))))));

        return result + std::log(x) - 0.5 * inv_x + series;
    }

} // namespace numericutils

}

#endif // INFOXTR_NUMERICUTILS_HPP
