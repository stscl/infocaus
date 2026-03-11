/***************************************************************
 *  File: distance.hpp
 *
 *  High performance distance measurements
 *  for scalar, vector and matrix data.
 *
 *  Supported distance methods:
 *      "euclidean"  : sqrt(sum((x - y)^2))
 *      "maximum"    : max(|x - y|)
 *      "manhattan"  : sum(|x - y|)
 *
 *  NA handling:
 *      When na_rm = true:
 *          NaN values are removed pairwise before calculation.
 *          If all elements are removed, result is NaN.
 *
 *  Functions:
 *      Dist(scalar, scalar)
 *      Dist(vector, scalar)
 *      Dist(vector, vector)
 *      Dist(matrix)
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 ***************************************************************/

#ifndef DIST_HPP
#define DIST_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

namespace Dist
{   
    /*
     * @enum DistanceMethod
     * @brief Enumerated type for supported distance metrics in state space projection.
     * 
     * This enum provides a type-safe, efficient way to specify distance calculation
     * methods without repeated string comparisons in performance-critical loops.
     * 
     * @var Euclidean
     *   L2 norm: sqrt(sum((x_i - y_i)^2)). Sensitive to large single-dimension differences.
     * 
     * @var Manhattan
     *   L1 norm: sum(|x_i - y_i|). More robust to outliers than Euclidean.
     * 
     * @var Maximum
     *   L-infinity norm: max(|x_i - y_i|). Captures worst-case dimensional deviation.
     * 
     * @var Invalid
     *   Sentinel value indicating an unrecognized or unsupported method string.
     * 
     * @note Stored as uint8_t for minimal memory footprint and optimal switch dispatch.
     */
    enum class DistanceMethod : uint8_t {
        Euclidean,
        Manhattan,
        Maximum,
        Invalid
    };

    /*
     * @brief Parses a distance method name string into the corresponding DistanceMethod enum.
     * 
     * This helper function converts user-facing string identifiers (e.g., "euclidean")
     * into the internal enum representation. It is designed to be called exactly once
     * at the entry point of high-performance routines, eliminating repeated string
     * comparisons in nested loops.
     * 
     * @param method The distance method name as a string. Accepted values:
     *               - "euclidean" : L2 distance
     *               - "manhattan" : L1 distance
     *               - "maximum"   : Chebyshev / L-infinity distance
     * 
     * @return The corresponding DistanceMethod enum value. Returns DistanceMethod::Invalid
     *         if the input string does not match any supported method.
     * 
     * @note Case-sensitive matching. Whitespace or alternative spellings will result in Invalid.
     * @warning Caller must validate the return value != Invalid before proceeding with computation.
     * 
     * @example
     *   auto method = parseDistanceMethod("manhattan");
     *   if (method == DistanceMethod::Invalid) {
     *       throw std::invalid_argument("Unknown distance metric");
     *   }
     */
    inline DistanceMethod parseDistanceMethod(const std::string& method) {
        if (method == "euclidean") return DistanceMethod::Euclidean;
        if (method == "manhattan") return DistanceMethod::Manhattan;
        if (method == "maximum")   return DistanceMethod::Maximum;
        return DistanceMethod::Invalid;
    }

    /***********************************************************
     * Scalar - Scalar
     ***********************************************************/
    inline double Dist(
        const double scalar1,
        const double scalar2)
    {
        if (std::isnan(scalar1) || std::isnan(scalar2))
            return std::numeric_limits<double>::quiet_NaN();

        return std::abs(scalar1 - scalar2);
    }

    /***********************************************************
     * Scalar - Vector
     * Result length equals vector length
     ***********************************************************/
    inline std::vector<double> Dist(
        const double scalar,
        const std::vector<double>& vec)
    {
        std::vector<double> result(vec.size(),
            std::numeric_limits<double>::quiet_NaN());

        if (std::isnan(scalar)) return result;

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (!std::isnan(vec[i]))
            {
                result[i] = std::abs(scalar - vec[i]);
            }
        }

        return result;
    }

    inline std::vector<double> Dist(
        const std::vector<double>& vec,
        const double scalar)
    {
        return Dist(scalar, vec);
    }

    /***********************************************************
    * Vector - Scalar
    * Scalar is internally expanded to vector length
    * Result is a single double distance value
    ***********************************************************/
    inline double Dist(
        const std::vector<double>& vec,
        const double scalar,
        std::string method = "euclidean",
        bool na_rm = true)
    {
        if (vec.empty() || std::isnan(scalar))
            return std::numeric_limits<double>::quiet_NaN();

        const DistanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == DistanceMethod::Invalid) {
            throw std::invalid_argument("Unsupported distance method: " + method);
        }

        double sum = 0.0;
        double maxv = 0.0;
        size_t n_valid = 0;

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (na_rm && std::isnan(vec[i]))
                continue;

            if (!na_rm && std::isnan(vec[i]))
                return std::numeric_limits<double>::quiet_NaN();

            double diff = vec[i] - scalar;

            switch (dist_method) {
                case DistanceMethod::Euclidean:
                    sum += diff * diff;
                    break;
                case DistanceMethod::Manhattan:
                    sum += std::abs(diff);
                    break;
                case DistanceMethod::Maximum:
                    {
                        double ad = std::abs(diff);
                        if (ad > maxv) maxv = ad;
                    }
                    break;
                default:
                    break; 
            }

            ++n_valid;
        }

        if (n_valid == 0)
            return std::numeric_limits<double>::quiet_NaN();

        if (dist_method == DistanceMethod::Euclidean)
            return std::sqrt(sum);
        else if (dist_method == DistanceMethod::Manhattan)
            return sum;
        else
            return maxv;  // maximum
    }

    inline double Dist(
        const double scalar,
        const std::vector<double>& vec,
        std::string method = "euclidean",
        bool na_rm = true)
    {
        return Dist(vec, scalar, method, na_rm);
    }

    /***********************************************************
     * Vector - Vector
     * Element-wise distance
     ***********************************************************/
    inline std::vector<double> Dist(
        const std::vector<double>& vec1,
        const std::vector<double>& vec2)
    {
        if (vec1.size() != vec2.size())
            throw std::invalid_argument("Vectors must have equal length.");

        std::vector<double> result(vec1.size(),
            std::numeric_limits<double>::quiet_NaN());

        for (size_t i = 0; i < vec1.size(); ++i)
        {
            if (!std::isnan(vec1[i]) && !std::isnan(vec2[i]))
            {
                result[i] = std::abs(vec1[i] - vec2[i]);
            }
        }

        return result;
    }
    
    /***********************************************************
     * Vector - Vector
     * Result is a single double distance value
     ***********************************************************/
    inline double Dist(
        const std::vector<double>& vec1,
        const std::vector<double>& vec2,
        std::string method = "euclidean",
        bool na_rm = true)
    {   
        if (vec1.empty() || vec2.empty() || vec1.size() != vec2.size())
            return std::numeric_limits<double>::quiet_NaN();

        const DistanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == DistanceMethod::Invalid) {
            throw std::invalid_argument("Unsupported distance method: " + method);
        }

        double sum = 0.0;
        double maxv = 0.0;
        size_t n_valid = 0;

        for (size_t i = 0; i < vec1.size(); ++i)
        {   
            bool element_has_na = std::isnan(vec1[i]) || std::isnan(vec2[i]);

            if (element_has_na && na_rm)
                continue;

            if (element_has_na && !na_rm)
                return std::numeric_limits<double>::quiet_NaN();

            double diff = vec1[i] - vec2[i];

            switch (dist_method) {
                case DistanceMethod::Euclidean:
                    sum += diff * diff;
                    break;
                case DistanceMethod::Manhattan:
                    sum += std::abs(diff);
                    break;
                case DistanceMethod::Maximum:
                    {
                        double ad = std::abs(diff);
                        if (ad > maxv) maxv = ad;
                    }
                    break;
                default:
                    break; 
            }

            ++n_valid;
        }

        if (n_valid == 0)
            return std::numeric_limits<double>::quiet_NaN();

        if (dist_method == DistanceMethod::Euclidean)
            return std::sqrt(sum);
        else if (dist_method == DistanceMethod::Manhattan)
            return sum;
        else
            return maxv;  // maximum
    }

    /***********************************************************
     * Matrix
     * Row-wise distance matrix
     ***********************************************************/
    inline std::vector<std::vector<double>> Dist(
        const std::vector<std::vector<double>>& mat,
        std::string method = "euclidean",
        bool na_rm = true)
    {   
        if (mat.empty()) return {};

        const DistanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == DistanceMethod::Invalid) {
            throw std::invalid_argument("Unsupported distance method: " + method);
        }

        const size_t n = mat.size();

        std::vector<std::vector<double>> distm(
            n,
            std::vector<double>(n,
                std::numeric_limits<double>::quiet_NaN()));

        for (size_t i = 0; i < n; ++i) {
            distm[i][i] = 0.0;

            for (size_t j = i+1; j < n; ++j) 
            { 
                // double distv = Dist(mat[i], mat[j], method, na_rm);
                double distv = 0.0;
                
                double sum = 0.0;
                double maxv = 0.0;
                size_t n_valid = 0;

                for (size_t ei = 0; ei < mat[i].size(); ++ei)
                {   
                    bool element_has_na = std::isnan(mat[i][ei]) || std::isnan(mat[j][ei]);

                    if (element_has_na && na_rm) continue;

                    if (element_has_na && !na_rm)
                    {
                        distv = std::numeric_limits<double>::quiet_NaN();
                        break;
                    } 

                    double diff = mat[i][ei] - mat[j][ei];

                    switch (dist_method) {
                        case DistanceMethod::Euclidean:
                            sum += diff * diff;
                            break;
                        case DistanceMethod::Manhattan:
                            sum += std::abs(diff);
                            break;
                        case DistanceMethod::Maximum:
                            {
                                double ad = std::abs(diff);
                                if (ad > maxv) maxv = ad;
                            }
                            break;
                        default:
                            break; 
                    }

                    ++n_valid;
                }

                if (n_valid == 0 || std::isnan(distv)) continue;

                if (dist_method == DistanceMethod::Euclidean)
                    distv = std::sqrt(sum);
                else if (dist_method == DistanceMethod::Manhattan)
                    distv = sum;
                else
                    distv = maxv;  // maximum

                distm[i][j] = distv;  // Correctly assign distance to upper triangle
                distm[j][i] = distv;  // Mirror the value to the lower triangle
                // distm[i][j] = distm[j][i] = dist(mat[i], mat[j], method, na_rm);
            }
        }

        return distm;
    }

    /***********************************************************
     * Matrix - Subset Distance
     * Compute distances from pred rows to lib rows
     *
     * Each element (i, j) equals
     *      dist(mat[pred[i]], mat[lib[j]])
     ***********************************************************/
    inline std::vector<std::vector<double>> Dist(
        const std::vector<std::vector<double>>& mat,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        std::string method = "euclidean",
        bool na_rm = true)
    {
        if (mat.empty()) return {};

        const DistanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == DistanceMethod::Invalid) {
            throw std::invalid_argument("Unsupported distance method: " + method);
        }

        const size_t n_rows = mat.size();

        // // Validate indices
        // for (size_t idx : lib)
        // {
        //     if (idx >= n_rows)
        //         throw std::out_of_range("lib index out of range.");
        // }

        // for (size_t idx : pred)
        // {
        //     if (idx >= n_rows)
        //         throw std::out_of_range("pred index out of range.");
        // }

        std::vector<std::vector<double>> distm(
            n_rows,
            std::vector<double>(n_rows,
                std::numeric_limits<double>::quiet_NaN()));

        for (size_t i = 0; i < pred.size(); ++i)
        {
            const size_t pi = pred[i];

            for (size_t j = 0; j < lib.size(); ++j)
            {
                const size_t lj = lib[j];

                if (pi == lj) 
                {
                    distm[pi][lj] = 0.0;
                    continue;
                }

                // distm[pi][lj] = Dist(
                //     mat[pi],
                //     mat[lj],
                //     method,
                //     na_rm);

                double distv = 0.0;

                double sum = 0.0;
                double maxv = 0.0;
                size_t n_valid = 0;

                for (size_t ei = 0; ei < mat[pi].size(); ++ei)
                {   
                    bool element_has_na = std::isnan(mat[pi][ei]) || std::isnan(mat[lj][ei]);

                    if (element_has_na && na_rm) continue;

                    if (element_has_na && !na_rm) 
                    {
                        distv = std::numeric_limits<double>::quiet_NaN();
                        break;
                    } 

                    double diff = mat[pi][ei] - mat[lj][ei];

                    switch (dist_method) {
                        case DistanceMethod::Euclidean:
                            sum += diff * diff;
                            break;
                        case DistanceMethod::Manhattan:
                            sum += std::abs(diff);
                            break;
                        case DistanceMethod::Maximum:
                            {
                                double ad = std::abs(diff);
                                if (ad > maxv) maxv = ad;
                            }
                            break;
                        default:
                            break; 
                    }

                    ++n_valid;
                }

                if (n_valid == 0 || std::isnan(distv)) continue;
                
                if (dist_method == DistanceMethod::Euclidean)
                    distv = std::sqrt(sum);
                else if (dist_method == DistanceMethod::Manhattan)
                    distv = sum;
                else
                    distv = maxv;  // maximum

                distm[pi][lj] = distv;
            }
        }

        return distm;
    }

} // namespace Dist

#endif // DIST_HPP
