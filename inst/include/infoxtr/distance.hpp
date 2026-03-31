/************************************************************************
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
 *  Matrix behavior:
 *      Distances can be computed either row-wise or column-wise.
 *
 *      byrow = true  (default)
 *          Each row is treated as an observation vector.
 *          Distance matrix is computed between rows.
 *
 *      byrow = false
 *          Each column is treated as an observation vector.
 *          Distance matrix is computed between columns.
 *
 *  Functions:
 *      distance(scalar, scalar)
 *      distance(vector, scalar)
 *      distance(vector, vector)
 *      distance(vector)                     -> full distance matrix
 *      distance(matrix)                     -> full distance matrix
 *      distance(matrix, lib, pred)          -> subset distance matrix
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 ************************************************************************/

#ifndef INFOXTR_DISTANCE_HPP
#define INFOXTR_DISTANCE_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

namespace infoxtr 
{

namespace distance
{   
    /*
     * @enum distanceMethod
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
    enum class distanceMethod : uint8_t {
        Euclidean,
        Manhattan,
        Maximum,
        Invalid
    };

    /*
     * @brief Parses a distance method name string into the corresponding distanceMethod enum.
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
     * @return The corresponding distanceMethod enum value. Returns distanceMethod::Invalid
     *         if the input string does not match any supported method.
     * 
     * @note Case-sensitive matching. Whitespace or alternative spellings will result in Invalid.
     * @warning Caller must validate the return value != Invalid before proceeding with computation.
     * 
     * @example
     *   auto method = parseDistanceMethod("manhattan");
     *   if (method == distanceMethod::Invalid) {
     *       throw std::invalid_argument("Unknown distance metric");
     *   }
     */
    inline distanceMethod parseDistanceMethod(const std::string& method) {
        if (method == "euclidean") return distanceMethod::Euclidean;
        if (method == "manhattan") return distanceMethod::Manhattan;
        if (method == "maximum")   return distanceMethod::Maximum;
        return distanceMethod::Invalid;
    }

    /***********************************************************
     * Scalar - Scalar
     ***********************************************************/
    inline double distance(
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
    inline std::vector<double> distance(
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

    inline std::vector<double> distance(
        const std::vector<double>& vec,
        const double scalar)
    {
        return distance(scalar, vec);
    }

    /***********************************************************
    * Vector - Scalar
    * Scalar is internally expanded to vector length
    * Result is a single double distance value
    ***********************************************************/
    inline double distance(
        const std::vector<double>& vec,
        const double scalar,
        std::string method = "euclidean",
        bool na_rm = true)
    {
        if (vec.empty() || std::isnan(scalar))
            return std::numeric_limits<double>::quiet_NaN();

        const distanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == distanceMethod::Invalid) {
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
                case distanceMethod::Euclidean:
                    sum += diff * diff;
                    break;
                case distanceMethod::Manhattan:
                    sum += std::abs(diff);
                    break;
                case distanceMethod::Maximum:
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

        if (dist_method == distanceMethod::Euclidean)
            return std::sqrt(sum);
        else if (dist_method == distanceMethod::Manhattan)
            return sum;
        else
            return maxv;  // maximum
    }

    inline double distance(
        const double scalar,
        const std::vector<double>& vec,
        std::string method = "euclidean",
        bool na_rm = true)
    {
        return distance(vec, scalar, method, na_rm);
    }

    /***********************************************************
     * Vector - Vector
     * Element-wise distance
     ***********************************************************/
    inline std::vector<double> distance(
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
    inline double distance(
        const std::vector<double>& vec1,
        const std::vector<double>& vec2,
        std::string method = "euclidean",
        bool na_rm = true)
    {   
        if (vec1.empty() || vec2.empty() || vec1.size() != vec2.size())
            return std::numeric_limits<double>::quiet_NaN();

        const distanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == distanceMethod::Invalid) {
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
                case distanceMethod::Euclidean:
                    sum += diff * diff;
                    break;
                case distanceMethod::Manhattan:
                    sum += std::abs(diff);
                    break;
                case distanceMethod::Maximum:
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

        if (dist_method == distanceMethod::Euclidean)
            return std::sqrt(sum);
        else if (dist_method == distanceMethod::Manhattan)
            return sum;
        else
            return maxv;  // maximum
    }

    /****************************************************************************
     * Vector Distance
     * Computes a full pairwise distance matrix from a numeric vector.
     ***************************************************************************/
    inline std::vector<std::vector<double>> distance(
        const std::vector<double>& vec)
    {
        if (vec.empty()) return {};

        const size_t n = vec.size();

        std::vector<std::vector<double>> distm(
            n,
            std::vector<double>(n,
                std::numeric_limits<double>::quiet_NaN()));

        for (size_t i = 0; i < n; ++i)
        {   
            if (std::isnan(vec[i])) continue;

            distm[i][i] = 0.0;

            for (size_t j = i + 1; j < n; ++j)
            {
                if (!std::isnan(vec[j]))
                {
                    double distv = std::abs(vec[i] - vec[j]);
                    distm[i][j] = distv;
                    distm[j][i] = distv;
                }
            }
        }

        return distm;
    }

    /****************************************************************************
     * Matrix Distance
     *
     * Computes a full pairwise distance matrix from a numeric matrix.
     *
     * Behavior depends on the byrow parameter:
     *
     * byrow = true  (default)
     *      Each row is treated as an observation vector.
     *      Distances are computed between rows.
     *      Result size: n_rows × n_rows
     *
     * byrow = false
     *      Each column is treated as an observation vector.
     *      Distances are computed between columns.
     *      Result size: n_cols × n_cols
     *
     * @param mat    Input numeric matrix stored as vector of rows
     * @param method Distance metric ("euclidean", "manhattan", "maximum")
     * @param na_rm  Remove NaN/NA values pairwise if true
     * @param byrow  If true compute row distances, otherwise column distances
     *
     * @return Symmetric distance matrix
     ***************************************************************************/
    inline std::vector<std::vector<double>> distance(
        const std::vector<std::vector<double>>& mat,
        std::string method = "euclidean",
        bool na_rm = true,
        bool byrow = true)
    {
        if (mat.empty()) return {};

        const distanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == distanceMethod::Invalid) {
            throw std::invalid_argument("Unsupported distance method: " + method);
        }

        const size_t n_rows = mat.size();
        const size_t n_cols = mat[0].size();

        const size_t n = byrow ? n_rows : n_cols;

        std::vector<std::vector<double>> distm(
            n,
            std::vector<double>(n,
                std::numeric_limits<double>::quiet_NaN()));

        for (size_t i = 0; i < n; ++i)
        {
            distm[i][i] = 0.0;

            for (size_t j = i + 1; j < n; ++j)
            {
                double sum = 0.0;
                double maxv = 0.0;
                size_t n_valid = 0;
                bool has_na = false;

                const size_t dim = byrow ? n_cols : n_rows;

                for (size_t k = 0; k < dim; ++k)
                {
                    double xi = byrow ? mat[i][k] : mat[k][i];
                    double xj = byrow ? mat[j][k] : mat[k][j];

                    bool element_has_na = std::isnan(xi) || std::isnan(xj);

                    if (element_has_na && na_rm) continue;

                    if (element_has_na && !na_rm)
                    {
                        has_na = true;
                        break;
                    }

                    double diff = xi - xj;

                    switch (dist_method) {
                        case distanceMethod::Euclidean:
                            sum += diff * diff;
                            break;

                        case distanceMethod::Manhattan:
                            sum += std::abs(diff);
                            break;

                        case distanceMethod::Maximum:
                        {
                            double ad = std::abs(diff);
                            if (ad > maxv) maxv = ad;
                            break;
                        }

                        default:
                            break;
                    }

                    ++n_valid;
                }

                if (has_na || n_valid == 0)
                    continue;

                double distv;

                if (dist_method == distanceMethod::Euclidean)
                    distv = std::sqrt(sum);
                else if (dist_method == distanceMethod::Manhattan)
                    distv = sum;
                else
                    distv = maxv;

                distm[i][j] = distv;
                distm[j][i] = distv;
            }
        }

        return distm;
    }

    /****************************************************************************************
     * Matrix Subset Distance
     *
     * Computes distances between selected rows or columns
     * of a matrix using index sets.
     *
     * Each element (pi, lj) equals
     *
     *      distance(mat[pred[i]], mat[lib[j]])
     *
     * Behavior depends on byrow:
     *
     * byrow = true
     *      Distances are computed between rows.
     *
     * byrow = false
     *      Distances are computed between columns.
     *
     * @param mat    Input numeric matrix
     * @param lib    Indices defining the library set
     * @param pred   Indices defining the prediction set
     * @param method Distance metric
     * @param na_rm  Remove NaN/NA values pairwise
     * @param byrow  If true operate on rows, otherwise columns
     *
     * @return A square matrix of size n_rows × n_rows (or n_cols × n_cols if byrow=false),
     *         with entries at [pred[i]][lib[j]] filled with computed distances.
     *         All other entries remain NaN.
     ***************************************************************************************/
    inline std::vector<std::vector<double>> distance(
        const std::vector<std::vector<double>>& mat,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        std::string method = "euclidean",
        bool na_rm = true,
        bool byrow = true)
    {
        if (mat.empty()) return {};

        const distanceMethod dist_method = parseDistanceMethod(method);
        if (dist_method == distanceMethod::Invalid) {
            throw std::invalid_argument("Unsupported distance method: " + method);
        }

        const size_t n_rows = mat.size();
        const size_t n_cols = mat[0].size();

        const size_t n = byrow ? n_rows : n_cols;

        std::vector<std::vector<double>> distm(
            n,
            std::vector<double>(n,
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

                double sum = 0.0;
                double maxv = 0.0;
                size_t n_valid = 0;
                bool has_na = false;

                const size_t dim = byrow ? n_cols : n_rows;

                for (size_t k = 0; k < dim; ++k)
                {
                    double xi = byrow ? mat[pi][k] : mat[k][pi];
                    double xj = byrow ? mat[lj][k] : mat[k][lj];

                    bool element_has_na = std::isnan(xi) || std::isnan(xj);

                    if (element_has_na && na_rm) continue;

                    if (element_has_na && !na_rm)
                    {
                        has_na = true;
                        break;
                    }

                    double diff = xi - xj;

                    switch (dist_method) {
                        case distanceMethod::Euclidean:
                            sum += diff * diff;
                            break;

                        case distanceMethod::Manhattan:
                            sum += std::abs(diff);
                            break;

                        case distanceMethod::Maximum:
                        {
                            double ad = std::abs(diff);
                            if (ad > maxv) maxv = ad;
                            break;
                        }

                        default:
                            break;
                    }

                    ++n_valid;
                }

                if (has_na || n_valid == 0)
                    continue;

                double distv;

                if (dist_method == distanceMethod::Euclidean)
                    distv = std::sqrt(sum);
                else if (dist_method == distanceMethod::Manhattan)
                    distv = sum;
                else
                    distv = maxv;

                distm[pi][lj] = distv;
            }
        }

        return distm;
    }

} // namespace distance

}

#endif // INFOXTR_DISTANCE_HPP
