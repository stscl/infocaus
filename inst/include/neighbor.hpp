/********************************************************************
 *  File: neighbor.hpp
 *
 *  High performance nearest neighbor search utilities
 *  for matrix based state space and distance matrix data.
 *
 *  Core methods:
 *      NN4Mat
 *      NN4DistMat
 *
 *  Description:
 *      Provide k-nearest neighbor index search for:
 *          1. Raw feature matrix
 *          2. Precomputed distance matrix
 *
 *      Supports full search and pred/lib restricted search.
 *
 *  Distance methods (for NN4Mat):
 *      "euclidean"  : sqrt(sum((x - y)^2))
 *      "maximum"    : max(|x - y|)
 *      "manhattan"  : sum(|x - y|)
 *
 *  Data layout:
 *      mat         : std::vector<std::vector<double>>
 *
 *                    byrow = true
 *                       mat[row][dimension]
 *                       Each row is an observation
 *
 *                    byrow = false
 *                       mat[dimension][column]
 *                       Each column is an observation
 *
 *      distmat     : std::vector<std::vector<double>>
 *                    distance matrix
 *
 *      lib         : library row indices
 *
 *      pred        : prediction row indices
 *
 *  Output:
 *      std::vector<std::vector<size_t>>
 *      Each row contains indices of nearest neighbors.
 *
 *  Notes:
 *      Self matching is excluded by default.
 *      If include_self = true, self index will be placed first.
 *
 *      Self inclusion does NOT count toward the k limit when k > 0.
 *
 *      If fewer than k valid neighbors are available, 
 *      all valid candidates are returned.
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 ********************************************************************/

#ifndef NEIGHBOR_HPP
#define NEIGHBOR_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include "numericutils.hpp"
#include "distance.hpp"

namespace neighbor
{

    /***********************************************************
     * NN4Mat
     *
     * Returns:
     *      neighbor index list for each observation
     *
     * Algorithm:
     *      1. For each observation i
     *      2. Compute distances to all other observations
     *      3. Exclude self unless include_self = true
     *      4. Partial sort to obtain k nearest
     ***********************************************************/
    inline std::vector<std::vector<size_t>> NN4Mat(
        const std::vector<std::vector<double>>& mat,
        size_t k,
        std::string method = "euclidean",
        bool include_self = false,
        bool byrow = true)
    {
        const distance::distanceMethod dist_method =
            distance::parseDistanceMethod(method);

        if (dist_method == distance::distanceMethod::Invalid)
        {
            throw std::invalid_argument(
                "Unsupported distance method: " + method);
        }

        const size_t n_obs =
            byrow ? mat.size() : mat[0].size();

        const size_t dim =
            byrow ? mat[0].size() : mat.size();

        std::vector<std::vector<size_t>> result(n_obs);

        for (size_t i = 0; i < n_obs; ++i)
        {
            std::vector<std::pair<double,size_t>> candidates;
            candidates.reserve(n_obs);

            for (size_t j = 0; j < n_obs; ++j)
            {
                if (i == j) continue;

                double sum = 0.0;
                double maxv = 0.0;
                size_t n_valid = 0;

                for (size_t d = 0; d < dim; ++d)
                {
                    double xi = byrow ? mat[i][d] : mat[d][i];
                    double xj = byrow ? mat[j][d] : mat[d][j];

                    if (std::isnan(xi) || std::isnan(xj))
                        continue;

                    double diff = xi - xj;

                    switch (dist_method)
                    {
                        case distance::distanceMethod::Euclidean:
                            sum += diff * diff;
                            break;

                        case distance::distanceMethod::Manhattan:
                            sum += std::abs(diff);
                            break;

                        case distance::distanceMethod::Maximum:
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

                if (n_valid == 0) continue;

                double distv;

                if (dist_method == distance::distanceMethod::Euclidean)
                    distv = std::sqrt(sum);
                else if (dist_method == distance::distanceMethod::Manhattan)
                    distv = sum;
                else
                    distv = maxv;

                candidates.emplace_back(distv, j);
            }

            std::vector<size_t> indices;
            indices.reserve(k);

            size_t effective_k = k;

            if (include_self)
            {
                indices.push_back(i);
                if (effective_k > 0) effective_k--;
            }

            size_t num_neighbors =
                std::min(effective_k, candidates.size());

            if (num_neighbors > 0)
            {
                std::partial_sort(
                    candidates.begin(),
                    candidates.begin() + num_neighbors,
                    candidates.end(),
                    [](const auto& a, const auto& b)
                    {
                        if (!numericutils::doubleNearlyEqual(a.first,b.first))
                            return a.first < b.first;
                        return a.second < b.second;
                    });

                for (size_t m = 0; m < num_neighbors; ++m)
                    indices.push_back(candidates[m].second);
            }

            result[i] = std::move(indices);
        }

        return result;
    }

    /***********************************************************
     * NN4Mat (pred/lib restricted version)
     *
     * Returns:
     *      neighbor index list aligned to full matrix size,
     *      only pred position has values
     ***********************************************************/
    inline std::vector<std::vector<size_t>> NN4Mat(
        const std::vector<std::vector<double>>& mat,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t k,
        std::string method = "euclidean",
        bool include_self = false,
        bool byrow = true)
    {
        const distance::distanceMethod dist_method =
            distance::parseDistanceMethod(method);

        if (dist_method == distance::distanceMethod::Invalid)
        {
            throw std::invalid_argument(
                "Unsupported distance method: " + method);
        }

        const size_t n_obs =
            byrow ? mat.size() : mat[0].size();

        const size_t dim =
            byrow ? mat[0].size() : mat.size();

        std::vector<std::vector<size_t>> result(n_obs);

        std::unordered_set<size_t> lib_set(lib.begin(), lib.end());

        for (size_t i : pred)
        {

            std::vector<std::pair<double,size_t>> candidates;
            candidates.reserve(lib.size());

            bool self_in_lib =
                lib_set.find(i) != lib_set.end();

            for (size_t j : lib)
            {
                if (i == j) continue;

                double sum = 0.0;
                double maxv = 0.0;
                size_t n_valid = 0;

                for (size_t d = 0; d < dim; ++d)
                {
                    double xi = byrow ? mat[i][d] : mat[d][i];
                    double xj = byrow ? mat[j][d] : mat[d][j];

                    if (std::isnan(xi) || std::isnan(xj))
                        continue;

                    double diff = xi - xj;

                    switch (dist_method)
                    {
                        case distance::distanceMethod::Euclidean:
                            sum += diff * diff;
                            break;

                        case distance::distanceMethod::Manhattan:
                            sum += std::abs(diff);
                            break;

                        case distance::distanceMethod::Maximum:
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

                if (n_valid == 0) continue;

                double distv;

                if (dist_method == distance::distanceMethod::Euclidean)
                    distv = std::sqrt(sum);
                else if (dist_method == distance::distanceMethod::Manhattan)
                    distv = sum;
                else
                    distv = maxv;

                candidates.emplace_back(distv, j);
            }

            std::vector<size_t> indices;
            indices.reserve(k);

            size_t effective_k = k;

            if (include_self && self_in_lib)
            {
                indices.push_back(i);
                if (effective_k > 0) effective_k--;
            }

            size_t num_neighbors =
                std::min(effective_k, candidates.size());

            if (num_neighbors > 0)
            {
                std::partial_sort(
                    candidates.begin(),
                    candidates.begin() + num_neighbors,
                    candidates.end(),
                    [](const auto& a, const auto& b)
                    {
                        if (!numericutils::doubleNearlyEqual(a.first,b.first))
                            return a.first < b.first;
                        return a.second < b.second;
                    });

                for (size_t m = 0; m < num_neighbors; ++m)
                    indices.push_back(candidates[m].second);
            }

            result[i] = std::move(indices);
        }

        return result;
    }

    /***********************************************************
     * NN4DistMat
     *
     * Returns:
     *      neighbor index list
     ***********************************************************/
    inline std::vector<std::vector<size_t>> NN4DistMat(
        const std::vector<std::vector<double>>& distmat,
        size_t k,
        bool include_self = false)
    {
      const size_t n = distmat.size();

      // Initialize result with empty vectors
      std::vector<std::vector<size_t>> result(n);

      for (size_t i = 0; i < n; ++i) {
        const auto& row = distmat[i];

        // if (std::isnan(row[i])) continue;

        std::vector<std::pair<double, size_t>> candidates;
        candidates.reserve(n);

        for (size_t j = 0; j < n; ++j) {

          if (i == j) continue;

          double d = row[j];
          if (!std::isnan(d)) {
            candidates.emplace_back(d, j);
          }
        }

        std::vector<size_t> indices;
        indices.reserve(k);
        size_t effective_k = k;

        // Handle self inclusion
        if (include_self) {
          indices.push_back(i);
          if (effective_k > 0) {
            effective_k -= 1;
          }
        }

        size_t num_neighbors = std::min(effective_k, candidates.size());

        if (num_neighbors > 0) {

          std::partial_sort(
            candidates.begin(),
            candidates.begin() + num_neighbors,
            candidates.end(),
            [](const std::pair<double, size_t>& a,
              const std::pair<double, size_t>& b) {
              if (!numericutils::doubleNearlyEqual(a.first, b.first)) {
                return a.first < b.first;
              } else {
                return a.second < b.second;
              }
            }
          );

          for (size_t m = 0; m < num_neighbors; ++m) {
            indices.push_back(candidates[m].second);
          }
        }

        result[i] = std::move(indices);
      }

      return result;
    }

    /***********************************************************
     * NN4DistMat (pred/lib restricted version)
     *
     * Returns:
     *      neighbor index list aligned to full matrix size,
     *      only pred position has values
     ***********************************************************/
    inline std::vector<std::vector<size_t>> NN4DistMat(
        const std::vector<std::vector<double>>& distmat,
        const std::vector<size_t>& lib,
        const std::vector<size_t>& pred,
        size_t k,
        bool include_self = false)
    { 
      const size_t n = distmat.size();

      // Initialize result with empty vectors
      std::vector<std::vector<size_t>> result(n);

      // Build lib membership lookup
      std::unordered_set<size_t> lib_set(lib.begin(), lib.end());

      // if (k > lib.size()) {
      //   throw std::invalid_argument("Invalid argument: k exceeds library set capacity " + std::to_string(lib.size()));
      // }

      for (size_t i : pred) {

        // if (i >= n || distmat[i].size() != n) continue;

        const auto& row = distmat[i];

        // if (std::isnan(row[i])) continue;

        std::vector<std::pair<double, size_t>> candidates;
        candidates.reserve(lib.size());

        bool self_in_lib = (lib_set.find(i) != lib_set.end());

        // Collect valid neighbors from lib excluding self
        for (size_t j : lib) {

          if (i == j) continue;

          double d = row[j];
          if (!std::isnan(d)) {
            candidates.emplace_back(d, j);
          }
        }

        std::vector<size_t> indices;
        indices.reserve(k);
        size_t effective_k = k;

        // Handle self inclusion
        if (include_self && self_in_lib) {
          indices.push_back(i);
          if (effective_k > 0) {
            effective_k -= 1;
          }
        }

        size_t num_neighbors = std::min(effective_k, candidates.size());

        if (num_neighbors > 0) {

          std::partial_sort(
            candidates.begin(),
            candidates.begin() + num_neighbors,
            candidates.end(),
            [](const std::pair<double, size_t>& a,
              const std::pair<double, size_t>& b) {
              if (!numericutils::doubleNearlyEqual(a.first, b.first)) {
                return a.first < b.first;
              } else {
                return a.second < b.second;
              }
            }
          );

          for (size_t m = 0; m < num_neighbors; ++m) {
            indices.push_back(candidates[m].second);
          }
        }

        result[i] = std::move(indices);
      }

      return result;
    }

} // namespace neighbor

#endif // NEIGHBOR_HPP
