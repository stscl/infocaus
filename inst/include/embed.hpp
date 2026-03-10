/*******************************************************************
 *  File: embed.hpp
 *
 *  Spatial embedding utilities for lattice/grid structures and 
 *  temporal-delay embedding utilities for time series.
 *
 *  This header provides high performance utilities for generating
 *  temporal embeddings, and for generating spatial embeddings on:
 *
 *    - Arbitrary lattices defined by neighbor lists
 *    - Regular 2D grids using Moore neighborhoods
 *
 *  Core capabilities:
 *
 *    - Multi lag neighbor expansion on graphs
 *    - Lagged value aggregation
 *    - Grid index mapping utilities
 *    - Temporal/Spatial embedding generation
 *
 *  Design principles:
 *
 *    - Type safety using std::size_t indices
 *    - No sentinel or magic values
 *    - Cache friendly memory layouts
 *    - Header only implementation
 *    - Robust NaN propagation
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 *******************************************************************/

#ifndef EMBED_HPP
#define EMBED_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <stdexcept>
#include <cstddef>
#include <utility>

namespace Embed 
{
    /* =============================================================
    * Type aliases
    * ============================================================= */

    using Index        = std::size_t;
    using NeighborList = std::vector<Index>;
    using NeighborMat  = std::vector<NeighborList>;
    using Vector       = std::vector<double>;
    using Matrix       = std::vector<Vector>;

    /* =============================================================
    * ------------------- LATTICE OPERATORS -----------------------
    * ============================================================= */

    /**
    * @brief Expand lattice neighbors to a given lag.
    *
    * This function computes lagged neighbor sets on a lattice graph.
    * Two expansion modes are supported:
    *
    *   cumulate = true:
    *       Returns the cumulative reachability closure
    *       R_lag(i) = nodes reachable from i within ≤ lag steps.
    *
    *   cumulate = false:
    *       Returns the exact lag shell
    *       S_lag(i) = nodes reachable in exactly lag steps,
    *                  excluding all nodes from previous lags.
    *
    * For lag = 0, each node returns itself.
    *
    * The algorithm recursively expands neighbor sets while preserving
    * sorted ordering and removing duplicates.
    *
    * Empty neighbor sets are represented as empty vectors.
    * No sentinel values are used.
    *
    * @param nb        Base adjacency list of the lattice
    * @param lag       Expansion lag (automatically clamped to n−1)
    * @param cumulate  Whether to return cumulative closure (true)
    *                  or exact lag shell (false)
    *
    * @return NeighborMat (std::vector<std::vector<size_t>>)
    *         Vector of lagged neighbor indices for each node
    *
    * @note
    *   - Safe for disconnected graphs
    *   - Stable ordering for deterministic results
    *   - lag ≥ n−1 produces saturated closure
    */
    inline NeighborMat LaggedNeighbors4Lattice(
        const NeighborMat& nb,
        size_t lag,
        bool cumulate = true
    ) {
        const size_t n = nb.size();
        NeighborMat result(n);

        const size_t max_lag = n - 1;
        lag = std::min(lag, max_lag);

        if (lag == 0) {
            for (size_t i = 0; i < n; ++i) {
                result[i] = { i };
            }
        } else {
            NeighborMat prev = LaggedNeighbors4Lattice(nb, lag - 1);
            
            if (cumulate) {
                for (size_t i = 0; i < n; ++i) {
                    if (prev[i].size() == n) {
                        result[i] = prev[i];
                    } else {
                        std::unordered_set<size_t> merged;
                        merged.reserve(prev[i].size() + nb[i].size());
                        for (size_t v : prev[i]) {
                            merged.insert(v);
                            for (size_t u : nb[v]) {
                                merged.insert(u);
                            }
                        }
                        result[i].assign(merged.begin(), merged.end());
                        std::sort(result[i].begin(), result[i].end());
                    }
                }
            } else {
                for (size_t i = 0; i < n; ++i) {
                    std::vector<size_t> newIndices;
                    for (size_t prev_nb : prev[i]) {
                        for (size_t cur_nb : nb[prev_nb]) {
                            if (!std::binary_search(prev[i].begin(), prev[i].end(), cur_nb)) {
                                newIndices.push_back(cur_nb);
                            }
                        }
                    }

                    if (newIndices.empty()) {
                        result[i] = { };
                    } else {
                        std::sort(newIndices.begin(), newIndices.end());
                        newIndices.erase(std::unique(newIndices.begin(), newIndices.end()),
                                        newIndices.end());
                        result[i] = std::move(newIndices);
                    }
                }
            }
        }

        return result;
    }

    /**
    * @brief Extract lagged values from a lattice vector.
    *
    * Each node collects values from its lagged neighbors (no recursively included).
    */
    inline Matrix LaggedValues4Lattice(
        const Vector& vec,
        const NeighborMat& nb,
        size_t lag
    ) {
        const size_t n = nb.size();
        Matrix out(n);

        const size_t max_lag = n - 1;
        lag = std::min(lag, max_lag);

        if (lag == 0) {
            for (size_t i = 0; i < n; ++i) {
                out[i] = { vec[i] };
            }
        } else {
            // Remove duplicates with previous lag (if lag > 1)
            NeighborMat prevNeighbors = LaggedNeighbors4Lattice(nb, lag-1);

            for (size_t i = 0; i < n; ++i) {
                // Convert previous lagged results to a set for fast lookup
                std::unordered_set<size_t> prevSet(prevNeighbors[i].begin(), prevNeighbors[i].end());
                // Remove duplicates from previous lagged results
                std::unordered_set<size_t> newIndices;
                for (size_t prev_nb : prevSet){
                    for (size_t cur_nb : nb[prev_nb]) {
                        if (prevSet.find(cur_nb) == prevSet.end()) {
                            newIndices.insert(cur_nb);
                        }
                    }
                }

                // Fill the lagged values
                if (newIndices.empty()) {
                    out[i] = { std::numeric_limits<double>::quiet_NaN() };
                } else {
                    out[i].reserve(newIndices.size());

                    // // Ignore order for value collection
                    // for (size_t j : newIndices) {
                    //     out[i].push_back(vec[j]);
                    // }

                    // Keep order for value collection
                    std::vector<size_t> sortedIndices(newIndices.begin(), newIndices.end());
                    std::sort(sortedIndices.begin(), sortedIndices.end());
                    for (size_t j : sortedIndices) {
                        out[i].push_back(vec[j]);
                    }
                }
            }
        }

        return out;
    }

    /**
    * @brief Generate lattice embeddings by averaging lagged neighbor values.
    *
    * Parameters:
    *   vec    Node values.
    *   nb     Neighbor list.
    *   E      Embedding dimension.
    *   tau    Lag step.
    *   style  0 include current state, otherwise exclude.
    *
    * Columns containing only NaN values are removed automatically.
    */
    inline Matrix GenLatticeEmbedding(
        const Vector& vec,
        const NeighborMat& nb,
        size_t E = 3,
        size_t tau = 1,
        size_t style = 1
    ) {
        const size_t n = vec.size();
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        size_t start = (style == 0 ? 0 : (tau == 0 ? 0 : tau));
        size_t step  = (tau == 0 ? 1 : tau);
        size_t end   = (tau == 0 ? E - 1 :
                    (style == 0 ? (E - 1) * tau : E * tau));
        end = std::min(end, n - 1);

        Matrix embed(n, Vector(E, NaN));
        std::unordered_map<size_t, NeighborMat> cache;
        cache.reserve(end + 1);

        if (start == 0){
            cache.emplace(start, LaggedNeighbors4Lattice(nb, start));
        } else {
            cache.emplace(start-1, LaggedNeighbors4Lattice(nb, start-1));
        }

        for (size_t lag = start; lag <= end; lag += 1) {
            if (cache.find(lag) == cache.end()){
                auto it = cache.find(lag-1);
                // assert(it != cache.end());
                const NeighborMat &prev = it->second;

                NeighborMat cur(n);
                for (size_t i = 0; i < n; ++i) {
                    if (prev[i].size() == n) {
                        cur[i] = prev[i];
                    } else {
                        std::unordered_set<size_t> merged;
                        merged.reserve(prev[i].size() + nb[i].size());
                        for (size_t v : prev[i]) {
                            merged.insert(v);
                            for (size_t u : nb[v]) {
                                merged.insert(u);
                            }
                        }
                        cur[i].assign(merged.begin(), merged.end());
                        std::sort(cur[i].begin(), cur[i].end());
                    }
                }

                cache.emplace(lag, std::move(cur));
            }
        }

        for (size_t lag = start; lag <= end; lag += step) {
            if (lag == 0) {
                for (size_t i = 0; i < n; ++i) {
                    embed[i][0] = vec[i];
                }
                continue;
            }

            const NeighborMat& cur  = cache.at(lag);
            const NeighborMat& prev = cache.at(lag - 1);

            const size_t col = (lag - start) / step;

            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                size_t cnt = 0;

                const auto& cur_set = cur[i];
                const auto& prev_set = prev[i];

                auto it_cur = cur_set.begin();
                auto it_prev = prev_set.begin();
                const auto it_prev_end = prev_set.end();

                // Compute set difference: cur \ prev (both sorted, prev ⊆ cur guaranteed)
                while (it_cur != cur_set.end() && it_prev != it_prev_end) {
                    if (*it_cur < *it_prev) {
                        const double v = vec[*it_cur];
                        if (!std::isnan(v)) { sum += v; ++cnt; }
                        ++it_cur;
                    } else if (*it_prev < *it_cur) {
                        ++it_prev;
                    } else {  // *it_cur == *it_prev → skip intersection
                        ++it_cur;
                        ++it_prev;
                    }
                }

                // Remaining elements in cur are guaranteed to be in the difference set
                // (since prev ⊆ cur and prev is exhausted)
                while (it_cur != cur_set.end()) {
                    const double v = vec[*it_cur];
                    if (!std::isnan(v)) { sum += v; ++cnt; }
                    ++it_cur;
                }

                if (cnt > 0) embed[i][col] = sum / cnt;
            }
        }

        /* Remove all-NaN columns */
        std::vector<size_t> validCols;
        for (size_t c = 0; c < embed.front().size(); ++c) {
            bool allNaN = true;
            for (size_t r = 0; r < embed.size(); ++r) {
                if (!std::isnan(embed[r][c])) {
                    allNaN = false;
                    break;
                }
            }
            if (!allNaN) validCols.push_back(c);
        }

        if (validCols.empty()) {
        throw std::invalid_argument(
            "No valid embeddings can be generated."
        );
        }
        if (validCols.size() == embed.front().size()) return embed;

        Matrix filtered(n);
        for (size_t r = 0; r < n; ++r) {
            filtered[r].reserve(validCols.size());
            for (size_t c : validCols) {
                filtered[r].push_back(embed[r][c]);
            }
        }
        return filtered;
    }

    /* =============================================================
    * -------------------- GRID OPERATORS -------------------------
    * ============================================================= */

    /**
    * @brief Convert 2D grid coordinates to linear index.
    */
    inline size_t GridIndex(size_t row, size_t col, size_t totalCols) {
        return row * totalCols + col;
    }

    /**
    * @brief Convert linear index to 2D grid coordinates.
    */
    inline std::pair<size_t, size_t> GridRowCol(size_t index, size_t totalCols) {
        return { index / totalCols, index % totalCols };
    }

    /**
    * @brief Flatten a grid matrix into a vector (row major).
    */
    inline Vector GridMatToVec(const Matrix& mat) {
        Vector out;
        if (mat.empty()) return out;

        const size_t rows = mat.size();
        const size_t cols = mat.front().size();
        out.reserve(rows * cols);

        for (const auto& row : mat) {
            out.insert(out.end(), row.begin(), row.end());
        }
        return out;
    }

    /**
    * @brief Reshape a vector into a grid matrix.
    */
    inline Matrix GridVecToMat(const Vector& vec, size_t nrow) {
        if (nrow == 0 || vec.size() % nrow != 0) {
            throw std::invalid_argument("GridVecToMat: incompatible dimensions.");
        }

        const size_t ncol = vec.size() / nrow;
        Matrix mat(nrow, Vector(ncol));

        for (size_t i = 0; i < nrow; ++i) {
            for (size_t j = 0; j < ncol; ++j) {
                mat[i][j] = vec[i * ncol + j];
            }
        }
        return mat;
    }

    /**
    * @brief Compute lagged Moore neighborhood values on a grid.
    */
    inline Matrix LaggedValues4Grid(
        const Matrix& mat,
        size_t lag
    ) {
        if (mat.empty() || mat.front().empty()) return {};

        const size_t rows = mat.size();
        const size_t cols = mat.front().size();
        const size_t total = rows * cols;
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        if (lag == 0) {
            Matrix out;
            out.reserve(total);
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    out.push_back({ mat[i][j] });
                }
            }
            return out;
        }

        // Since mat is non-empty and mat[0] is non-empty, rows >= 1 and cols >= 1
        const size_t maxLag = std::max(rows, cols) - 1;
        if (lag > maxLag) {
        lag = maxLag;
        }

        std::vector<std::pair<int,int>> offsets;
        for (int dx = -static_cast<int>(lag); dx <= static_cast<int>(lag); ++dx) {
            for (int dy = -static_cast<int>(lag); dy <= static_cast<int>(lag); ++dy) {
                if (std::max(std::abs(dx), std::abs(dy)) == static_cast<int>(lag)) {
                    offsets.emplace_back(dx, dy);
                }
            }
        }

        Matrix out(total, Vector(offsets.size(), NaN));

        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                const size_t id = GridIndex(r, c, cols);
                for (size_t k = 0; k < offsets.size(); ++k) {
                    const int nr = static_cast<int>(r) + offsets[k].first;
                    const int nc = static_cast<int>(c) + offsets[k].second;
                    if (nr >= 0 && nr < static_cast<int>(rows) &&
                        nc >= 0 && nc < static_cast<int>(cols)) {
                        out[id][k] = mat[nr][nc];
                    }
                }
            }
        }
        return out;
    }

    /**
    * @brief Generate grid spatial embeddings.
    */
    inline Matrix GenGridEmbedding(
        const Matrix& mat,
        size_t E = 3,
        size_t tau = 1,
        size_t style = 1
    ) {
        if (mat.empty() || mat.front().empty() || E == 0) return {};

        const size_t rows = mat.size();
        const size_t cols = mat.front().size();
        const size_t total = rows * cols;
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        /* -------------------------------------------------------
        Compute maximum meaningful lag for this grid
        Moore neighborhood cannot exceed this radius
        ------------------------------------------------------- */
        const size_t maxLag = std::max(rows, cols) - 1;

        /* -------------------------------------------------------
        Clamp embedding dimension E so that all generated lags
        are guaranteed to be within [0, maxLag] and overflow-free
        ------------------------------------------------------- */
        size_t maxE = 0;

        if (tau == 0) {
        /* lag(k) = k, k in [0, E-1] */
        maxE = maxLag + 1;
        }
        else if (style == 0) {
        /* lag(k) = k * tau, k in [1, E-1] */
        maxE = maxLag / tau + 1;
        }
        else {
        /* lag(k) = (k + 1) * tau, k in [0, E-1] */
        maxE = maxLag / tau;
        }

        if (maxE == 0) return {};
        if (E > maxE) E = maxE;

        Matrix embed(total, Vector(E, NaN));

        // /*
        //  * Repeats the grid data with a lag offset. The implementation is straightforward
        //  * but involves redundant computations for large datasets.
        //  * Kept here as a reference to the original approach.
        //  */
        // auto fill_column = [&](size_t col, size_t lag) {
        //   Matrix lagged = LaggedValues4Grid(mat, lag);
        //   for (size_t i = 0; i < lagged.size(); ++i) {
        //     double sum = 0.0;
        //     size_t cnt = 0;
        //     for (double v : lagged[i]) {
        //       if (!std::isnan(v)) {
        //         sum += v;
        //         ++cnt;
        //       }
        //     }
        //     if (cnt > 0) embed[i][col] = sum / cnt;
        //   }
        // };

        /* -------------------------
        Cache offsets per lag
        ------------------------- */
        std::unordered_map<size_t, std::vector<std::pair<int,int>>> offsetCache;
        offsetCache.reserve(E + 1);

        auto get_offsets = [&](size_t lag)
        -> const std::vector<std::pair<int,int>>&
        {
            auto it = offsetCache.find(lag);
            if (it != offsetCache.end()) return it->second;

            std::vector<std::pair<int,int>> offsets;
            if (lag == 0) {
            offsets.emplace_back(0, 0);
            } else {
            const int L = static_cast<int>(lag);
            for (int dx = -L; dx <= L; ++dx) {
                for (int dy = -L; dy <= L; ++dy) {
                if (std::max(std::abs(dx), std::abs(dy)) == L) {
                    offsets.emplace_back(dx, dy);
                }
                }
            }
            }

            auto res = offsetCache.emplace(lag, std::move(offsets));
            return res.first->second;
        };

        /* -------------------------
        Fill one embedding column
        ------------------------- */
        auto fill_column = [&](size_t col, size_t lag) {
        const auto& offsets = get_offsets(lag);

        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
            const size_t id = GridIndex(r, c, cols);

            double sum = 0.0;
            size_t cnt = 0;

            for (const auto& off : offsets) {
                const int nr = static_cast<int>(r) + off.first;
                const int nc = static_cast<int>(c) + off.second;
                if (nr >= 0 && nr < static_cast<int>(rows) &&
                    nc >= 0 && nc < static_cast<int>(cols)) {
                double v = mat[nr][nc];
                if (!std::isnan(v)) {
                    sum += v;
                    ++cnt;
                }
                }
            }

            if (cnt > 0) {
                embed[id][col] = sum / cnt;
            }
            // else keep NaN
            }
        }
        };

        /* -------------------------
        Generate embeddings
        ------------------------- */
        if (tau == 0) {
            Vector flat = GridMatToVec(mat);
            for (size_t i = 0; i < total; ++i) embed[i][0] = flat[i];
            for (size_t k = 1; k < E; ++k) fill_column(k, k);
        }
        else if (style == 0) {
            Vector flat = GridMatToVec(mat);
            for (size_t i = 0; i < total; ++i) embed[i][0] = flat[i];
            for (size_t k = 1; k < E; ++k) fill_column(k, k * tau);
        }
        else {
            for (size_t k = 0; k < E; ++k) fill_column(k, (k + 1) * tau);
        }

        /* remove all-NaN columns */
        std::vector<size_t> validCols;
        for (size_t c = 0; c < embed.front().size(); ++c) {
            bool allNaN = true;
            for (size_t r = 0; r < embed.size(); ++r) {
                if (!std::isnan(embed[r][c])) {
                    allNaN = false;
                    break;
                }
            }
            if (!allNaN) validCols.push_back(c);
        }

        if (validCols.empty()) {
        throw std::invalid_argument(
            "No valid embeddings can be generated."
        );
        }
        if (validCols.size() == embed.front().size()) return embed;

        Matrix filtered(total);
        for (size_t r = 0; r < total; ++r) {
            filtered[r].reserve(validCols.size());
            for (size_t c : validCols) {
                filtered[r].push_back(embed[r][c]);
            }
        }
        return filtered;
    }

    /* ========================================================
    * ------------------- TS OPERATORS -----------------------
    * ======================================================== */

    /**
    * @brief Generate time-delay embeddings for a univariate time series.
    *
    * This function reconstructs the state space of a scalar time series
    * using time-delay embedding with dimension E and lag tau.
    *
    * - When tau = 0, embedding uses lags of 0, 1, ..., E-1.  (original behavior)
    * - When tau > 0 and style = 1, embedding uses lags of tau, 2*tau, ..., E*tau.
    * - When tau > 0 and style = 0, embedding uses lags of 0, tau, 2*tau, ..., (E-1)*tau.
    *
    * Example:
    * Input: vec = {1, 2, 3, 4, 5}, E = 3, tau = 0
    * Output:
    * 1    NaN    NaN
    * 2    1      NaN
    * 3    2      1
    * 4    3      2
    * 5    4      3
    *
    * All values are pre-initialized to NaN. (Elements are filled only when
    * sufficient non-NaN lagged values are available. *Previously bound,
    * now abandoned*) Columns containing only NaN values are removed before
    * returning. If no valid embedding columns remain (due to short input
    * and large E/tau), an exception is thrown.
    *
    * @param vec The input time series as a vector of doubles.
    * @param E Embedding dimension.
    * @param tau Time lag.
    * @param style Lag style when tau > 0:
    *        - style = 1: tau, 2*tau, ..., E*tau
    *        - style = 0: 0, tau, 2*tau, ..., (E-1)*tau
    * @return A 2D vector (matrix) with valid embeddings (rows × cols).
    */
    inline Matrix GenTSEmbedding(
        const Vector& vec,
        size_t E = 3,
        size_t tau = 1,
        size_t style = 0
    ) {
    const size_t N = vec.size();

    if (E == 0) {
        throw std::invalid_argument("Embedding dimension E must be >= 1.");
    }

    // Compute the maximum required lag before embedding
    size_t max_lag;
    if (tau == 0) {
        max_lag = E - 1;
    } else if (style == 0) {
        max_lag = (E - 1) * tau;
    } else { // style == 1
        max_lag = E * tau;
    }

    // Pre-check: if the largest required lag exceeds available data
    if (max_lag >= N) {
        throw std::invalid_argument(
            "Embedding parameters require a lag larger than available data length."
        );
    }

    const double NaN = std::numeric_limits<double>::quiet_NaN();

    // Preallocate embedding matrix: N rows, E columns
    Matrix emb(N, std::vector<double>(E, NaN));

    for (size_t t = 0; t < N; ++t) {
        for (size_t j = 0; j < E; ++j) {
        size_t lag;

        if (tau == 0) {
            lag = j; // Original behavior: 0, 1, ..., E-1
        } else if (style == 0) {
            lag = j * tau;         // 0, tau, 2*tau, ..., (E-1)*tau
        } else { // style == 1;
            lag = (j + 1) * tau;   // tau, 2*tau, ..., E*tau
        }

        if(t >= lag) {
            emb[t][j] = vec[t - lag];
            // else leave NaN
        }
        }
    }

    // Check which columns contain at least one non-NaN value
    std::vector<size_t> keep;
    keep.reserve(E);
    for (size_t j = 0; j < E; ++j) {
        for (size_t i = 0; i < N; ++i) {
        if (!std::isnan(emb[i][j])) {
            keep.push_back(j);
            break;
        }
        }
    }

    if (keep.empty()) {
        throw std::invalid_argument(
            "No valid embeddings can be generated."
        );
    }
    if (keep.size() == E) return emb;

    // Create cleaned matrix with only columns having valid data
    Matrix cleaned(N);
    for (size_t i = 0; i < N; ++i) {
        cleaned[i].reserve(keep.size());
        for (size_t j : keep) {
        cleaned[i].push_back(emb[i][j]);
        }
    }

    return cleaned;
    }

} // namespace Embed

#endif // EMBED_HPP
