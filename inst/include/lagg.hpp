/*******************************************************************
 *  File: lagg.hpp
 *
 *  Lag utilities for matrices.
 *
 *  Provides lagged aggregation operators for:
 *
 *    - Arbitrary lattice graphs
 *    - Regular grids
 *    - Time series
 *
 *  Each column of the input matrix is treated as a variable.
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 *******************************************************************/

#ifndef LAGG_HPP
#define LAGG_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <stdexcept>
#include <cstddef>

namespace lagg
{

    using Index        = std::size_t;
    using NeighborList = std::vector<Index>;
    using NeighborMat  = std::vector<NeighborList>;
    using Vector       = std::vector<double>;
    using Matrix       = std::vector<Vector>;

    static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    inline size_t gridIndex(size_t r, size_t c, size_t ncol)
    {
        return r * ncol + c;
    }

    /* ============================================================
     * LATTICE LAG
     * ============================================================ */

    inline Matrix lagg(
        const Matrix& mat,
        const NeighborMat& nb,
        size_t lag
    )
    {
        const size_t n = nb.size();
        const size_t p = mat.front().size();

        Matrix out(n, Vector(p, NaN));

        if (lag == 0) {
            return mat;
        }

        for (size_t i = 0; i < n; ++i)
        {
            std::vector<char> visited(n, 0);
            std::vector<size_t> frontier;
            std::vector<size_t> next;

            frontier.push_back(i);
            visited[i] = 1;

            for (size_t step = 1; step <= lag; ++step)
            {
                next.clear();

                for (size_t v : frontier)
                {
                    for (size_t u : nb[v])
                    {
                        if (!visited[u])
                        {
                            visited[u] = 1;
                            next.push_back(u);
                        }
                    }
                }

                if (step == lag)
                {
                    for (size_t j = 0; j < p; ++j)
                    {
                        double sum = 0.0;
                        size_t cnt = 0;

                        for (size_t id : next)
                        {
                            double v = mat[id][j];
                            if (!std::isnan(v))
                            {
                                sum += v;
                                ++cnt;
                            }
                        }

                        if (cnt > 0)
                            out[i][j] = sum / cnt;
                    }
                }

                frontier.swap(next);

                if (frontier.empty())
                    break;
            }
        }

        return out;
    }

    /* ============================================================
     *  GRID LAG
     * ============================================================ */

    inline Matrix lagg(
        const Matrix& mat,
        size_t nrows,
        size_t lag
    )
    {
        const size_t N = mat.size();
        const size_t p = mat.front().size();

        if (nrows == 0 || N % nrows != 0)
            throw std::invalid_argument("Invalid grid dimensions.");

        const size_t ncols = N / nrows;

        Matrix out(N, Vector(p, NaN));

        if (lag == 0)
            return mat;

        std::vector<std::pair<int,int>> offsets;

        const int L = static_cast<int>(lag);

        for (int dx = -L; dx <= L; ++dx)
        {
            for (int dy = -L; dy <= L; ++dy)
            {
                if (std::max(std::abs(dx), std::abs(dy)) == L)
                    offsets.emplace_back(dx, dy);
            }
        }

        for (size_t r = 0; r < nrows; ++r)
        {
            for (size_t c = 0; c < ncols; ++c)
            {
                const size_t id = gridIndex(r, c, ncols);

                for (size_t j = 0; j < p; ++j)
                {
                    double sum = 0.0;
                    size_t cnt = 0;

                    for (const auto& off : offsets)
                    {
                        int nr = static_cast<int>(r) + off.first;
                        int nc = static_cast<int>(c) + off.second;

                        if (nr >= 0 && nr < (int)nrows &&
                            nc >= 0 && nc < (int)ncols)
                        {
                            size_t nid = gridIndex(nr, nc, ncols);
                            double v = mat[nid][j];

                            if (!std::isnan(v))
                            {
                                sum += v;
                                ++cnt;
                            }
                        }
                    }

                    if (cnt > 0)
                        out[id][j] = sum / cnt;
                }
            }
        }

        return out;
    }

    /* ============================================================
     *  TIME SERIES LAG
     * ============================================================ */

    inline Matrix lagg(
        const Matrix& mat,
        size_t lag
    )
    {
        const size_t n = mat.size();
        const size_t p = mat.front().size();

        Matrix out(n, Vector(p, NaN));

        if (lag == 0)
            return mat;

        for (size_t t = lag; t < n; ++t)
        {
            for (size_t j = 0; j < p; ++j)
            {
                double v = mat[t - lag][j];
                if (!std::isnan(v))
                    out[t][j] = v;
            }
        }

        return out;
    }

} // namespace lagg

#endif // LAGG_HPP
