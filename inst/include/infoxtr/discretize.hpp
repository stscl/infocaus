/**********************************************************************
 *  File: discretize.hpp
 *
 *  Discretization utilities for continuous numeric vectors.
 *
 *  Methods implemented
 *      ---------------------------------------------
 *      sd        Standard deviation discretization
 *      equal     Equal interval discretization
 *      geometric Geometric interval discretization
 *      quantile  Quantile discretization
 *      manual    Manual discretization
 *      natural   Jenks natural breaks
 *      headtail  Head/Tail breaks
 *
 *  Input
 *      const std::vector<double>& vec
 *
 *  Output
 *      std::vector<uint64_t>
 *          0  -> NaN / NA values
 *          1..n -> discretized classes
 *
 *  Behaviour
 *      NaN values are ignored when computing breaks
 *      but assigned class 0 in the output.
 *
 *  RNG
 *      Sampling uses std::mt19937_64
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 **********************************************************************/

#ifndef INFOXTR_DISCRETIZE_HPP
#define INFOXTR_DISCRETIZE_HPP

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <limits>
#include <string>
#include <iostream>
#include <stdexcept>
#include "infoxtr/numericutils.hpp"

namespace infoxtr
{

namespace discretize
{

    /***********************************************************
     * Utility helpers
     ***********************************************************/

    inline std::vector<double> remove_nan(const std::vector<double>& v)
    {
        std::vector<double> out;
        out.reserve(v.size());

        for (double x : v)
        {
            if (!std::isnan(x))
                out.push_back(x);
        }

        if (out.empty())
            throw std::invalid_argument(
                "[Discretize] Input vector contains no valid numeric values (all values are NaN).");

        return out;
    }

    inline double mean(const std::vector<double>& v)
    {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    }

    inline double stddev(const std::vector<double>& v)
    {   
        if (v.size() <= 1)
            return 0.0;

        double m = mean(v);
        double s = 0.0;

        for (double x : v)
            s += (x - m) * (x - m);

        return std::sqrt(s / (v.size() - 1));
    }

    inline double min_val(const std::vector<double>& v)
    {
        return *std::min_element(v.begin(), v.end());
    }

    inline double max_val(const std::vector<double>& v)
    {
        return *std::max_element(v.begin(), v.end());
    }

    /***********************************************************
     * Standard deviation discretization
     ***********************************************************/
    inline std::vector<uint64_t> sdDisc(
        const std::vector<double>& vec,
        size_t n,
        bool right_closed = true)
    {
        auto x = remove_nan(vec);
        double m = mean(x);
        double sd = stddev(x);

        std::vector<uint64_t> res(vec.size(), 0);

        if (infoxtr::numericutils::doubleNearlyEqual(sd, 0.0)) 
        {
            for (size_t i = 0; i < vec.size(); ++i)
                if (!std::isnan(vec[i])) res[i] = 1;
            return res;
        }

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;
            
            double z = (vec[i] - m) / sd + n / 2.0;

            long idx;
            if (right_closed)
                idx = std::ceil(z);
            else
                idx = std::floor(z) + 1;
            idx = std::max<long>(1, std::min<long>(idx, static_cast<long>(n)));

            res[i] = static_cast<uint64_t>(idx);
        }

        return res;
    }

    /***********************************************************
     * Equal interval discretization
     ***********************************************************/
    inline std::vector<uint64_t> equalDisc(
        const std::vector<double>& vec,
        size_t n,
        bool right_closed = true)
    {
        auto x = remove_nan(vec);
        double minx = min_val(x);
        double maxx = max_val(x);

        double interval = (maxx - minx) / n;

        std::vector<uint64_t> res(vec.size(), 0);

        if (infoxtr::numericutils::doubleNearlyEqual(interval, 0.0)) 
        {
            for (size_t i = 0; i < vec.size(); ++i)
                if (!std::isnan(vec[i])) res[i] = 1;
            return res;
        }

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;
            
            double val = (vec[i] - minx) / interval;

            long idx;
            if (right_closed)
                idx = std::ceil(val);
            else
                idx = std::floor(val) + 1;
            idx = std::max<long>(1, std::min<long>(idx, static_cast<long>(n)));

            res[i] = static_cast<uint64_t>(idx);
        }

        return res;
    }

    /***********************************************************
     * Geometric discretization
     ***********************************************************/
    inline std::vector<uint64_t> geometricDisc(
        const std::vector<double>& vec,
        size_t n,
        bool right_closed = true)
    {
        auto x = remove_nan(vec);
        double minx = min_val(x);
        double maxx = max_val(x);

        if (minx <= 0)
            throw std::invalid_argument(
                "[Discretize] geometricDisc requires strictly positive data");

        double factor = std::pow(maxx / minx, 1.0 / n);

        std::vector<uint64_t> res(vec.size(), 0);

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;

            long idx;
            if (right_closed)
                idx = std::ceil(std::log(vec[i] / minx) / std::log(factor));
            else
                idx = std::floor(std::log(vec[i] / minx) / std::log(factor)) + 1;
            idx = std::max<long>(1, std::min<long>(idx, static_cast<long>(n)));

            res[i] = static_cast<uint64_t>(idx);
        }

        return res;
    }

    /***********************************************************
     * Quantile discretization
     ***********************************************************/
    inline std::vector<uint64_t> quantileDisc(
        const std::vector<double>& vec,
        size_t n,
        bool right_closed = true)
    {
        auto x = remove_nan(vec);

        std::vector<double> sorted = x;
        std::sort(sorted.begin(), sorted.end());

        std::vector<double> q(n + 1);

        for (size_t i = 0; i <= n; ++i)
        {
            double pos = i * (sorted.size() - 1) / static_cast<double>(n);
            size_t lo = static_cast<size_t>(std::floor(pos));
            size_t hi = std::min(lo + 1, sorted.size() - 1);
            double frac = pos - lo;
            q[i] = sorted[lo] * (1 - frac) + sorted[hi] * frac;
        }

        std::vector<uint64_t> res(vec.size(), 0);

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;

            bool assigned = false;

            for (size_t j = 0; j < n; ++j)
            {   
                bool in_bin = right_closed ? (vec[i] <= q[j + 1])
                                        : (vec[i] < q[j + 1]);
                if (in_bin)
                {
                    res[i] = static_cast<uint64_t>(j + 1);
                    assigned = true;
                    break;
                }
            }

            if (!assigned)
                res[i] = static_cast<uint64_t>(n);

        }

        return res;
    }

    /***********************************************************
     * Manual breakpoints discretization
     ***********************************************************/
    inline std::vector<uint64_t> manualDisc(
        const std::vector<double>& vec,
        const std::vector<double>& breakpoints,
        bool right_closed = true)
    {
        if (breakpoints.empty())
            throw std::invalid_argument("[Discretize] manualDisc: breakpoints cannot be empty");

        auto x = remove_nan(vec);

        std::vector<double> bp = breakpoints;

        std::sort(bp.begin(), bp.end());

        std::vector<uint64_t> res(vec.size(), 0);

        size_t n = bp.size() + 1;

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;

            bool assigned = false;

            for (size_t j = 0; j < bp.size(); ++j)
            {   
                bool classify_val = right_closed ? (vec[i] <= bp[j])
                                                : (vec[i] < bp[j]);
                if (classify_val)
                {
                    res[i] = static_cast<uint64_t>(j + 1);
                    assigned = true;
                    break;
                }
            }

            if (!assigned)
                res[i] = static_cast<uint64_t>(n);
        }

        return res;
    }

    /***********************************************************
     * Jenks natural breaks (core algorithm)
     ***********************************************************/
    inline std::vector<double> jenksBreaks(
        std::vector<double> data,
        size_t n_classes)
    {
        std::sort(data.begin(), data.end());

        size_t n = data.size();

        std::vector<std::vector<size_t>> lower(n + 1,
            std::vector<size_t>(n_classes + 1));

        std::vector<std::vector<double>> var(n + 1,
            std::vector<double>(n_classes + 1,
            std::numeric_limits<double>::infinity()));

        for (size_t i = 1; i <= n_classes; ++i)
        {
            lower[1][i] = 1;
            var[1][i] = 0.0;
        }

        for (size_t l = 2; l <= n; ++l)
        {
            double sum = 0.0;
            double sumsq = 0.0;

            for (size_t m = 1; m <= l; ++m)
            {
                size_t idx = l - m;
                double val = data[idx];

                sum += val;
                sumsq += val * val;

                double variance = sumsq - (sum * sum) / m;

                if (idx != 0)
                {
                    for (size_t j = 2; j <= n_classes; ++j)
                    {
                        double v = variance + var[idx][j - 1];

                        if (var[l][j] >= v)
                        {
                            lower[l][j] = idx + 1;
                            var[l][j] = v;
                        }
                    }
                }
            }

            lower[l][1] = 1;
            var[l][1] = sumsq - (sum * sum) / l;
        }

        std::vector<double> breaks(n_classes - 1);

        size_t k = n;

        for (size_t c = n_classes - 1; c > 0; --c)
        {
            size_t break_idx = lower[k][c + 1];
            if (break_idx < 2 || break_idx - 2 >= data.size()) {
                breaks[c - 1] = (c == n_classes - 1) ? data.back() : data.front();
            } else {
                breaks[c - 1] = data[break_idx - 2];
            }
            k = break_idx - 1;
        }

        return breaks;
    }

    /***********************************************************
     * Natural breaks discretization
     ***********************************************************/
    inline std::vector<uint64_t> naturalDisc(
        const std::vector<double>& vec,
        size_t n,
        size_t sample_begin = 3000,
        double sample_prob = 0.15,
        uint64_t seed = 123456789,
        bool right_closed = true)
    {
        auto x = remove_nan(vec);

        if (x.size() > sample_begin)
        {
            std::mt19937_64 rng(seed);

            size_t sample_size =
                std::round(x.size() * sample_prob);

            sample_size = std::max<size_t>(1, sample_size);

            std::vector<size_t> idx(x.size());
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);

            std::vector<double> data;
            data.reserve(sample_size);
            for (size_t i = 0; i < sample_size; ++i)
                data.push_back(x[idx[i]]);
            x = std::move(data);
        }

        auto breaks = jenksBreaks(x, n);

        std::vector<uint64_t> res(vec.size(), 0);

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;

            bool assigned = false;

            for (size_t j = 0; j < breaks.size(); ++j)
            {   
                bool classify_val = right_closed ? (vec[i] <= breaks[j])
                                                : (vec[i] < breaks[j]);
                if (classify_val)
                {
                    res[i] = static_cast<uint64_t>(j + 1);
                    assigned = true;
                    break;
                }
            }

            if (!assigned)
                res[i] = static_cast<uint64_t>(n);
        }

        return res;
    }

    /***********************************************************
     * Head / Tail breaks discretization
     ***********************************************************/
    inline std::vector<uint64_t> htDisc(
        const std::vector<double>& vec,
        double threshold = 0.4,
        size_t iter_step = 100,
        bool right_closed = true)
    {
        std::vector<double> x = remove_nan(vec);

        std::vector<uint64_t> result(vec.size(), 0);

        if (x.empty()) return result;

        std::vector<double> head = x;

        std::vector<double> breaks;
        breaks.push_back(min_val(x));

        for (size_t i = 0; i < iter_step; ++i)
        {
            double mu = mean(head);

            breaks.push_back(mu);

            std::vector<double> new_head;

            for (double v : head)
            {
                if (v > mu)
                {
                    new_head.push_back(v);
                }
            }

            if (new_head.empty()) break;

            double prop =
                static_cast<double>(new_head.size()) /
                static_cast<double>(head.size());

            head = new_head;

            if (prop >= threshold || head.size() <= 1) break;
        }

        breaks.push_back(max_val(x));

        std::sort(breaks.begin(), breaks.end());

        breaks.erase(
            std::unique(breaks.begin(), breaks.end()),
            breaks.end());

        if (breaks.size() < 2) 
        {
            for (size_t i = 0; i < vec.size(); ++i)
                if (!std::isnan(vec[i])) result[i] = 1;
            return result;
        }

        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (std::isnan(vec[i])) continue;

            bool assigned = false;

            for (size_t j = 0; j < breaks.size() - 1; ++j)
            {
                if (right_closed ? vec[i] <= breaks[j + 1]
                                : vec[i] < breaks[j + 1])
                {
                    result[i] = static_cast<uint64_t>(j + 1);
                    assigned = true;
                    break;
                }
            }

            if (!assigned)
                result[i] = static_cast<uint64_t>(breaks.size() - 1);
        }

        return result;
    }

    /***********************************************************
     * Unified interface
     ***********************************************************/
    inline std::vector<uint64_t> discretize(
        const std::vector<double>& vec,
        const std::string& method,
        size_t n = 5,
        size_t sample_begin = 3000,
        double sample_prob = 0.15,
        uint64_t seed = 123456789,
        double threshold = 0.4,
        size_t iter_step = 100,
        const std::vector<double>& breakpoints = {},
        bool right_closed = true)
    {
        if (method == "sd")
            return sdDisc(vec, n, right_closed);
        else if (method == "equal")
            return equalDisc(vec, n, right_closed);
        else if (method == "geometric")
            return geometricDisc(vec, n, right_closed);
        else if (method == "manual")
            return manualDisc(vec, breakpoints, right_closed);
        else if (method == "quantile")
            return quantileDisc(vec, n, right_closed);
        else if (method == "natural" || method == "jenks")
            return naturalDisc(vec, n, sample_begin, sample_prob, seed, right_closed);
        else if (method == "headtail"|| method == "headtails")  
            return htDisc(vec, threshold, iter_step, right_closed);
        else
            throw std::invalid_argument("Unknown discretization method");
    }

} // namespace discretize

}

#endif // INFOXTR_DISCRETIZE_HPP
