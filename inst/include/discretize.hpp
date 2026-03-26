/**********************************************************************
 *  File: discretize.hpp
 *
 *  Discretization utilities for continuous numeric vectors.
 *
 *  Methods implemented
 *      sd        Standard deviation discretization
 *      equal     Equal interval discretization
 *      geometric Geometric interval discretization
 *      quantile  Quantile discretization
 *      manual    Manual breakpoints
 *      natural   Jenks natural breaks
 *      headtail  Head/Tail breaks
 *
 *  Input
 *      const std::vector<double>& vec
 *
 *  Output
 *      std::vector<size_t>
 *          0  -> NaN values
 *          1..n -> discretized classes
 *
 *  Behaviour
 *      NaN values are ignored when computing breaks
 *      but assigned class 0 in the output.
 *
 *  RNG
 *      Sampling uses std::mt19937_64
 *
 **********************************************************************/

#ifndef DISCRETIZE_HPP
#define DISCRETIZE_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <iostream>
#include <stdexcept>

namespace Disc
{

/***********************************************************
 * Utility helpers
 ***********************************************************/

inline bool is_nan(double v)
{
    return std::isnan(v);
}

inline std::vector<double> remove_nan(const std::vector<double>& v, bool& has_nan)
{
    std::vector<double> out;
    out.reserve(v.size());

    for (double x : v)
    {
        if (std::isnan(x))
        {
            has_nan = true;
            continue;
        }
        out.push_back(x);
    }

    return out;
}

inline double mean(const std::vector<double>& v)
{
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

inline double stddev(const std::vector<double>& v)
{
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
inline std::vector<size_t> sdDisc(
    const std::vector<double>& vec,
    size_t n)
{
    bool has_nan = false;
    auto x = remove_nan(vec, has_nan);

    if (has_nan)
        std::cerr << "Warning: NaN values detected, assigned to class 0\n";

    double m = mean(x);
    double sd = stddev(x);

    std::vector<size_t> res(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (is_nan(vec[i]))
        {
            res[i] = 0;
            continue;
        }

        double diff = (vec[i] - m) / sd;
        long idx = std::floor((diff + n / 2.0) / n * n);

        idx = std::max<long>(1, std::min<long>(idx, n));

        res[i] = static_cast<size_t>(idx);
    }

    return res;
}

/***********************************************************
 * Equal interval discretization
 ***********************************************************/
inline std::vector<size_t> equalDisc(
    const std::vector<double>& vec,
    size_t n)
{
    bool has_nan = false;
    auto x = remove_nan(vec, has_nan);

    if (has_nan)
        std::cerr << "Warning: NaN values detected, assigned to class 0\n";

    double minx = min_val(x);
    double maxx = max_val(x);

    double interval = (maxx - minx) / n;

    std::vector<size_t> res(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (is_nan(vec[i]))
        {
            res[i] = 0;
            continue;
        }

        long idx = std::ceil((vec[i] - minx) / interval);
        idx = std::max<long>(1, std::min<long>(idx, n));

        res[i] = static_cast<size_t>(idx);
    }

    return res;
}

/***********************************************************
 * Geometric discretization
 ***********************************************************/
inline std::vector<size_t> geometricDisc(
    const std::vector<double>& vec,
    size_t n)
{
    bool has_nan = false;
    auto x = remove_nan(vec, has_nan);

    if (has_nan)
        std::cerr << "Warning: NaN values detected, assigned to class 0\n";

    double minx = min_val(x);
    double maxx = max_val(x);

    double factor = std::pow(maxx / minx, 1.0 / n);

    std::vector<size_t> res(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (is_nan(vec[i]))
        {
            res[i] = 0;
            continue;
        }

        long idx =
            std::floor(std::log(vec[i] / minx) / std::log(factor)) + 1;

        idx = std::max<long>(1, std::min<long>(idx, n));

        res[i] = static_cast<size_t>(idx);
    }

    return res;
}

/***********************************************************
 * Quantile discretization
 ***********************************************************/
inline std::vector<size_t> quantileDisc(
    const std::vector<double>& vec,
    size_t n)
{
    bool has_nan = false;
    auto x = remove_nan(vec, has_nan);

    if (has_nan)
        std::cerr << "Warning: NaN values detected, assigned to class 0\n";

    std::vector<double> sorted = x;
    std::sort(sorted.begin(), sorted.end());

    std::vector<double> q(n + 1);

    for (size_t i = 0; i <= n; ++i)
    {
        size_t idx = i * (sorted.size() - 1) / n;
        q[i] = sorted[idx];
    }

    std::vector<size_t> res(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (is_nan(vec[i]))
        {
            res[i] = 0;
            continue;
        }

        for (size_t j = 0; j < n; ++j)
        {
            if (vec[i] <= q[j + 1])
            {
                res[i] = j + 1;
                break;
            }
        }
    }

    return res;
}

/***********************************************************
 * Jenks natural breaks
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
        var[1][i] = 0;
    }

    for (size_t l = 2; l <= n; ++l)
    {
        double sum = 0;
        double sumsq = 0;

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
        breaks[c - 1] = data[lower[k][c + 1] - 2];
        k = lower[k][c + 1] - 1;
    }

    return breaks;
}

/***********************************************************
 * Natural breaks discretization
 ***********************************************************/
inline std::vector<size_t> naturalDisc(
    const std::vector<double>& vec,
    size_t n,
    double sample_prob = 0.15,
    uint64_t seed = 123456789)
{
    bool has_nan = false;
    auto x = remove_nan(vec, has_nan);

    if (has_nan)
        std::cerr << "Warning: NaN values detected, assigned to class 0\n";

    std::vector<double> data = x;

    if (x.size() > 3000)
    {
        std::mt19937_64 rng(seed);

        size_t sample_size =
            std::round(x.size() * sample_prob);

        sample_size = std::max<size_t>(1, sample_size);

        std::vector<size_t> idx(x.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::shuffle(idx.begin(), idx.end(), rng);

        data.resize(sample_size);

        for (size_t i = 0; i < sample_size; ++i)
            data[i] = x[idx[i]];
    }

    auto breaks = jenksBreaks(data, n);

    std::vector<size_t> res(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (is_nan(vec[i]))
        {
            res[i] = 0;
            continue;
        }

        bool assigned = false;

        for (size_t j = 0; j < breaks.size(); ++j)
        {
            if (vec[i] < breaks[j])
            {
                res[i] = j + 1;
                assigned = true;
                break;
            }
        }

        if (!assigned)
            res[i] = n;
    }

    return res;
}

/***********************************************************
 * Unified interface
 ***********************************************************/
inline std::vector<size_t> Disc(
    const std::vector<double>& vec,
    const std::string& method,
    size_t n = 5,
    double sample_prob = 1.0,
    uint64_t seed = 123456)
{
    if (method == "sd")
        return sdDisc(vec, n);

    if (method == "equal")
        return equalDisc(vec, n);

    if (method == "geometric")
        return geometricDisc(vec, n);

    if (method == "quantile")
        return quantileDisc(vec, n);

    if (method == "natural")
        return naturalDisc(vec, n, sample_prob, seed);

    throw std::invalid_argument("Unknown discretization method");
}

}

#endif
