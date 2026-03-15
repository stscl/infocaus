/**********************************************************************
 *  File: infotheo.hpp
 *
 *  High performance information theoretic measurements
 *  for discrete state data (uint64_t encoded).
 *
 *  Data layout:
 *      Series = std::vector<uint64_t>
 *      Matrix = std::vector<std::vector<uint64_t>>   // mat[var][obs]
 *
 *  Functions:
 *      Entropy
 *      JE   Joint Entropy
 *      CE   Conditional Entropy
 *      MI   Mutual Information
 *      CMI  Conditional Mutual Information
 *
 *  Author: Wenbo Lyu
 *  License: GPL-3
 **********************************************************************/

#ifndef INFOTHEO_HPP
#define INFOTHEO_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

namespace InfoTheo
{

using Series = std::vector<uint64_t>;
using Matrix = std::vector<std::vector<uint64_t>>;

/***********************************************************
 * Utilities
 ***********************************************************/
inline double convert_log_base(double x, double base)
{
    if (x <= 0.0) return 0.0;

    if (!(base > 0.0) || std::abs(base - 1.0) < 1e-12)
        throw std::invalid_argument("Log base must be positive and not equal to 1.");

    return x / std::log(base);
}

/***********************************************************
 * Strong 64bit mix
 ***********************************************************/
static inline uint64_t mix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/***********************************************************
 * Combine multiple states into one hash
 ***********************************************************/
inline uint64_t hash_combine(uint64_t h, uint64_t v)
{
    return h ^ (mix64(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
}

/***********************************************************
 * Entropy
 ***********************************************************/
inline double Entropy(
    const Series& series,
    double base = 2.0,
    bool na_rm = true)
{
    if (series.empty())
        return std::numeric_limits<double>::quiet_NaN();

    std::unordered_map<uint64_t, size_t> freq;
    freq.reserve(series.size() * 1.3);

    size_t n_valid = 0;

    for (uint64_t x : series)
    {
        if (na_rm && x == 0)
            continue;

        ++freq[x];
        ++n_valid;
    }

    if (n_valid == 0)
        return std::numeric_limits<double>::quiet_NaN();

    double h = 0.0;

    for (const auto& kv : freq)
    {
        double p = static_cast<double>(kv.second) / n_valid;
        h -= p * std::log(p);
    }

    return convert_log_base(h, base);
}

/***********************************************************
 * Joint Entropy
 ***********************************************************/
inline double JE(
    const Matrix& mat,
    const std::vector<size_t>& vars,
    double base = 2.0,
    bool na_rm = true)
{
    if (mat.empty() || vars.empty())
        return std::numeric_limits<double>::quiet_NaN();

    const size_t n_obs = mat[0].size();
    const size_t n_cols = mat.size();

    std::unordered_set<size_t> valid_vars;
    valid_vars.reserve(vars.size());

    for (size_t v : vars)
        if (v < n_cols)
            valid_vars.insert(v);

    std::vector<size_t> clean_vars(valid_vars.begin(), valid_vars.end());

    std::unordered_map<uint64_t, size_t> freq;
    freq.reserve(n_obs * 1.3);

    size_t n_valid = 0;

    for (size_t i = 0; i < n_obs; ++i)
    {
        uint64_t key = 0;
        bool skip = false;

        for (size_t v : clean_vars)
        {
            uint64_t val = mat[v][i];

            if (na_rm && val == 0)
            {
                skip = true;
                break;
            }

            key = hash_combine(key, val);
        }

        if (skip) continue;

        ++freq[key];
        ++n_valid;
    }

    if (n_valid == 0)
        return std::numeric_limits<double>::quiet_NaN();

    double h = 0.0;

    for (const auto& kv : freq)
    {
        double p = static_cast<double>(kv.second) / n_valid;
        h -= p * std::log(p);
    }

    return convert_log_base(h, base);
}

/***********************************************************
 * Conditional Entropy
 ***********************************************************/
inline double CE(
    const Matrix& mat,
    const std::vector<size_t>& target,
    const std::vector<size_t>& conds,
    double base = 2.0,
    bool na_rm = true)
{
    if (mat.empty() || target.empty())
        return std::numeric_limits<double>::quiet_NaN();

    if (conds.empty())
        throw std::invalid_argument("conds cannot be empty");

    std::vector<size_t> tc = conds;
    tc.insert(tc.end(), target.begin(), target.end());

    return JE(mat, tc, base, na_rm)
         - JE(mat, conds, base, na_rm);
}

/***********************************************************
 * Mutual Information
 ***********************************************************/
inline double MI(
    const Matrix& mat,
    const std::vector<size_t>& target,
    const std::vector<size_t>& interact,
    double base = 2.0,
    bool na_rm = true)
{
    if (mat.empty() || target.empty())
        return std::numeric_limits<double>::quiet_NaN();

    if (interact.empty())
        throw std::invalid_argument("interact cannot be empty");

    std::vector<size_t> ti = interact;
    ti.insert(ti.end(), target.begin(), target.end());

    return JE(mat, target, base, na_rm) +
           JE(mat, interact, base, na_rm) -
           JE(mat, ti, base, na_rm);
}

/***********************************************************
 * Conditional Mutual Information
 ***********************************************************/
inline double CMI(
    const Matrix& mat,
    const std::vector<size_t>& target,
    const std::vector<size_t>& interact,
    const std::vector<size_t>& conds,
    double base = 2.0,
    bool na_rm = true)
{
    if (mat.empty() || target.empty())
        return std::numeric_limits<double>::quiet_NaN();

    if (interact.empty() || conds.empty())
        throw std::invalid_argument("interact / conds cannot be empty");

    std::vector<size_t> ct  = conds;
    std::vector<size_t> ci  = conds;
    std::vector<size_t> cti = conds;

    ct.insert(ct.end(), target.begin(), target.end());
    ci.insert(ci.end(), interact.begin(), interact.end());
    cti.insert(cti.end(), target.begin(), target.end());
    cti.insert(cti.end(), interact.begin(), interact.end());

    return JE(mat, ct, base, na_rm)
         + JE(mat, ci, base, na_rm)
         - JE(mat, conds, base, na_rm)
         - JE(mat, cti, base, na_rm);
}

} // namespace InfoTheo

#endif // INFOTHEO_HPP
