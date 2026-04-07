/******************************************************************************
 * File: surd.hpp
 *
 * Synergistic–Unique–Redundant Decomposition (SURD) of mutual
 * information between a target variable and multiple source variables.
 *
 * The SURD framework decomposes the mutual information structure into:
 *
 *   0 = Redundant information
 *       Information shared among multiple variables.
 *
 *   1 = Unique information
 *       Information uniquely provided by a single variable.
 *
 *   2 = Synergistic information
 *       Information that only emerges from combinations of variables.
 *
 *   3 = Information loss
 *       Remaining uncertainty in the target after conditioning on all sources.
 *
 * The algorithm operates as follows:
 *
 *   1. Generate all combinations of source variables up to max_order.
 *   2. Compute mutual information I(Y ; X_set) for each combination.
 *   3. Group combinations by order (number of variables).
 *   4. Sort each group by mutual information.
 *   5. Perform ladder-style decomposition:
 *
 *        Order 1:
 *           Incremental differences determine redundant vs unique
 *           contributions.
 *
 *        Higher orders:
 *           Information exceeding the maximum lower-order mutual
 *           information is attributed to synergy.
 *
 *   6. Optionally normalize contributions so they sum to one.
 *   7. Append information loss as the final component.
 *
 * Two data types are supported:
 *
 *   Discrete data
 *      Mutual information computed using joint entropy estimators.
 *
 *   Continuous data
 *      Mutual information estimated via the KSG k-nearest neighbor method.
 *
 * Input data format:
 *
 *   Row 0  : target variable
 *   Row 1+ : source variables
 *
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ******************************************************************************/

#ifndef INFOXTR_SURD_HPP
#define INFOXTR_SURD_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infoxtr/combn.hpp"
#include "infoxtr/numericutils.hpp"
#include "infoxtr/infotheo.hpp"
#include <RcppThread.h>

namespace infoxtr 
{

namespace surd 
{

    using DiscMat = std::vector<std::vector<uint64_t>>;

    /***********************************************************
     * Result structure
     ***********************************************************/
    struct SURDRes
    {
        std::vector<std::vector<size_t>> unique_vars;
        std::vector<double>              unique_vals;

        std::vector<std::vector<size_t>> redundant_vars;
        std::vector<double>              redundant_vals;

        std::vector<std::vector<size_t>> synergy_vars;
        std::vector<double>              synergy_vals;

        std::vector<std::vector<size_t>> mi_vars;
        std::vector<double>              mi_vals;

        double info_leak = 0.0;
    };

    /***************************************************************
     * Synergistic-Unique-Redundant Decomposition for Discrete Data
     ***************************************************************/

inline SURDRes surd(
    const DiscMat& mat,
    size_t max_order = std::numeric_limits<size_t>::max(),
    size_t threads = 1,
    double base = 2.0,
    bool normalize = false)
{   
    if (threads == 0) threads = 1;
    size_t hw = std::thread::hardware_concurrency();
    if (hw > 0) threads = std::min(threads, hw);

    if (mat.size() < 2)
        throw std::invalid_argument("SURD requires >=2 variables");

    const size_t n_vars    = mat.size();
    const size_t n_sources = n_vars - 1;

    max_order = std::min(max_order , n_sources);

    const double log_base = std::log(base);

    /***********************************************************
     * Construct variable combination vector
     ***********************************************************/
    std::vector<size_t> ag_idx(n_sources);
    std::iota(ag_idx.begin(), ag_idx.end(), 1);

    const std::vector<std::vector<size_t>> combs =
        infoxtr::combn::genSubsets(ag_idx, max_order);

    const size_t n_combs = combs.size();
    
    // Initialize containers for result
    std::vector<double> info(n_combs , 0.0);
    std::vector<double> I_unique(n_sources , 0.0);
    std::vector<std::vector<size_t>> unique_vars(n_sources);

    for (size_t i = 0; i < n_sources; i++)
        unique_vars[i] = {i + 1};
    
    std::vector<double> I_R(n_combs , 0.0);
    std::vector<double> I_S(n_combs , 0.0);

    /***********************************************************
     * Joint table
     ***********************************************************/
    infoxtr::infotheo::JointTable jt =
        infoxtr::infotheo::joint_table(mat);

    const size_t k        = jt.k;
    const size_t n_states = jt.counts.size();

    size_t total_n = 0;
    for (auto c : jt.counts)
        total_n += c;

    /***********************************************************
     * Target states
     ***********************************************************/
    std::vector<uint64_t> s_vals;
    std::vector<size_t>   s_counts;

    for (size_t i = 0; i < n_states; i++)
    {
        uint64_t s = jt.states[i * k];

        auto it = std::find(s_vals.begin(), s_vals.end(), s);

        if (it == s_vals.end())
        {
            s_vals.push_back(s);
            s_counts.push_back(jt.counts[i]);
        }
        else
        {
            size_t idx = it - s_vals.begin();
            s_counts[idx] += jt.counts[i];
        }
    }

    const size_t n_s = s_vals.size();

    /***********************************************************
     * state -> s index
     ***********************************************************/
    std::vector<size_t> state_s_index(n_states);

    for (size_t i = 0; i < n_states; i++)
    {
        uint64_t s = jt.states[i * k];

        for (size_t j = 0; j < n_s; j++)
            if (s_vals[j] == s)
                state_s_index[i] = j;
    }

    /***********************************************************
     * Precompute projection + px groups
     ***********************************************************/
    struct PxGroup
    {
        std::vector<uint64_t> state;
        size_t count = 0;
        std::vector<size_t> rows;
    };

    std::vector<std::vector<PxGroup>> px_groups(n_combs);

    if (threads <= 1) 
    {   
        for (size_t ci = 0; ci < n_combs; ci++)
        {
            auto & subset = combs[ci];
            const size_t kk = subset.size();

            std::vector<std::vector<uint64_t>> proj(n_states);

            for (size_t i = 0; i < n_states; i++)
            {
                size_t base = i * k;

                proj[i].resize(kk);

                for (size_t j = 0; j < kk; j++)
                    proj[i][j] = jt.states[base + subset[j]];
            }

            std::vector<size_t> order(n_states);

            for (size_t i = 0; i < n_states; i++)
                order[i] = i;

            std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b)
            {
                for (size_t j = 0; j < kk; j++)
                {
                    if (proj[a][j] < proj[b][j]) return true;
                    if (proj[a][j] > proj[b][j]) return false;
                }
                return false;
            });

            PxGroup g;
            g.state = proj[order[0]];
            g.count = jt.counts[order[0]];
            g.rows.push_back(order[0]);

            for (size_t i = 1; i < n_states; i++)
            {
                bool same = true;

                for (size_t j = 0; j < kk; j++)
                    if (proj[order[i]][j] != g.state[j])
                        same = false;

                if (same)
                {
                    g.count += jt.counts[order[i]];
                    g.rows.push_back(order[i]);
                }
                else
                {
                    px_groups[ci].push_back(g);

                    g.state = proj[order[i]];
                    g.count = jt.counts[order[i]];
                    g.rows.clear();
                    g.rows.push_back(order[i]);
                }
            }

            px_groups[ci].push_back(g);
        }
    } 
    else  
    {
        RcppThread::parallelFor(0, n_combs, [&](size_t ci) {
            auto & subset = combs[ci];
            const size_t kk = subset.size();

            std::vector<std::vector<uint64_t>> proj(n_states);

            for (size_t i = 0; i < n_states; i++)
            {
                size_t base = i * k;

                proj[i].resize(kk);

                for (size_t j = 0; j < kk; j++)
                    proj[i][j] = jt.states[base + subset[j]];
            }

            std::vector<size_t> order(n_states);

            for (size_t i = 0; i < n_states; i++)
                order[i] = i;

            std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b)
            {
                for (size_t j = 0; j < kk; j++)
                {
                    if (proj[a][j] < proj[b][j]) return true;
                    if (proj[a][j] > proj[b][j]) return false;
                }
                return false;
            });

            PxGroup g;
            g.state = proj[order[0]];
            g.count = jt.counts[order[0]];
            g.rows.push_back(order[0]);

            for (size_t i = 1; i < n_states; i++)
            {
                bool same = true;

                for (size_t j = 0; j < kk; j++)
                    if (proj[order[i]][j] != g.state[j])
                        same = false;

                if (same)
                {
                    g.count += jt.counts[order[i]];
                    g.rows.push_back(order[i]);
                }
                else
                {
                    px_groups[ci].push_back(g);

                    g.state = proj[order[i]];
                    g.count = jt.counts[order[i]];
                    g.rows.clear();
                    g.rows.push_back(order[i]);
                }
            }

            px_groups[ci].push_back(g);
        }, threads);
    } 

    /***********************************************************
     * Loop target states
     ***********************************************************/
    for (size_t si = 0; si < n_s; si++)
    {
        double p_s =
            static_cast<double>(s_counts[si]) / total_n;

        std::vector<double> I_s(n_combs , 0.0);
        
        if (threads <= 1)
        {
            for (size_t ci = 0; ci < n_combs; ci++)
            {
                double sum = 0.0;

                for (auto & g : px_groups[ci])
                {
                    size_t psx_count = 0;

                    for (auto r : g.rows)
                        if (state_s_index[r] == si)
                            psx_count += jt.counts[r];

                    if (psx_count == 0)
                        continue;

                    double psx =
                        static_cast<double>(psx_count) / total_n;

                    double px =
                        static_cast<double>(g.count) / total_n;

                    sum += psx * std::log(psx / (p_s * px));
                }

                double pointwise = sum / p_s / log_base;

                I_s[ci] = pointwise;

                // accumulate MI 
                info[ci] += p_s * pointwise;
            }
        }
        else 
        {
            RcppThread::parallelFor(0, n_combs, [&](size_t ci) {
                double sum = 0.0;

                for (auto & g : px_groups[ci])
                {
                    size_t psx_count = 0;

                    for (auto r : g.rows)
                        if (state_s_index[r] == si)
                            psx_count += jt.counts[r];

                    if (psx_count == 0)
                        continue;

                    double psx =
                        static_cast<double>(psx_count) / total_n;

                    double px =
                        static_cast<double>(g.count) / total_n;

                    sum += psx * std::log(psx / (p_s * px));
                }

                double pointwise = sum / p_s / log_base;

                I_s[ci] = pointwise;

                // accumulate MI 
                info[ci] += p_s * pointwise;
            }, threads);
        }

        /**************************************************
         * SURD decomposition
         **************************************************/

        struct Node
        {
            size_t idx;
            size_t len;
            double val;
        };

        std::vector<Node> nodes;
        nodes.reserve(n_combs);

        for (size_t i = 0; i < n_combs; i++)
        {
            Node n;
            n.idx = i;
            n.len = combs[i].size();
            n.val = I_s[i];
            nodes.push_back(n);
        }

        // sort once 
        std::sort(nodes.begin(), nodes.end(),
        [](const Node & a, const Node & b)
        {
            return a.val < b.val;
        });

        // find max subset length 
        size_t max_len = 0;

        for (auto & n : nodes)
            if (n.len > max_len)
                max_len = n.len;

        /**************************************************
         * SURD monotonic filter
         **************************************************/

        for (size_t l = 1; l < max_len; l++)
        {
            double Il1max = -std::numeric_limits<double>::infinity();

            for (auto & n : nodes)
                if (n.len == l)
                    if (n.val > Il1max)
                        Il1max = n.val;

            for (auto & n : nodes)
                if (n.len == l + 1)
                    if (n.val < Il1max)
                        n.val = 0.0;
        }

        /**************************************************
         * SURD information layers
         **************************************************/

        std::vector<size_t> red_vars(n_sources);
        for(size_t i = 0; i < n_sources; i++)
            red_vars[i] = i+1;

        double prev = 0.0;

        for (auto & n : nodes)
        {
            double delta = n.val - prev;

            if (delta < 0)
                delta = 0;

            double info_add = delta * p_s;

            const auto & subset = combs[n.idx];

            if (subset.size() == 1)
            {
                // redundant information
                for (size_t ri = 0; ri < n_combs; ri++)
                {
                    if (combs[ri] == red_vars)
                    {
                        I_R[ri] += info_add;
                        break;
                    }
                }

                auto it = std::find(
                    red_vars.begin(),
                    red_vars.end(),
                    subset[0]);

                if (it != red_vars.end())
                    red_vars.erase(it);
            }
            else
            {
                // synergy information
                I_S[n.idx] += info_add;
            }

            if (n.val > prev)
                prev = n.val;
        }
    }

    SURDRes result;

    /**************************************************
     * Redundant + Unique
     **************************************************/

    for (size_t i = 0; i < n_combs; i++)
    {
        if (I_R[i] > 0)
        {
            if (combs[i].size() == 1)
            {
                result.unique_vars.push_back(combs[i]);
                result.unique_vals.push_back(I_R[i]);
            }
            else
            {
                result.redundant_vars.push_back(combs[i]);
                result.redundant_vals.push_back(I_R[i]);
            }
        }
    }

    /**************************************************
     * Synergy
     **************************************************/

    for (size_t i = 0; i < n_combs; i++)
    {
        if (combs[i].size() > 1 && I_S[i] > 0)
        {
            result.synergy_vars.push_back(combs[i]);
            result.synergy_vals.push_back(I_S[i]);
        }
    }

    for (size_t i = 0; i < n_combs; i++)
    {
        result.mi_vars.push_back(combs[i]);
        result.mi_vals.push_back(info[i]);
    }

    // Information leak
    std::vector<size_t> target_idx = {0};
    double leak = infoxtr::infotheo::ce(mat, target_idx, ag_idx, base, false) 
                  / infoxtr::infotheo::je(mat, target_idx, base, false);
    result.info_leak = std::max(0.0, std::min(1.0, leak));

    /***********************************************************
     * Normalize (optional)
     ***********************************************************/
    if (normalize)
    {
        double max_mi = 0.0;

        for (double v : result.mi_vals)
            if (v > max_mi)
                max_mi = v;

        if (max_mi <= 0.0)
            max_mi = 1.0;

        for (auto & v : result.unique_vals)
            v /= max_mi;

        for (auto & v : result.redundant_vals)
            v /= max_mi;

        for (auto & v : result.synergy_vals)
            v /= max_mi;

        // for (auto & v : result.mi_vals)
        //     v /= max_mi;
    }

    return result;
}

// inline SURDRes surd_pointwise_fast(
//     const DiscMat& mat,
//     size_t max_order = std::numeric_limits<size_t>::max(),
//     size_t threads = 1,
//     double base = 2.0,
//     bool normalize = false)
// {
//     SURDRes result;

//     if (mat.size() < 2)
//         throw std::invalid_argument("SURD requires >=2 variables");

//     const size_t n_vars    = mat.size();
//     const size_t n_sources = n_vars - 1;

//     max_order = std::min(max_order, n_sources);

//     const double log_base = std::log(base);

//     /***********************************************************
//      * Generate subsets
//      ***********************************************************/
//     std::vector<std::vector<size_t>> combs;

//     for (size_t mask = 1; mask < (1ULL << n_sources); mask++)
//     {
//         std::vector<size_t> subset;

//         for (size_t j = 0; j < n_sources; j++)
//             if (mask & (1ULL << j))
//                 subset.push_back(j + 1);

//         if (subset.size() <= max_order)
//             combs.push_back(subset);
//     }

//     const size_t n_combs = combs.size();

//     std::vector<double> info(n_combs, 0.0);

//     /*********************************************************************
//      * Joint table
//      *********************************************************************/
//     infoxtr::infotheo::JointTable jt = infoxtr::infotheo::joint_table(mat);

//     const size_t k        = jt.k;
//     const size_t n_states = jt.counts.size();

//     size_t total_n = 0;
//     for (auto c : jt.counts)
//         total_n += c;

//     /***********************************************************
//      * Target states
//      ***********************************************************/
//     std::vector<uint64_t> s_vals;
//     std::vector<size_t>   s_counts;

//     for (size_t i = 0; i < n_states; i++)
//     {
//         uint64_t s = jt.states[i * k];

//         auto it = std::find(s_vals.begin(), s_vals.end(), s);

//         if (it == s_vals.end())
//         {
//             s_vals.push_back(s);
//             s_counts.push_back(jt.counts[i]);
//         }
//         else
//         {
//             size_t idx = it - s_vals.begin();
//             s_counts[idx] += jt.counts[i];
//         }
//     }

//     const size_t n_s = s_vals.size();

//     /***********************************************************
//      * state -> s index
//      ***********************************************************/
//     std::vector<size_t> state_s_index(n_states);

//     for (size_t i = 0; i < n_states; i++)
//     {
//         uint64_t s = jt.states[i * k];

//         for (size_t j = 0; j < n_s; j++)
//             if (s_vals[j] == s)
//                 state_s_index[i] = j;
//     }

//     /***********************************************************
//      * Precompute projections and px
//      ***********************************************************/
//     std::vector<std::vector<uint64_t>> proj_all(n_combs);
//     std::vector<std::vector<uint64_t>> px_states_all(n_combs);
//     std::vector<std::vector<size_t>>   px_counts_all(n_combs);

//     for (size_t ci = 0; ci < n_combs; ci++)
//     {
//         auto & subset = combs[ci];

//         const size_t kk = subset.size();

//         auto & proj = proj_all[ci];
//         proj.reserve(n_states * kk);

//         for (size_t i = 0; i < n_states; i++)
//         {
//             size_t base = i * k;

//             for (auto v : subset)
//                 proj.push_back(jt.states[base + v]);
//         }

//         std::vector<size_t> order(n_states);

//         for (size_t i = 0; i < n_states; i++)
//             order[i] = i;

//         std::sort(order.begin(), order.end(),
//         [&](size_t a, size_t b)
//         {
//             size_t ia = a * kk;
//             size_t ib = b * kk;

//             for (size_t j = 0; j < kk; j++)
//             {
//                 uint64_t va = proj[ia + j];
//                 uint64_t vb = proj[ib + j];

//                 if (va < vb) return true;
//                 if (va > vb) return false;
//             }

//             return false;
//         });

//         auto equal = [&](size_t a, size_t b)
//         {
//             size_t ia = a * kk;
//             size_t ib = b * kk;

//             for (size_t j = 0; j < kk; j++)
//                 if (proj[ia + j] != proj[ib + j])
//                     return false;

//             return true;
//         };

//         auto & px_states = px_states_all[ci];
//         auto & px_counts = px_counts_all[ci];

//         size_t run = jt.counts[order[0]];

//         for (size_t i = 1; i < n_states; i++)
//         {
//             if (equal(order[i], order[i - 1]))
//             {
//                 run += jt.counts[order[i]];
//             }
//             else
//             {
//                 size_t idx = order[i - 1] * kk;

//                 for (size_t j = 0; j < kk; j++)
//                     px_states.push_back(proj[idx + j]);

//                 px_counts.push_back(run);

//                 run = jt.counts[order[i]];
//             }
//         }

//         size_t idx_last = order.back() * kk;

//         for (size_t j = 0; j < kk; j++)
//             px_states.push_back(proj[idx_last + j]);

//         px_counts.push_back(run);
//     }

//     /***********************************************************
//      * Loop target states
//      ***********************************************************/
//     for (size_t si = 0; si < n_s; si++)
//     {
//         const double p_s =
//             static_cast<double>(s_counts[si]) / total_n;

//         std::vector<double> I_s(n_combs, 0.0);

//         for (size_t ci = 0; ci < n_combs; ci++)
//         {
//             const auto & subset = combs[ci];
//             const size_t kk = subset.size();

//             const auto & proj = proj_all[ci];
//             const auto & px_states = px_states_all[ci];
//             const auto & px_counts = px_counts_all[ci];

//             std::vector<uint64_t> proj_s;
//             std::vector<size_t>   weight_s;

//             for (size_t i = 0; i < n_states; i++)
//             {
//                 if (state_s_index[i] != si)
//                     continue;

//                 size_t base = i * kk;

//                 for (size_t j = 0; j < kk; j++)
//                     proj_s.push_back(proj[base + j]);

//                 weight_s.push_back(jt.counts[i]);
//             }

//             if (proj_s.empty())
//                 continue;

//             const size_t ns = weight_s.size();

//             std::vector<size_t> order(ns);

//             for (size_t i = 0; i < ns; i++)
//                 order[i] = i;

//             std::sort(order.begin(), order.end(),
//             [&](size_t a, size_t b)
//             {
//                 size_t ia = a * kk;
//                 size_t ib = b * kk;

//                 for (size_t j = 0; j < kk; j++)
//                 {
//                     uint64_t va = proj_s[ia + j];
//                     uint64_t vb = proj_s[ib + j];

//                     if (va < vb) return true;
//                     if (va > vb) return false;
//                 }

//                 return false;
//             });

//             auto equal = [&](size_t a, size_t b)
//             {
//                 size_t ia = a * kk;
//                 size_t ib = b * kk;

//                 for (size_t j = 0; j < kk; j++)
//                     if (proj_s[ia + j] != proj_s[ib + j])
//                         return false;

//                 return true;
//             };

//             std::vector<uint64_t> psx_states;
//             std::vector<size_t>   psx_counts;

//             size_t run = weight_s[order[0]];

//             for (size_t i = 1; i < ns; i++)
//             {
//                 if (equal(order[i], order[i - 1]))
//                 {
//                     run += weight_s[order[i]];
//                 }
//                 else
//                 {
//                     size_t idx = order[i - 1] * kk;

//                     for (size_t j = 0; j < kk; j++)
//                         psx_states.push_back(proj_s[idx + j]);

//                     psx_counts.push_back(run);

//                     run = weight_s[order[i]];
//                 }
//             }

//             size_t idx = order.back() * kk;

//             for (size_t j = 0; j < kk; j++)
//                 psx_states.push_back(proj_s[idx + j]);

//             psx_counts.push_back(run);

//             double sum = 0.0;

//             for (size_t i = 0; i < psx_counts.size(); i++)
//             {
//                 double psx =
//                     static_cast<double>(psx_counts[i]) / total_n;

//                 size_t base = i * kk;

//                 for (size_t j = 0; j < px_counts.size(); j++)
//                 {
//                     bool match = true;

//                     for (size_t t = 0; t < kk; t++)
//                         if (px_states[j * kk + t] != psx_states[base + t])
//                             match = false;

//                     if (match)
//                     {
//                         double px =
//                             static_cast<double>(px_counts[j]) / total_n;

//                         sum += psx *
//                                std::log(psx / (p_s * px));

//                         break;
//                     }
//                 }
//             }

//             double pointwise = sum / p_s / log_base;

//             I_s[ci] = pointwise;
//         }

//         /**************************************************
//          * SURD decomposition (Python aligned)
//          **************************************************/

//         struct Node
//         {
//             size_t idx;
//             size_t len;
//             double val;
//         };

//         std::vector<Node> nodes;
//         nodes.reserve(n_combs);

//         for (size_t i = 0; i < n_combs; i++)
//         {
//             Node nd;
//             nd.idx = i;
//             nd.len = combs[i].size();
//             nd.val = I_s[i];
//             nodes.push_back(nd);
//         }

//         /* sort by information */
//         std::sort(nodes.begin(), nodes.end(),
//         [](const Node & a, const Node & b)
//         {
//             return a.val < b.val;
//         });

//         /* find max subset size */
//         size_t max_len = 0;
//         for (auto & n : nodes)
//             if (n.len > max_len)
//                 max_len = n.len;

//         /* SURD monotonicity filter */
//         for (size_t l = 1; l < max_len; l++)
//         {
//             double Il1max = -1e300;

//             for (auto & n : nodes)
//                 if (n.len == l)
//                     if (n.val > Il1max)
//                         Il1max = n.val;

//             for (auto & n : nodes)
//                 if (n.len == l + 1)
//                     if (n.val < Il1max)
//                         n.val = 0.0;
//         }

//         /* re-sort (python does argsort again) */
//         std::sort(nodes.begin(), nodes.end(),
//         [](const Node & a, const Node & b)
//         {
//             return a.val < b.val;
//         });

//         /* diff + distribute */
//         double prev = 0.0;

//         for (auto & n : nodes)
//         {
//             double delta = n.val - prev;

//             if (delta < 0)
//                 delta = 0;

//             info[n.idx] += delta * p_s;

//             if (n.val > prev)
//                 prev = n.val;
//         }
//     }

//     /***********************************************************
//      * Store results
//      ***********************************************************/
//     for (size_t i = 0; i < n_combs; i++)
//     {
//         result.values.push_back(info[i]);

//         size_t kk = combs[i].size();

//         if (kk == 1) result.types.push_back(0);
//         else if (kk == 2) result.types.push_back(1);
//         else result.types.push_back(2);

//         result.var_indices.push_back(combs[i]);
//     }

//     return result;
// }
    
} // namespace surd

}

#endif // INFOXTR_SURD_HPP
