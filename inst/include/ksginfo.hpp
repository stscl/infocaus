/**********************************************************************
 *  File: ksginfo.hpp
 *
 *  Continuous information theoretic measurements
 *  using k-nearest neighbor estimators.
 *
 *  Algorithms:
 *
 *      Kozachenko–Leonenko entropy estimator
 *      Kraskov–Stögbauer–Grassberger MI estimator
 *
 *  Data layout:
 *      Series = std::vector<double>
 *      Matrix = std::vector<std::vector<double>> // mat[var][obs]
 *
 *  Estimator variants (alg parameter):
 *
 *      Entropy / Joint Entropy / Conditional Entropy
 *      ---------------------------------------------
 *      alg = 0
 *            Kozachenko–Leonenko entropy estimator
 *
 *      alg = 1
 *            Bias corrected KL estimator
 *            (adds +1/k correction)
 *
 *      Mutual Information / Conditional Mutual Information
 *      ---------------------------------------------------
 *      alg = 0
 *            Kraskov–Stögbauer–Grassberger estimator I (KSG1)
 *
 *      alg = 1
 *            Kraskov–Stögbauer–Grassberger estimator II (KSG2)
 *
 *  Distance backend:
 *      Dist::Dist (Chebyshev metric)
 *
 *  Special functions:
 *      NumericUtils::Digamma
 *
 *  Author: Wenbo Lyu (Github: @SpatLyu)
 *  License: GPL-3
 **********************************************************************/

#ifndef KSGINFO_HPP
#define KSGINFO_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "distance.hpp"
#include "numericutils.hpp"

namespace KSGInfo
{

    using Series = std::vector<double>;
    using Matrix = std::vector<std::vector<double>>;

    /***********************************************************
     * Utility: subset matrix
     ***********************************************************/
    inline Matrix subset(
        const Matrix& mat,
        const std::vector<size_t>& vars)
    {
        Matrix out;
        out.reserve(vars.size());

        for (size_t v : vars)
            out.push_back(mat[v]);

        return out;
    }

    /***********************************************************
     * Entropy (Kozachenko–Leonenko)
     ***********************************************************/
    inline double Entropy(
        const Series& series,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0)
    {
        const size_t n = series.size();

        Matrix vec;
        vec.reserve(1);
        vec.emplace_back(series);

        auto dist = Dist::Dist(vec,"maximum",true,false);

        double avg = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            std::vector<double> row = dist[i];

            if (i < row.size())
                row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbour count");

            std::nth_element(
                row.begin(),
                row.begin() + (k - 1),
                row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);
            
            avg += (NumericUtils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                    ? 0.0 : std::log(eps * 2.0);
        }

        avg /= static_cast<double>(n);

        double H = NumericUtils::Digamma(n)
                 - NumericUtils::Digamma(k)
                 + avg;
        
        if (alg == 1)
            H += 1.0 / k;

        if (!NumericUtils::doubleNearlyEqual(base,std::exp(1.0)))
            H /= std::log(base);

        return H;
    }

    /***********************************************************
     * Joint Entropy
     ***********************************************************/
    inline double JE(
        const Matrix& mat,
        const std::vector<size_t>& vars,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0)
    {
        Matrix sub = subset(mat,vars);

        const size_t d = sub.size();
        const size_t n = sub[0].size();

        auto dist = Dist::Dist(sub,"maximum",true,false);

        double avg = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            std::vector<double> row = dist[i];

            if (i < row.size())
                row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbour count");

            std::nth_element(
                row.begin(),
                row.begin() + (k - 1),
                row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            avg += (NumericUtils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                    ? 0.0 : std::log(eps * 2.0);
        }

        avg /= static_cast<double>(n);

        double H = NumericUtils::Digamma(n)
                 - NumericUtils::Digamma(k)
                 + d * avg;

        if (alg == 1)
            H += 1.0 / k;

        if (!NumericUtils::doubleNearlyEqual(base,std::exp(1.0)))
            H /= std::log(base);

        return H;
    }

    /***********************************************************
     * Conditional Entropy
     ***********************************************************/
    inline double CE(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& cond,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0)
    {
        std::vector<size_t> tc = cond;
        tc.insert(tc.end(),target.begin(),target.end());

        return JE(mat,tc,k,alg,base) - JE(mat,cond,k,alg,base);
    }

    /***********************************************************
     * Mutual Information (KSG estimator)
     ***********************************************************/
    inline double MI(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& interact,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0,
        bool normalize = false)
    {
        std::vector<size_t> xy = target;
        xy.insert(xy.end(), interact.begin(), interact.end());

        auto d_xy = Dist::Dist(subset(mat,xy),"maximum",true,false);
        auto d_x  = Dist::Dist(subset(mat,target),"maximum",true,false);
        auto d_y  = Dist::Dist(subset(mat,interact),"maximum",true,false);

        const size_t n = d_xy.size();
        const size_t d = xy.size();

        double sum = 0.0;
        double avg_log_eps = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            std::vector<double> row = d_xy[i];

            if (i < row.size())
                row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbour count");

            std::nth_element(row.begin(),row.begin()+k-1,row.end());
            
            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            avg_log_eps += (NumericUtils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                            ? 0.0 : std::log(eps * 2.0);

            size_t nx = 0, ny = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i == j) continue;

                if (alg == 0)
                {
                    if (!std::isnan(d_x[i][j]) && d_x[i][j] < eps) nx++;
                    if (!std::isnan(d_y[i][j]) && d_y[i][j] < eps) ny++;
                }
                else 
                {
                    if (!std::isnan(d_x[i][j]) && d_x[i][j] <= eps) nx++;
                    if (!std::isnan(d_y[i][j]) && d_y[i][j] <= eps) ny++;
                } 
            }

            if (alg == 0)
                sum += NumericUtils::Digamma(nx+1)
                     + NumericUtils::Digamma(ny+1);
            else
                sum += NumericUtils::Digamma(nx)
                     + NumericUtils::Digamma(ny);
        }

        avg_log_eps /= n;

        double mi;

        if (alg == 0)
            mi = NumericUtils::Digamma(k)
               + NumericUtils::Digamma(n)
               - sum / n;
        else
            mi = NumericUtils::Digamma(k)
               - 1.0 / k
               + NumericUtils::Digamma(n)
               - sum / n;

        mi = std::max(0.0, mi);

        if (!normalize) 
        {
            if (!NumericUtils::doubleNearlyEqual(base,std::exp(1.0)))
                mi /= std::log(base);

            return mi;
        } 

        double hxy = NumericUtils::Digamma(n)
                   - NumericUtils::Digamma(k)
                   + d * avg_log_eps;
        if (alg == 1) hxy += 1.0 / k;

        if (hxy <= 0) {
            if (!NumericUtils::doubleNearlyEqual(base,std::exp(1.0)))
                mi /= std::log(base);

            return mi;
        } 

        return mi / hxy;
    }

    /***********************************************************
     * Conditional Mutual Information
     ***********************************************************/
    inline double CMI(
        const Matrix& mat,
        const std::vector<size_t>& target,
        const std::vector<size_t>& interact,
        const std::vector<size_t>& conds,
        size_t k = 3,
        size_t alg = 0,
        double base = 2.0,
        bool normalize = false)
    {
        std::vector<size_t> xyz = conds;
        xyz.insert(xyz.end(), target.begin(), target.end());
        xyz.insert(xyz.end(), interact.begin(), interact.end());

        std::vector<size_t> xy = target;
        xy.insert(xy.end(), interact.begin(), interact.end());

        std::vector<size_t> xz = conds;
        xz.insert(xz.end(), target.begin(), target.end());

        std::vector<size_t> yz = conds;
        yz.insert(yz.end(), interact.begin(), interact.end());

        auto d_xyz = Dist::Dist(subset(mat,xyz),"maximum",true,false);
        auto d_xz  = Dist::Dist(subset(mat,xz),"maximum",true,false);
        auto d_yz  = Dist::Dist(subset(mat,yz),"maximum",true,false);
        auto d_z   = Dist::Dist(subset(mat,conds),"maximum",true,false);

        const size_t n = d_xyz.size();
        const size_t d = xy.size();

        double sum = 0.0;
        double avg_log_eps = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            std::vector<double> row = d_xyz[i];

            if (i < row.size())
                row[i] = std::numeric_limits<double>::quiet_NaN();

            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](double v){ return std::isnan(v); }),
                row.end());

            if (row.size() < k)
                throw std::runtime_error("k larger than valid neighbour count");

            std::nth_element(row.begin(),row.begin()+k-1,row.end());

            double eps = row[k-1];
            // double eps = std::max(row[k-1], 1e-15);

            avg_log_eps += (NumericUtils::doubleNearlyEqual(eps*2.0, 0.0) || eps < 0) 
                            ? 0.0 : std::log(eps * 2.0);

            size_t nxz = 0, nyz = 0, nz = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i == j) continue;

                if (alg == 0)
                {
                    if (!std::isnan(d_xz[i][j]) && d_xz[i][j] < eps) nxz++;
                    if (!std::isnan(d_yz[i][j]) && d_yz[i][j] < eps) nyz++;
                    if (!std::isnan(d_z[i][j])  && d_z[i][j]  < eps) nz++;
                }
                else
                {
                    if (!std::isnan(d_xz[i][j]) && d_xz[i][j] <= eps) nxz++;
                    if (!std::isnan(d_yz[i][j]) && d_yz[i][j] <= eps) nyz++;
                    if (!std::isnan(d_z[i][j])  && d_z[i][j]  <= eps) nz++;    
                } 
            }

            if (alg == 0)
                sum += NumericUtils::Digamma(nxz+1)
                     + NumericUtils::Digamma(nyz+1)
                     - NumericUtils::Digamma(nz+1);
            else
                sum += NumericUtils::Digamma(nxz)
                     + NumericUtils::Digamma(nyz)
                     - NumericUtils::Digamma(nz);
        }

        avg_log_eps /= n;

        double cmi = NumericUtils::Digamma(k) - sum / n;
        if (alg == 1) cmi -= 1.0 / k;

        if (!normalize)
        {
            if (!NumericUtils::doubleNearlyEqual(base,std::exp(1.0)))
                cmi /= std::log(base);

            return cmi;
        } 

        double hxy_z = NumericUtils::Digamma(n)
                     - NumericUtils::Digamma(k)
                     + d * avg_log_eps;
        if (alg == 1) hxy_z += 1.0 / k;

        if (hxy_z <= 0)
        {
            if (!NumericUtils::doubleNearlyEqual(base,std::exp(1.0)))
                cmi /= std::log(base);

            return cmi;
        } 

        return cmi / hxy_z;
    }

} // namespace KSGInfo

#endif // KSGINFO_HPP
