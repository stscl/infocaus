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

enum class KSGAlgorithm
{
    Alg1 = 1,
    Alg2 = 2
};

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
    size_t k,
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
        std::vector<double> row;

        for (size_t j = 0; j < n; ++j)
        {
            if (i == j) continue;
            row.push_back(dist[i][j]);
        }

        std::nth_element(row.begin(),row.begin()+k-1,row.end());

        double eps = row[k-1];

        avg += std::log(eps);
    }

    avg /= n;

    double H =
        NumericUtils::Digamma(n)
        - NumericUtils::Digamma(k)
        + avg
        + std::log(2.0);

    return H;
}

/***********************************************************
Joint Entropy
***********************************************************/
inline double JE(
    const Matrix& mat,
    const std::vector<size_t>& vars,
    size_t k = 3)
{
    Matrix sub = subset(mat,vars);

    const size_t d = sub.size();
    const size_t n = sub[0].size();

    auto dist = Dist::Dist(sub,"maximum",true,false);

    double avg = 0.0;

    for (size_t i=0;i<n;++i)
    {
        std::vector<double> row;

        for (size_t j=0;j<n;++j)
        {
            if (i==j) continue;
            row.push_back(dist[i][j]);
        }

        std::nth_element(row.begin(),row.begin()+k-1,row.end());

        double eps = row[k-1];

        avg += std::log(eps);
    }

    avg /= n;

    double H =
        NumericUtils::Digamma(n)
        - NumericUtils::Digamma(k)
        + d * avg
        + d * std::log(2.0);

    return H;
}

/***********************************************************
Conditional Entropy
***********************************************************/
inline double CE(
    const Matrix& mat,
    const std::vector<size_t>& target,
    const std::vector<size_t>& cond,
    size_t k = 3)
{
    std::vector<size_t> tc = cond;
    tc.insert(tc.end(),target.begin(),target.end());

    return JE(mat,tc,k) - JE(mat,cond,k);
}

/***********************************************************
Mutual Information (KSG estimator)
***********************************************************/
inline double MI(
    const Matrix& mat,
    const std::vector<size_t>& target,
    const std::vector<size_t>& interact,
    const std::vector<size_t>& xvars,
    const std::vector<size_t>& yvars,
    size_t k = 3,
    KSGAlgorithm alg = KSGAlgorithm::Alg1,
    bool normalize = false)
{
    std::vector<size_t> ti = target;
    ti.insert(ti.end(), interact.begin(), interact.end());

    auto d_xy = Dist::Dist(subset(mat,ti),"maximum",true,false);
    auto d_x  = Dist::Dist(subset(mat,target),"maximum",true,false);
    auto d_y  = Dist::Dist(subset(mat,interact),"maximum",true,false);

    const size_t n = d_xy.size();

    double sum = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        std::vector<double> row;

        for (size_t j = 0; j<n; ++j)
        {
            if (i==j) continue;
            row.push_back(d_xy[i][j]);
        }

        std::nth_element(row.begin(),row.begin()+k-1,row.end());

        double eps = row[k-1];

        size_t nx=0, ny=0;

        for (size_t j=0;j<n;++j)
        {
            if (i==j) continue;

            if (d_x[i][j] <= eps) nx++;
            if (d_y[i][j] <= eps) ny++;
        }

        if (alg == KSGAlgorithm::Alg1)
        {
            sum +=
                NumericUtils::Digamma(nx+1)
              + NumericUtils::Digamma(ny+1);
        }
        else
        {
            sum +=
                NumericUtils::Digamma(nx)
              + NumericUtils::Digamma(ny);
        }
    }

    double mi;

    if (alg == KSGAlgorithm::Alg1)
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
        return mi;

    double hxy = JE(mat,xy,k);

    if (hxy <= 0)
        return mi;

    return mi / hxy;
}

/***********************************************************
Conditional Mutual Information
***********************************************************/
inline double CMI(
    const Matrix& mat,
    const std::vector<size_t>& xvars,
    const std::vector<size_t>& yvars,
    const std::vector<size_t>& zvars,
    size_t k = 3,
    KSGAlgorithm alg = KSGAlgorithm::Alg1,
    bool normalize = false)
{
    std::vector<size_t> xyz = zvars;
    xyz.insert(xyz.end(),xvars.begin(),xvars.end());
    xyz.insert(xyz.end(),yvars.begin(),yvars.end());

    std::vector<size_t> xz = zvars;
    xz.insert(xz.end(),xvars.begin(),xvars.end());

    std::vector<size_t> yz = zvars;
    yz.insert(yz.end(),yvars.begin(),yvars.end());

    auto d_xyz = Dist::Dist(subset(mat,xyz),"maximum",true,false);
    auto d_xz  = Dist::Dist(subset(mat,xz),"maximum",true,false);
    auto d_yz  = Dist::Dist(subset(mat,yz),"maximum",true,false);
    auto d_z   = Dist::Dist(subset(mat,zvars),"maximum",true,false);

    const size_t n = d_xyz.size();

    double sum = 0.0;

    for (size_t i=0;i<n;++i)
    {
        std::vector<double> row;

        for (size_t j=0;j<n;++j)
        {
            if (i==j) continue;
            row.push_back(d_xyz[i][j]);
        }

        std::nth_element(row.begin(),row.begin()+k-1,row.end());

        double eps = row[k-1];

        size_t nxz=0, nyz=0, nz=0;

        for (size_t j=0;j<n;++j)
        {
            if (i==j) continue;

            if (d_xz[i][j] <= eps) nxz++;
            if (d_yz[i][j] <= eps) nyz++;
            if (d_z[i][j]  <= eps) nz++;
        }

        sum +=
            NumericUtils::Digamma(nxz+1)
          + NumericUtils::Digamma(nyz+1)
          - NumericUtils::Digamma(nz+1);
    }

    double cmi =
        NumericUtils::Digamma(k)
      - sum / n;

    if (!normalize)
        return cmi;

    double ce_xy_z = CE(mat,xvars,yvars,k);

    if (ce_xy_z <= 0)
        return cmi;

    return cmi / ce_xy_z;
}

} // namespace KSGInfo

#endif // KSGINFO_HPP
