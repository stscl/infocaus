#include <vector>
#include <cmath>
#include <cstdint>
#include <string>
#include <limits>
#include <numeric>
#include <algorithm>
#include <unordered_map> 
#include "infotheo.hpp"
#include "ksginfo.hpp"
#include "DataTrans.h"

// Wrapper function to calculate Shannon entropy for discrete data
// [[Rcpp::export(rng = false)]]
double RcppDiscEntropy(SEXP series,
                       double base = 2.0,
                       bool na_rm = true)
{
    std::vector<uint64_t> s;

    switch (TYPEOF(series))
    {
        case INTSXP:
        {   
            Rcpp::IntegerVector v(series);
            s.reserve(v.size());

            std::vector<int> uniq;
            uniq.reserve(v.size());

            for (int i = 0; i < v.size(); ++i)
                if (!Rcpp::IntegerVector::is_na(v[i]))
                    uniq.push_back(v[i]);

            std::sort(uniq.begin(), uniq.end());
            uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

            std::unordered_map<int, uint64_t> dict;
            for (uint64_t i = 0; i < uniq.size(); ++i)
                dict[uniq[i]] = i+1;

            for (int i = 0; i < v.size(); ++i)
            {
                if (Rcpp::IntegerVector::is_na(v[i]))
                {
                    s.push_back( 0 );
                }
                else
                {
                    s.push_back( dict[v[i]] );
                }
            }
            
            break;
        }

        case REALSXP:
        {
            Rcpp::NumericVector v(series);
            s.reserve(v.size());

            std::vector<double> uniq;
            uniq.reserve(v.size());

            for (int i = 0; i < v.size(); ++i)
                if (!Rcpp::NumericVector::is_na(v[i]))
                    uniq.push_back(v[i]);

            std::sort(uniq.begin(), uniq.end());
            uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

            std::unordered_map<double, uint64_t> dict;
            for (uint64_t i = 0; i < uniq.size(); ++i)
                dict[uniq[i]] = i+1;

            for (int i = 0; i < v.size(); ++i)
            {
                if (Rcpp::NumericVector::is_na(v[i]))
                {
                    s.push_back( 0 );
                }
                else
                {
                    s.push_back( dict[v[i]] );
                }
            }

            break;
        }

        case STRSXP:
        {
            Rcpp::CharacterVector v(series);
            s.reserve(v.size());

            std::vector<std::string> uniq;
            uniq.reserve(v.size());

            for (int i = 0; i < v.size(); ++i)
                if (!Rcpp::CharacterVector::is_na(v[i]))
                    uniq.push_back(std::string(v[i]));

            std::sort(uniq.begin(), uniq.end());
            uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

            std::unordered_map<std::string, uint64_t> dict;
            for (uint64_t i = 0; i < uniq.size(); ++i)
                dict[uniq[i]] = i+1;

            for (int i = 0; i < v.size(); ++i)
            {
                if (Rcpp::CharacterVector::is_na(v[i]))
                {
                    s.push_back( 0 );
                }
                else
                {
                    s.push_back( dict[std::string(v[i])] );
                }
            }

            break;
        }

        default:
            Rcpp::stop("Input must be Integer, Numeric, or Character vector.");
    }

    return InfoTheo::Entropy(s, base, na_rm);
}

// Wrapper function to calculate Shannon entropy for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContEntropy(const Rcpp::NumericVector& vec,
                       int k = 3, 
                       int alg = 0,
                       double base = 2.0)
{
    std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
    
    return KSGInfo::Entropy(
        vec_std, static_cast<size_t>(std::abs(k)), 
        static_cast<size_t>(std::abs(alg)), base);
}

// Wrapper function to calculate joint entropy for discrete data
// [[Rcpp::export(rng = false)]]
double RcppDiscJE(SEXP mat,
                  const Rcpp::IntegerVector& vars,
                  double base = 2.0,
                  bool na_rm = true)
{
    InfoTheo::Matrix m = pat_r2std(mat,false);
    std::vector<size_t> v = Rcpp::as<std::vector<size_t>>(vars);

    const size_t n_cols = m.size();
    for (auto& idx : v) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Column index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return InfoTheo::JE(m, v, base, na_rm);
}

// Wrapper function to calculate joint entropy for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContJE(const Rcpp::NumericMatrix& mat,
                  const Rcpp::IntegerVector& vars,
                  int k = 3, 
                  int alg = 0,
                  double base = 2.0)
{
    std::vector<std::vector<double>> m = mat_r2std(mat, false);
    std::vector<size_t> v = Rcpp::as<std::vector<size_t>>(vars);

    const size_t n_cols = m.size();
    for (auto& idx : v) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Column index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return KSGInfo::JE(m, v, static_cast<size_t>(std::abs(k)), 
                       static_cast<size_t>(std::abs(alg)), base);
}

// Wrapper function to calculate conditional entropy for discrete data
// [[Rcpp::export(rng = false)]]
double RcppDiscCE(SEXP mat,
                  const Rcpp::IntegerVector& target,
                  const Rcpp::IntegerVector& conds,
                  double base = 2.0,
                  bool na_rm = true)
{
    InfoTheo::Matrix m = pat_r2std(mat,false);

    std::vector<size_t> t = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> c = Rcpp::as<std::vector<size_t>>(conds);

    const size_t n_cols = m.size();
    for (auto& idx : t) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : c) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Conds index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }

    return InfoTheo::CE(m, t, c, base, na_rm);
}

// Wrapper function to calculate conditional entropy for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContCE(const Rcpp::NumericMatrix& mat,
                  const Rcpp::IntegerVector& target,
                  const Rcpp::IntegerVector& conds,
                  int k = 3, 
                  int alg = 0,
                  double base = 2.0)
{
    std::vector<std::vector<double>> m = mat_r2std(mat, false);

    std::vector<size_t> t = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> c = Rcpp::as<std::vector<size_t>>(conds);

    const size_t n_cols = m.size();
    for (auto& idx : t) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : c) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Conds index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return KSGInfo::CE(m, t, c, static_cast<size_t>(std::abs(k)), 
                       static_cast<size_t>(std::abs(alg)), base);
}

// Wrapper function to calculate mutual information for discrete data
// [[Rcpp::export(rng = false)]]
double RcppDiscMI(SEXP mat,
                  const Rcpp::IntegerVector& target,
                  const Rcpp::IntegerVector& interact,
                  double base = 2.0,
                  bool na_rm = true,
                  bool normalize = false)
{
    InfoTheo::Matrix m = pat_r2std(mat,false);

    std::vector<size_t> t = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> i = Rcpp::as<std::vector<size_t>>(interact);

    const size_t n_cols = m.size();
    for (auto& idx : t) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : i) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }

    return InfoTheo::MI(m, t, i, base, na_rm, normalize);
}

// Wrapper function to calculate mutual information for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContMI(const Rcpp::NumericMatrix& mat,
                  const Rcpp::IntegerVector& target,
                  const Rcpp::IntegerVector& interact,
                  int k = 3, 
                  int alg = 0,
                  double base = 2.0,
                  bool normalize = false)
{
    std::vector<std::vector<double>> m = mat_r2std(mat, false);

    std::vector<size_t> t = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> i = Rcpp::as<std::vector<size_t>>(interact);

    const size_t n_cols = m.size();
    for (auto& idx : t) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : i) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return KSGInfo::MI(m, t, i, static_cast<size_t>(std::abs(k)), 
                       static_cast<size_t>(std::abs(alg)), base, normalize);
}

// Wrapper function to calculate conditional mutual information for discrete data
// [[Rcpp::export(rng = false)]]
double RcppDiscCMI(SEXP mat,
                   const Rcpp::IntegerVector& target,
                   const Rcpp::IntegerVector& interact,
                   const Rcpp::IntegerVector& conds,
                   double base = 2.0,
                   bool na_rm = true,
                   bool normalize = false)
{
    InfoTheo::Matrix m = pat_r2std(mat,false);

    std::vector<size_t> t = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> i = Rcpp::as<std::vector<size_t>>(interact);
    std::vector<size_t> c = Rcpp::as<std::vector<size_t>>(conds);

    const size_t n_cols = m.size();
    for (auto& idx : t) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : i) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : c) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Conds index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }

    return InfoTheo::CMI(m, t, i, c, base, na_rm, normalize);
}

// Wrapper function to calculate conditional mutual information for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContCMI(const Rcpp::NumericMatrix& mat,
                   const Rcpp::IntegerVector& target,
                   const Rcpp::IntegerVector& interact,
                   const Rcpp::IntegerVector& conds,
                   int k = 3, 
                   int alg = 0,
                   double base = 2.0,
                   bool normalize = false)
{
    std::vector<std::vector<double>> m = mat_r2std(mat, false);

    std::vector<size_t> t = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> i = Rcpp::as<std::vector<size_t>>(interact);
    std::vector<size_t> c = Rcpp::as<std::vector<size_t>>(conds);

    const size_t n_cols = m.size();
    for (auto& idx : t) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : i) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : c) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Conds index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return KSGInfo::CMI(m, t, i, c, static_cast<size_t>(std::abs(k)), 
                       static_cast<size_t>(std::abs(alg)), base, normalize);
}
