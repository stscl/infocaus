#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <numeric>
#include <cstdint>
#include <algorithm>
#include <unordered_map> 
#include "infoxtr.h"

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

    return infoxtr::infotheo::entropy(s, base, na_rm);
}

// Wrapper function to calculate Shannon entropy for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContEntropy(const Rcpp::NumericVector& vec,
                       int k = 3, 
                       int alg = 0,
                       double base = 2.0)
{
    std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
    
    return infoxtr::ksginfo::entropy(
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
    infoxtr::infotheo::Matrix m = infoxtr::convert::pat_r2std(mat,false);
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
    
    return infoxtr::infotheo::je(m, v, base, na_rm);
}

// Wrapper function to calculate joint entropy for continuous data
// [[Rcpp::export(rng = false)]]
double RcppContJE(const Rcpp::NumericMatrix& mat,
                  const Rcpp::IntegerVector& vars,
                  int k = 3, 
                  int alg = 0,
                  double base = 2.0)
{
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);
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
    
    return infoxtr::ksginfo::je(
                m, v, static_cast<size_t>(std::abs(k)), 
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
    infoxtr::infotheo::Matrix m = infoxtr::convert::pat_r2std(mat,false);

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

    return infoxtr::infotheo::ce(m, t, c, base, na_rm);
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
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

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
    
    return infoxtr::ksginfo::ce(
                m, t, c, static_cast<size_t>(std::abs(k)), 
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
    infoxtr::infotheo::Matrix m = infoxtr::convert::pat_r2std(mat,false);

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

    return infoxtr::infotheo::mi(m, t, i, base, na_rm, normalize);
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
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

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
    
    return infoxtr::ksginfo::mi(
                m, t, i, static_cast<size_t>(std::abs(k)), 
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
    infoxtr::infotheo::Matrix m = infoxtr::convert::pat_r2std(mat,false);

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

    return infoxtr::infotheo::cmi(m, t, i, c, base, na_rm, normalize);
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
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

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
    
    return infoxtr::ksginfo::cmi(
                m, t, i, c, static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(alg)), base, normalize);
}

// Wrapper function to calculate transfer entropy for discrete time series data
// [[Rcpp::export(rng = false)]]
double RcppDiscTE(SEXP mat,
                  const Rcpp::IntegerVector& target,
                  const Rcpp::IntegerVector& agent,
                  int lag_p = 3,
                  int lag_q = 3,
                  double base = 2.0,
                  bool na_rm = true,
                  bool normalize = false,
                  bool lag_single = false)
{
    infoxtr::infotheo::Matrix m = infoxtr::convert::pat_r2std(mat, false);

    std::vector<size_t> tg = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> ag = Rcpp::as<std::vector<size_t>>(agent);

    const size_t n_cols = m.size();
    for (auto& idx : tg) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : ag) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }

    return infoxtr::transferentropy::transferentropy(
                m, tg, ag, static_cast<size_t>(std::abs(lag_p)), 
                static_cast<size_t>(std::abs(lag_q)), base, 
                na_rm, normalize, lag_single);
}

// Wrapper function to calculate transfer entropy for continuous time series data
// [[Rcpp::export(rng = false)]]
double RcppContTE(const Rcpp::NumericMatrix& mat,
                  const Rcpp::IntegerVector& target,
                  const Rcpp::IntegerVector& agent,
                  int lag_p = 3,
                  int lag_q = 3,
                  int k = 3, 
                  int alg = 0,
                  double base = 2.0,
                  bool normalize = false,
                  bool lag_single = false)
{
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

    std::vector<size_t> tg = Rcpp::as<std::vector<size_t>>(target);
    std::vector<size_t> ag = Rcpp::as<std::vector<size_t>>(agent);

    const size_t n_cols = m.size();
    for (auto& idx : tg) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    for (auto& idx : ag) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Interact index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    
    return infoxtr::transferentropy::transferentropy(
                m, tg, ag, 
                static_cast<size_t>(std::abs(lag_p)), 
                static_cast<size_t>(std::abs(lag_q)), 
                static_cast<size_t>(std::abs(k)), 
                static_cast<size_t>(std::abs(alg)), 
                base, normalize, lag_single);
}

// Wrapper function to preform SURD decomposition for discrete data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppDiscSURD(SEXP mat,
                        int max_order = 3,
                        int threads = 1,
                        double base = 2.0,
                        bool normalize = true)
{
    infoxtr::surd::DiscMat m = infoxtr::convert::pat_r2std(mat, false);

    infoxtr::surd::SURDRes res = infoxtr::surd::surd(
        m, static_cast<size_t>(std::abs(max_order)),
        static_cast<size_t>(std::abs(threads)), base, normalize);

    const size_t k = res.size();

    Rcpp::NumericVector values(k);
    Rcpp::CharacterVector types(k);
    Rcpp::CharacterVector names(k);

    for (size_t i = 0; i < k; ++i)
    {
        values[i] = res.values[i];

        // variable name
        if (res.types[i] == 3)
        {
            // InfoLeak uses all sources
            std::string nm = "InfoLeak";
            names[i] = nm;
            types[i] = "InfoLeak";
            continue;
        }

        const auto& vars = res.var_indices[i];

        std::string nm;

        for (size_t j = 0; j < vars.size(); ++j)
        {
            if (j > 0)
                nm += "_";

            nm += "V";
            nm += std::to_string(vars[j]);
        }

        names[i] = nm;

        switch (res.types[i])
        {
            case 0:
                types[i] = "R";
                break;
            case 1:
                types[i] = "U";
                break;
            case 2:
                types[i] = "S";
                break;
            default:
                types[i] = "Unknown";
        }
    }

    // values.attr("names") = names;

    return Rcpp::List::create(
        Rcpp::Named("vars")  = names,
        Rcpp::Named("types")  = types,
        Rcpp::Named("values") = values
    );
}

// Wrapper function to preform SURD decomposition for continuous data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppContSURD(const Rcpp::NumericMatrix& mat,
                        int max_order = 3,
                        int k = 3,
                        int alg = 0,
                        int threads = 1,
                        double base = 2.0,
                        bool normalize = true)
{
    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

    infoxtr::surd::SURDRes res = infoxtr::surd::surd(
        m, static_cast<size_t>(std::abs(max_order)),
        static_cast<size_t>(std::abs(k)),
        static_cast<size_t>(std::abs(alg)),
        static_cast<size_t>(std::abs(threads)), 
        base, normalize);

    const size_t n_vals = res.size();

    Rcpp::NumericVector values(n_vals);
    Rcpp::CharacterVector types(n_vals);
    Rcpp::CharacterVector names(n_vals);

    for (size_t i = 0; i < n_vals; ++i)
    {
        values[i] = res.values[i];

        // variable name
        if (res.types[i] == 3)
        {
            // InfoLeak uses all sources
            std::string nm = "InfoLeak";
            names[i] = nm;
            types[i] = "InfoLeak";
            continue;
        }

        const auto& vars = res.var_indices[i];

        std::string nm;

        for (size_t j = 0; j < vars.size(); ++j)
        {
            if (j > 0)
                nm += "_";

            nm += "V";
            nm += std::to_string(vars[j]);
        }

        names[i] = nm;

        switch (res.types[i])
        {
            case 0:
                types[i] = "R";
                break;
            case 1:
                types[i] = "U";
                break;
            case 2:
                types[i] = "S";
                break;
            default:
                types[i] = "Unknown";
        }
    }

    // values.attr("names") = names;

    return Rcpp::List::create(
        Rcpp::Named("vars")  = names,
        Rcpp::Named("types")  = types,
        Rcpp::Named("values") = values
    );
}
