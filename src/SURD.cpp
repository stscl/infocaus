#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include "infoxtr.h"

// Wrapper function to preform synergistic–unique–redundant decomposition
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSURD(const Rcpp::NumericMatrix& mat,
                    const Rcpp::IntegerVector& target,
                    const Rcpp::IntegerVector& agent,
                    int lag = 1,
                    int max_order = 3,
                    int threads = 1,
                    double base = 2.0,
                    bool normalize = true,
                    const Rcpp::IntegerVector& bin = Rcpp::IntegerVector::create(5),
                    const Rcpp::CharacterVector& method = Rcpp::CharacterVector::create("equal"),
                    Rcpp::Nullable<Rcpp::List> nb = R_NilValue,
                    Rcpp::Nullable<int> nrows = R_NilValue)
{   
    const size_t n_cols = static_cast<size_t>(mat.ncol());
    const size_t n_obs = static_cast<size_t>(mat.nrow());
    
    size_t tg_idx = target[0];
    if (tg_idx < 1 || tg_idx > n_cols) {
        Rcpp::stop("Target index %d out of bounds [1, %d]", 
                   static_cast<int>(tg_idx), 
                   static_cast<int>(n_cols));
    }
    tg_idx -= 1; // to 0-based
    
    std::vector<size_t> ag = Rcpp::as<std::vector<size_t>>(agent);
    for (auto& idx : ag) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Agent index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    std::sort(ag.begin(), ag.end());
    ag.erase(
        std::unique(ag.begin(), ag.end()),
        ag.end()
    );
    if (ag.empty())
        Rcpp::stop("Agent indices should not be empty");

    // Construct discrete data matrix
    std::vector<std::vector<uint64_t>> pm(
        ag.size() + 1, std::vector<uint64_t>(n_obs, 0)
    );
    
    // Preserve original values in target variable and discretize it
    std::vector<double> vec(n_obs);
    for (size_t r = 0; r < n_obs; ++r)
    {
        vec[r] = mat(r, tg_idx);
    }
    pm[0] = infoxtr::discretize::discretize(
        vec, method, static_cast<size_t>(std::abs(n))
    );
    
    // Generate lagged values for agent variables
    std::vector<std::vector<double>> lagged_values;

    std::vector<std::vector<double>> cppMat(
        n_obs, std::vector<double>(ag.size())
    );
    for (size_t j = 0; j < ag.size(); ++j)
    {   
        for (size_t r = 0; r < n_obs; ++r)
        {
            cppMat[r][j] = mat(r, ag[j]);
        }
    }

    if (nb.isNotNull()) 
    {
        // Convert Rcpp::List to std::vector<std::vector<size_t>>
        std::vector<std::vector<size_t>> nb_std = infoxtr::convert::nb2std(nb.get());
        lagged_values = infoxtr::lagg::lagg(
            cppMat, nb_std, static_cast<size_t>(std::abs(lag)), false);
    } 
    else if (nrows.isNotNull())
    {
        lagged_values = infoxtr::lagg::lagg(
            cppMat, 
            static_cast<size_t>(std::abs(Rcpp::as<int>(nrows))), 
            static_cast<size_t>(std::abs(lag)), false);
    }
    else  
    {
        lagged_values = infoxtr::lagg::lagg(
            cppMat, static_cast<size_t>(std::abs(lag)), false);
    }

    // Discrete lagged values for agent variables
    for (size_t j = 0; j < lagged_values.size(); ++j)
    {   
        pm[j + 1] = infoxtr::discretize::discretize(
            lagged_values[j], method, static_cast<size_t>(std::abs(n))
        );
    }

    infoxtr::surd::SURDRes res = infoxtr::surd::surd(
        pm, static_cast<size_t>(std::abs(max_order)),
        static_cast<size_t>(std::abs(threads)), base, normalize);

    std::vector<std::string> names;
    std::vector<std::string> types;
    std::vector<double> values;

    auto make_name = [](const std::vector<size_t>& vars)
    {
        std::string nm;

        for (size_t j = 0; j < vars.size(); ++j)
        {
            if (j > 0)
                nm += "_";

            nm += "V";
            nm += std::to_string(vars[j]);
        }

        return nm;
    };

    /**************************************************
     * Unique
     **************************************************/

    for (size_t i = 0; i < res.unique_vals.size(); ++i)
    {
        names.push_back(make_name(res.unique_vars[i]));
        types.push_back("U");
        values.push_back(res.unique_vals[i]);
    }

    /**************************************************
     * Redundant
     **************************************************/

    for (size_t i = 0; i < res.redundant_vals.size(); ++i)
    {
        names.push_back(make_name(res.redundant_vars[i]));
        types.push_back("R");
        values.push_back(res.redundant_vals[i]);
    }

    /**************************************************
     * Synergy
     **************************************************/

    for (size_t i = 0; i < res.synergy_vals.size(); ++i)
    {
        names.push_back(make_name(res.synergy_vars[i]));
        types.push_back("S");
        values.push_back(res.synergy_vals[i]);
    }

    /**************************************************
     * InfoLeak
     **************************************************/

    names.push_back("InfoLeak");
    types.push_back("InfoLeak");
    values.push_back(res.info_leak);

    Rcpp::CharacterVector names_r(names.begin(), names.end());
    Rcpp::CharacterVector types_r(types.begin(), types.end());
    Rcpp::NumericVector values_r(values.begin(), values.end());

    return Rcpp::List::create(
        Rcpp::Named("vars")  = names_r,
        Rcpp::Named("types")  = types_r,
        Rcpp::Named("values") = values_r
    );
}
