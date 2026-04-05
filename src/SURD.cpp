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
                    int n = 5,
                    int max_order = 3,
                    int threads = 1,
                    double base = 2.0,
                    bool normalize = true,
                    const std::string& method = "equal",
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
        size_t col_id = ag[j];

        for (size_t r = 0; r < n_obs; ++r)
        {
            cppMat[r][col_id] = mat(r, col_id);
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
