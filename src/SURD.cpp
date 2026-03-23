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
    if (tg_idx < 1 || tg_idx > n_cols) 
    {
        Rcpp::stop("Target index %d out of bounds [1, %d]", 
                   static_cast<int>(tg_idx), 
                   static_cast<int>(n_cols));
    }
    tg_idx -= 1; // to 0-based

    std::vector<size_t> ag_raw = Rcpp::as<std::vector<size_t>>(agent);
    for (auto& idx : ag_raw) 
    {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Agent index %d out of bounds [1, %d]",
                    static_cast<int>(idx),
                    static_cast<int>(n_cols));
        }
        idx -= 1;
    }
    size_t nag_raw = ag_raw.size();

    std::vector<size_t> bin_vec = Rcpp::as<std::vector<size_t>>(bin);
    std::vector<std::string> method_vec = Rcpp::as<std::vector<std::string>>(method);

    // Expand bin and method (target + agents)
    std::vector<size_t> bin_expanded(nag_raw + 1);
    std::vector<std::string> method_expanded(nag_raw + 1);

    // ---- bin ----
    if (bin_vec.size() == 1) 
    {
        std::fill(bin_expanded.begin(), bin_expanded.end(), bin_vec[0]);
    } 
    else if (bin_vec.size() == 2) 
    {
        bin_expanded[0] = bin_vec[0];
        std::fill(bin_expanded.begin() + 1, bin_expanded.end(), bin_vec[1]);
    } 
    else 
    {
        bin_expanded[0] = bin_vec[0];
        for (size_t i = 1; i < nag_raw + 1; ++i) 
        {
            size_t idx = (i - 1) % (bin_vec.size() - 1);
            bin_expanded[i] = bin_vec[idx + 1];
        }
    }

    // ---- method ----
    if (method_vec.size() == 1) 
    {
        std::fill(method_expanded.begin(), method_expanded.end(), method_vec[0]);
    } 
    else if (method_vec.size() == 2) 
    {
        method_expanded[0] = method_vec[0];
        std::fill(method_expanded.begin() + 1, method_expanded.end(), method_vec[1]);
    } 
    else 
    {
        method_expanded[0] = method_vec[0];
        for (size_t i = 1; i < nag_raw + 1; ++i) 
        {
            size_t idx = (i - 1) % (method_vec.size() - 1);
            method_expanded[i] = method_vec[idx + 1];
        }
    }

    // Create ordering index
    std::vector<size_t> order(ag_raw.size(), 0);
    std::iota(order.begin(), order.end(), 0);

    // Sort indices by ag_raw
    std::sort(order.begin(), order.end(),
            [&](size_t i, size_t j) {
                return ag_raw[i] < ag_raw[j];
            });

    // Apply ordering
    std::vector<size_t> ag_sorted;
    std::vector<size_t> bin_sorted;
    std::vector<std::string> method_sorted;

    for (size_t idx : order) 
    {
        ag_sorted.push_back(ag_raw[idx]);
        bin_sorted.push_back(bin_expanded[idx + 1]);
        method_sorted.push_back(method_expanded[idx + 1]);
    }

    // Deduplicate while keeping alignment
    std::vector<size_t> ag;
    std::vector<size_t> bin_final;
    std::vector<std::string> method_final;

    for (size_t i = 0; i < ag_sorted.size(); ++i) 
    {
        if (i == 0 || ag_sorted[i] != ag_sorted[i - 1]) 
        {
            ag.push_back(ag_sorted[i]);
            bin_final.push_back(bin_sorted[i]);
            method_final.push_back(method_sorted[i]);
        }
    }

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
        vec, method_expanded[0], bin_expanded[0]
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
            lagged_values[j], method_final[j], bin_final[j]
        );
    }

    infoxtr::surd::SURDRes res = infoxtr::surd::surd(
        pm, 
        static_cast<size_t>(std::abs(max_order)),
        static_cast<size_t>(std::abs(threads)), 
        std::abs(base), normalize);

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
