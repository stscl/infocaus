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
                    Rcpp::Nullable<int> nrow = R_NilValue)
{   
    const size_t n_cols = mat.size();

    if (tg[0] < 1 || tg[0] > n_cols) {
            Rcpp::stop("Target index %d out of bounds [1, %d]", 
                       static_cast<int>(tg[0]), 
                       static_cast<int>(n_cols));
        }
    const size_t tg_idx = tg[0] - 1; // to 0-based
    
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
        Rcpp::stop("Agent index should not be empty");

    std::vector<std::vector<double>> m = infoxtr::convert::mat_r2std(mat, false);

    // Convert Rcpp::List to std::vector<std::vector<size_t>>
    std::vector<std::vector<size_t>> nb_std;
    if (nb.isNotNull()) {
        nb_std = infoxtr::convert::nb2std(nb);
    }

    // Process each selected variable
    for (size_t j = 0; j < n_vars; ++j)
    {
        size_t col_id = v[j];

        // Extract column vector from R matrix
        std::vector<double> vec(n_obs);
        for (size_t r = 0; r < n_obs; ++r)
        {
            vec[r] = mat(r, col_id);
        }

        // Generate lattice embedding
        std::vector<std::vector<double>> embeddings =
            Embed::GenLatticeEmbedding(
                vec,
                nb_std,
                static_cast<size_t>(std::abs(E[j])),
                static_cast<size_t>(std::abs(tau[j])),
                static_cast<size_t>(std::abs(style[j]))
            );

        // Convert continuous embedding -> symbolic patterns
        pm[j] = SymDync::GenSymbolicPattern(
            embeddings,
            relative,
            na_rm
        );
    }

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
