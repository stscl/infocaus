#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <cstdint>
#include <algorithm>
#include "infoxtr.h"

// Wrapper function to quantify interventional causality by knockoff operation
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppKOCMI(
    const Rcpp::NumericMatrix& mat,
    const Rcpp::IntegerVector& target,
    const Rcpp::IntegerVector& agent,
    const Rcpp::IntegerVector& conds,
    const Rcpp::NumericMatrix& knockoff,
    Rcpp::Nullable<Rcpp::NumericMatrix> null_knockoff = R_NilValue,
    const std::string& type = "cont",
    int nboots = 10000,
    int k = 3,
    int alg = 0,
    int threads = 1,
    int seed = 123456789,
    double base = 2.0,
    const std::string& method = "equal",
    bool contain_null = true)
{   
    if (contain_null && !null_knockoff.isNotNull())
    {
        Rcpp::stop("[KOCMI] When `contain_null` is true, the `null_knockoff` matrix for the source variable must be supplied.");
    }

    if (nboots < 0) 
    {
        Rcpp::warning("`nboots` is negative, using absolute value: %d", nboots);
    }

    const size_t n_cols = static_cast<size_t>(mat.ncol());
    const size_t n_obs = static_cast<size_t>(mat.nrow());
    
    size_t tg_idx = target[0];
    if (tg_idx < 1 || tg_idx > n_cols) {
        Rcpp::stop("Target index %d out of bounds [1, %d]", 
                   static_cast<int>(tg_idx), 
                   static_cast<int>(n_cols));
    }
    tg_idx -= 1; // to 0-based

    size_t ag_idx = agent[0];
    if (ag_idx < 1 || ag_idx > n_cols) {
        Rcpp::stop("Agent index %d out of bounds [1, %d]", 
                   static_cast<int>(ag_idx), 
                   static_cast<int>(n_cols));
    }
    ag_idx -= 1; // to 0-based
    
    std::vector<size_t> cg = Rcpp::as<std::vector<size_t>>(conds);
    for (auto& idx : cg) {
        if (idx < 1 || idx > n_cols) {
            Rcpp::stop("Conditioning index %d out of bounds [1, %d]", 
                       static_cast<int>(idx), 
                       static_cast<int>(n_cols));
        }
        idx -= 1;  // to 0-based
    }
    std::sort(cg.begin(), cg.end());
    cg.erase(
        std::unique(cg.begin(), cg.end()),
        cg.end()
    );
    if (cg.empty())
        Rcpp::stop("Conditioning indices should not be empty");

    // Convert original values from R to C++
    std::vector<double> tg_std(n_obs);
    std::vector<double> ag_std(n_obs);
    for (size_t r = 0; r < n_obs; ++r)
    {
        tg_std[r] = mat(r, tg_idx);
        ag_std[r] = mat(r, ag_idx);
    }
    std::vector<std::vector<double>> cg_mat(
        n_obs, std::vector<double>(cg.size())
    );
    for (size_t j = 0; j < cg.size(); ++j)
    {   
        for (size_t r = 0; r < n_obs; ++r)
        {
            cg_mat[j][r] = mat(r, cg[j]);
        }
    }
    
    std::vector<std::vector<double>> nkm;
    if (null_knockoff.isNotNull())
    {
        nkm = infoxtr::convert::mat_r2std(null_knockoff.get(), false);
    }
    std::vector<std::vector<double>> km = infoxtr::convert::mat_r2std(knockoff, false);
    
    // Initialize result container
    infoxtr::kocmi::KOCMIRes res;
   
    if (type == "cont")
    {
        res = infoxtr::kocmi::kocmi(
            tg_std, ag_std, cg_mat, km, nkm,
            static_cast<size_t>(std::abs(nboots)),
            static_cast<size_t>(std::abs(k)),
            static_cast<size_t>(std::abs(alg)),
            static_cast<size_t>(std::abs(threads)),
            static_cast<uint64_t>(std::abs(seed)),
            contain_null
        );
    }
    else  
    {  
        std::vector<uint64_t> tg_vec = infoxtr::discretize::discretize(
            tg_std, method, static_cast<size_t>(std::abs(k))
        );

        std::vector<uint64_t> ag_vec = infoxtr::discretize::discretize(
            ag_std, method, static_cast<size_t>(std::abs(k))
        );

        std::vector<std::vector<uint64_t>> cg_discm(
            cg.size(), std::vector<uint64_t>(n_obs)
        );
        for (size_t j = 0; j < cg.size(); ++j)
        {   
            cg_discm[j] = infoxtr::discretize::discretize(
                cg_mat[j], method, static_cast<size_t>(std::abs(k))
            );    
        }

        std::vector<std::vector<uint64_t>> kdiscm(
            km.size(), std::vector<uint64_t>(n_obs)
        );
        for (size_t j = 0; j < km.size(); ++j)
        {   
            kdiscm[j] = infoxtr::discretize::discretize(
                km[j], method, static_cast<size_t>(std::abs(k))
            );    
        }

        std::vector<std::vector<uint64_t>> nkdiscm(
            nkm.size(), std::vector<uint64_t>(n_obs)
        );
        for (size_t j = 0; j < nkm.size(); ++j)
        {   
            nkdiscm[j] = infoxtr::discretize::discretize(
                nkm[j], method, static_cast<size_t>(std::abs(k))
            );    
        }

        res = infoxtr::kocmi::kocmi(
            tg_vec, ag_vec, cg_discm, kdiscm, nkdiscm,
            static_cast<size_t>(std::abs(nboots)),
            static_cast<size_t>(std::abs(threads)),
            static_cast<uint64_t>(std::abs(seed)),
            contain_null, base
        );
    }

    return Rcpp::NumericVector::create(
        Rcpp::Named("t_stat") = res.t_stat,
        Rcpp::Named("p_value") = res.p_value
    );
}
