#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <numeric>
#include <algorithm>
#include "discretize.hpp"
#include "DataTrans.h"

// Wrapper function to discretize continuous numeric vector
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppDisc(
    const Rcpp::NumericVector& vec,
    int n = 5,
    const std::string& method = "natural",
    int sample_begin = 3000,
    double sample_prob = 0.15,
    int seed = 123456789,
    double threshold = 0.4,
    int iter_step = 100,
    Rcpp::Nullable<Rcpp::NumericVector> breakpoints = R_NilValue,
    bool right_closed = true
) {
    // Convert Rcpp::NumericVector to std::vector<double>
    std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
    
    std::vector<double> bp_std;
    if (breakpoints.isNotNull()) {
        bp_std = Rcpp::as<std::vector<double>>(breakpoints);
    }
    
    // Check the discretization method
    static const std::vector<std::string> methods = {
        "sd",
        "equal",
        "geometric",
        "quantile",
        "manual",
        "natural",
        "jenks",
        "headtail",
        "headtails"
    };

    if (std::find(methods.begin(), methods.end(), method) == methods.end())
    {
        Rcpp::stop(
            "Unknown discretization method. Available methods are: "
            "sd, equal, geometric, quantile, manual, natural (or jenks) and headtail(s)"
        );
    }

    // Call the discretization function
    std::vector<uint64_t> discv = 
        discretize::discretize(
            vec_std, method, 
            static_cast<size_t>(std::abs(n)), 
            static_cast<size_t>(std::abs(sample_begin)),
            sample_prob,
            static_cast<uint64_t>(std::abs(seed)),
            threshold,
            static_cast<size_t>(std::abs(iter_step)),
            bp_std, right_closed);

    // Convert the result back to Rcpp::IntegerVector
    return Rcpp::wrap(discv);
}
