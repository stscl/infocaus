#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include "lag.hpp"
#include "DataTrans.h"

// Wrapper function to calculate spatial lag value for vector spatial data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenLatticeLag(const Rcpp::NumericMatrix& mat,
                                      const Rcpp::List& nb, 
                                      int lag = 1) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> cppMat = mat_r2std(mat, byrow = true);

  // Convert Rcpp::List to std::vector<std::vector<size_t>>
  std::vector<std::vector<size_t>> nb_std = nb2std(nb);

  // Calculate lagged values
  std::vector<std::vector<double>> lagged_values =
    Lag::GenLatticeLag(cppMat, nb_std, static_cast<size_t>(std::abs(lag)));

  return mat_r2std(lagged_values, byrow = true);
}
