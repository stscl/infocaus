#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include "infoxtr.h"

// Wrapper function to calculate spatial lag value for vector spatial data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenLatticeLag(const Rcpp::NumericMatrix& mat,
                                      const Rcpp::List& nb, 
                                      int lag = 1) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> cppMat = infoxtr::convert::mat_r2std(mat, true);

  // Convert Rcpp::List to std::vector<std::vector<size_t>>
  std::vector<std::vector<size_t>> nb_std = infoxtr::convert::nb2std(nb);

  // Calculate lagged values
  std::vector<std::vector<double>> lagged_values =
    infoxtr::lagg::lagg(cppMat, nb_std, static_cast<size_t>(std::abs(lag)));

  return infoxtr::convert::mat_std2r(lagged_values, true);
}

// Wrapper function to calculate spatial lag value for raster spatial data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenGridLag(const Rcpp::NumericMatrix& mat,
                                   int nrows, int lag = 1) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> cppMat = infoxtr::convert::mat_r2std(mat, true);

  // Calculate lagged values
  std::vector<std::vector<double>> lagged_values =
    infoxtr::lagg::lagg(cppMat, 
                        static_cast<size_t>(std::abs(nrows)), 
                        static_cast<size_t>(std::abs(lag)));

  return infoxtr::convert::mat_std2r(lagged_values, true);
}

// Wrapper function to calculate temporal lag value for time series data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenTSLag(const Rcpp::NumericMatrix& mat,
                                 int lag = 1) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> cppMat = infoxtr::convert::mat_r2std(mat, true);

  // Calculate lagged values
  std::vector<std::vector<double>> lagged_values =
    infoxtr::lagg::lagg(cppMat, static_cast<size_t>(std::abs(lag)));

  return infoxtr::convert::mat_std2r(lagged_values, true);
}
