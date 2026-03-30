#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <numeric>
#include <algorithm>
#include "distance.hpp"
#include "DataTrans.h"

// Wrapper function to compute the distance between two vectors
// [[Rcpp::export(rng = false)]]
double RcppDist4Vec(
    const Rcpp::NumericVector& v1,
    const Rcpp::NumericVector& v2,
    std::string method = "euclidean",
    bool na_rm = true
) {
    // Convert Rcpp::NumericVector to std::vector<double>
    std::vector<double> v1_std = Rcpp::as<std::vector<double>>(v1);
    std::vector<double> v2_std = Rcpp::as<std::vector<double>>(v2);

    // Call the distance function
    double distv = distance::distance(v1_std, v2_std, method, na_rm);

    return distv;
}

// Wrapper function to compute a row-wise distance matrix for an input matrix
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppDist4Mat(
    const Rcpp::NumericMatrix& mat,
    std::string method = "euclidean",
    bool na_rm = true,
    bool byrow = true
) {
    // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
    std::vector<std::vector<double>> cppMat = mat_r2std(mat, true);

    // Call the distance function
    std::vector<std::vector<double>> distm = distance::distance(cppMat, method, na_rm, byrow);

    // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix and return
    return mat_std2r(distm, true);
}

// Wrapper function to compute a row-wise distance matrix for an input matrix subset
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppDist4MatSub(
    const Rcpp::NumericMatrix& mat,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    std::string method = "euclidean",
    bool na_rm = true,
    bool byrow = true
) {
    // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
    std::vector<std::vector<double>> cppMat = mat_r2std(mat, true);

    const int n_obs = byrow ? mat.nrow() : mat.ncol(); 

    // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
    std::vector<size_t> lib_std;
    lib_std.reserve(lib.size());
    for (int i = 0; i < lib.size(); ++i) {
        if (lib[i] < 1 || lib[i] > n_obs) {
            Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
        }
        lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }

    std::vector<size_t> pred_std;
    pred_std.reserve(pred.size());
    for (int i = 0; i < pred.size(); ++i) {
        if (pred[i] < 1 || pred[i] > n_obs) {
            Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
        }
        pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }

    // Call the distance function
    std::vector<std::vector<double>> distm = distance::distance(
        cppMat, lib_std, pred_std, method, na_rm, byrow);

    // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix and return
    return mat_std2r(distm, true);
}
