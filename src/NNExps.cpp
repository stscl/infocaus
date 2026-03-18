#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <numeric>
#include <algorithm>
#include "neighbour.hpp"
#include "DataTrans.h"

// Wrapper function to compute the nearest neighbours for an input feature matrix
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppNN4Mat(
    const Rcpp::NumericMatrix& mat,
    int k,
    std::string method = "euclidean",
    bool include_self = false
) {
    // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
    int numRows = mat.nrow();
    int numCols = mat.ncol();
    std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            cppMat[r][c] = mat(r, c);
        }
    }

    // Call the neighbpurbood function
    std::vector<std::vector<size_t>> neighbours = NN::NN4Mat(
        cppMat, static_cast<size_t>(std::abs(k)), method, include_self);

    // Return nb object (List in R side)
    return std2nb(neighbours);
}

// Wrapper function to compute the nearest neighbours for an input feature matrix subset
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppNN4MatSub(
    const Rcpp::NumericMatrix& mat,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    int k,
    std::string method = "euclidean",
    bool include_self = false
) {
    // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
    int numRows = mat.nrow();
    int numCols = mat.ncol();
    std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            cppMat[r][c] = mat(r, c);
        }
    }

    // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
    std::vector<size_t> lib_std;
    lib_std.reserve(lib.size());
    for (int i = 0; i < lib.size(); ++i) {
        if (lib[i] < 1 || lib[i] > numRows) {
            Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
        }
        lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }

    std::vector<size_t> pred_std;
    pred_std.reserve(pred.size());
    for (int i = 0; i < pred.size(); ++i) {
        if (pred[i] < 1 || pred[i] > numRows) {
            Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
        }
        pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }

    // Call the neighbpurbood function
    std::vector<std::vector<size_t>> neighbours = NN::NN4Mat(
        cppMat, lib_std, pred_std, static_cast<size_t>(std::abs(k)), method, include_self);

    // Return nb object (List in R side)
    return std2nb(neighbours);
}

// Wrapper function to compute the nearest neighbours for an input distance matrix
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppNN4DistMat(
    const Rcpp::NumericMatrix& distmat,
    int k,
    bool include_self = false
) {
    // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
    int numRows = distmat.nrow();
    int numCols = distmat.ncol();
    std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            cppMat[r][c] = distmat(r, c);
        }
    }

    // Call the neighbpurbood function
    std::vector<std::vector<size_t>> neighbours = NN::NN4DistMat(
        cppMat, static_cast<size_t>(std::abs(k)), include_self);

    // Return nb object (List in R side)
    return std2nb(neighbours);
}

// Wrapper function to compute the nearest neighbours for an input distance matrix subset
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppNN4DistMatSub(
    const Rcpp::NumericMatrix& distmat,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    int k,
    bool include_self = false
) {
    // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
    int numRows = distmat.nrow();
    int numCols = distmat.ncol();
    std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            cppMat[r][c] = distmat(r, c);
        }
    }

    // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
    std::vector<size_t> lib_std;
    lib_std.reserve(lib.size());
    for (int i = 0; i < lib.size(); ++i) {
        if (lib[i] < 1 || lib[i] > numRows) {
            Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
        }
        lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }

    std::vector<size_t> pred_std;
    pred_std.reserve(pred.size());
    for (int i = 0; i < pred.size(); ++i) {
        if (pred[i] < 1 || pred[i] > numRows) {
            Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
        }
        pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }

    // Call the neighbpurbood function
    std::vector<std::vector<size_t>> neighbours = NN::NN4DistMat(
        cppMat, lib_std, pred_std, static_cast<size_t>(std::abs(k)), include_self);

    // Return nb object (List in R side)
    return std2nb(neighbours);
}
