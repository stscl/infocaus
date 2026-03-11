#ifndef DataTrans_H
#define DataTrans_H

#include <vector>
#include <cstdint>
#include <string>
#include <limits>
#include <numeric>
#include <algorithm>
#include <unordered_map> 
#include <Rcpp.h>

/********************************************************************
 *
 *  Spatial Neighbour Structure Conversion Utilities
 *
 *  These functions convert between:
 *
 *      R representation:
 *          Rcpp::List
 *          Each element is an IntegerVector
 *          Indices are 1 based (R convention)
 *
 *      C++ representation:
 *          std::vector<std::vector<size_t>>
 *          Indices are 0 based (C++ convention)
 *
 *  This structure corresponds to the common "nb" object used in
 *  spatial statistics to represent adjacency lists.
 *
 *  Example in R:
 *
 *      nb[[1]] = c(2, 3)
 *      nb[[2]] = c(1)
 *      nb[[3]] = c(1)
 *
 *  Meaning:
 *      Spatial unit 1 is neighbor with 2 and 3
 *      Spatial unit 2 is neighbor with 1
 *      Spatial unit 3 is neighbor with 1
 *
 *  Conversion rules:
 *
 *      R → C++
 *          - Convert 1 based indices to 0 based
 *          - Store in std::vector<std::vector<size_t>>
 *
 *      C++ → R
 *          - Convert 0 based indices to 1 based
 *          - Return Rcpp::List of IntegerVector
 *
 *  Assumptions:
 *
 *      - Input R list must contain at least two spatial units
 *      - Each element of the list must be an IntegerVector
 *      - No structural validation of symmetry is performed
 *
 ********************************************************************/

// Function to convert Rcpp::List to std::vector<std::vector<size_t>> (the `nb` object)
std::vector<std::vector<size_t>> nb2std(const Rcpp::List& nb);

// Function to convert std::vector<std::vector<size_t>> (the `nb` object) to Rcpp::List
Rcpp::List std2nb(const std::vector<std::vector<size_t>>& nb);

/********************************************************************
 *
 *  Matrix Conversion Utilities (R <-> C++)
 *
 *  These functions convert between:
 *
 *      R representation:
 *          Rcpp::NumericMatrix
 *
 *      C++ representation:
 *          std::vector<std::vector<double>>
 *
 *  The orientation of the conversion is controlled by the
 *  `byrow` argument.
 *
 *  When byrow = true
 *
 *      R matrix rows correspond to elements of the outer vector.
 *
 *      R:
 *          [ r11, r12 ]
 *          [ r21, r22 ]
 *          [ r31, r32 ]
 *
 *      C++:
 *          {
 *              {r11, r12},
 *              {r21, r22},
 *              {r31, r32}
 *          }
 *
 *
 *  When byrow = false
 *
 *      R matrix columns correspond to elements of the outer vector.
 *
 *      R:
 *          [ r11, r12 ]
 *          [ r21, r22 ]
 *          [ r31, r32 ]
 *
 *      C++:
 *          {
 *              {r11, r21, r31},
 *              {r12, r22, r32}
 *          }
 *
 *
 *  Notes
 *
 *      - No copying beyond necessary allocation is performed.
 *      - Matrix dimensions must be non zero.
 *      - No assumption is made about missing values (NA are kept).
 *
 ********************************************************************/

// Function to convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
std::vector<std::vector<double>> mat_r2std(
    const Rcpp::NumericMatrix& mat,
    bool byrow = true
);

// Function to convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
Rcpp::NumericMatrix mat_std2r(
    const std::vector<std::vector<double>>& mat,
    bool byrow = true
);

#endif // DataTrans_H
