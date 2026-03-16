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
std::vector<std::vector<size_t>> nb2std(const Rcpp::List& nb) {
  // Get the number of elements in the nb object
  size_t n = static_cast<size_t>(nb.size());
  if (n <= 1) {
    Rcpp::stop("The nb object must contain at least two spatial units (got %d)", n);
  }
  
  // Create a std::vector<std::vector<size_t>> to store the result
  std::vector<std::vector<size_t>> result(n);

  // Iterate over each element in the nb object
  for (size_t i = 0; i < n; ++i) {
    // Get the current element (should be an integer vector)
    Rcpp::IntegerVector current_nb = nb[i];
    size_t cur_num_nb = static_cast<size_t>(current_nb.size());

    // Create a vector<size_t> to store the current subset of elements
    std::vector<size_t> current_subset;
    current_subset.reserve(cur_num_nb);

    // Iterate over each element in the current subset
    for (size_t j = 0; j < cur_num_nb; ++j) {
      // Subtract one from each element to convert from R's 1-based indexing to C++'s 0-based indexing
      current_subset.push_back(current_nb[j] - 1);
    }

    // Add the current subset to the result
    result[i] = current_subset;
  }

  return result;
}

// Function to convert std::vector<std::vector<size_t>> (the `nb` object) to Rcpp::List
Rcpp::List std2nb(const std::vector<std::vector<size_t>>& nb) {
  size_t n = nb.size();
  Rcpp::List result(n);

  for (size_t i = 0; i < n; ++i) {
    const auto& neighbors = nb[i];
    Rcpp::IntegerVector r_neighbors(neighbors.size());
    for (size_t j = 0; j < neighbors.size(); ++j) {
      r_neighbors[j] = static_cast<int>(neighbors[j] + 1);
    }
    result[i] = r_neighbors;
  }

  return result;
}

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
 *          [ r11 ] [ r12 ]
 *          [ r21 ] [ r22 ]
 *          [ r31 ] [ r32 ]
 *
 *      C++:
 *          {
 *              {r1},
 *              {r2},
 *              {r3}
 *          }
 *
 *
 *  When byrow = false
 *
 *      R matrix columns correspond to elements of the outer vector.
 *
 *      R:
 *            c1  c2  c3
 *        r1
 *        r2
 *        r3
 *
 *      C++:
 *          {
 *              {c1},
 *              {c2},
 *              {c3}
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

std::vector<std::vector<double>> mat_r2std(const Rcpp::NumericMatrix& mat,
                                           bool byrow = true)

Rcpp::NumericMatrix mat_std2r(const std::vector<std::vector<double>>& mat,
                              bool byrow = true)