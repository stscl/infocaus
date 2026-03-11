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
) {

  size_t nrow = static_cast<size_t>(mat.nrow());
  size_t ncol = static_cast<size_t>(mat.ncol());

  if (nrow == 0 || ncol == 0) {
    Rcpp::stop("Input matrix must have positive dimensions.");
  }

  std::vector<std::vector<double>> result;

  if (byrow) {

    // Each row becomes one vector
    result.resize(nrow);

    for (size_t i = 0; i < nrow; ++i) {
      result[i].reserve(ncol);

      for (size_t j = 0; j < ncol; ++j) {
        result[i].push_back(mat(i, j));
      }
    }

  } else {

    // Each column becomes one vector
    result.resize(ncol);

    for (size_t j = 0; j < ncol; ++j) {
      result[j].reserve(nrow);

      for (size_t i = 0; i < nrow; ++i) {
        result[j].push_back(mat(i, j));
      }
    }

  }

  return result;
}

// Function to convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
Rcpp::NumericMatrix mat_std2r(
    const std::vector<std::vector<double>>& mat,
    bool byrow = true
) {

  size_t outer = mat.size();

  if (outer == 0) {
    Rcpp::stop("Input matrix container is empty.");
  }

  size_t inner = mat[0].size();

  if (inner == 0) {
    Rcpp::stop("Matrix rows/columns must contain elements.");
  }

  if (byrow) {

    // std rows -> R rows
    Rcpp::NumericMatrix result(outer, inner);

    for (size_t i = 0; i < outer; ++i) {

    //   if (mat[i].size() != inner) {
    //     Rcpp::stop("Inconsistent row length in input matrix.");
    //   }

      for (size_t j = 0; j < inner; ++j) {
        result(i, j) = mat[i][j];
      }
    }

    return result;

  } else {

    // std rows -> R columns
    Rcpp::NumericMatrix result(inner, outer);

    for (size_t j = 0; j < outer; ++j) {

    //   if (mat[j].size() != inner) {
    //     Rcpp::stop("Inconsistent column length in input matrix.");
    //   }

      for (size_t i = 0; i < inner; ++i) {
        result(i, j) = mat[j][i];
      }
    }

    return result;
  }
}

/********************************************************************
 *  pat_r2std
 *
 *  Convert an R matrix (Integer / Numeric / Character)
 *  into std::vector<std::vector<uint64_t>>.
 *
 *  Structure:
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
 *  Design:
 *      - Scan matrix once to collect global unique values
 *      - Sort uniques
 *      - Assign id 1..uniq
 *
 *  NA handling:
 *      NA encoded as {0}
 *
 ********************************************************************/
std::vector<std::vector<uint64_t>> pat_r2std(SEXP x, bool byrow = true)
{
    if (!Rf_isMatrix(x))
        Rcpp::stop("Input must be a matrix.");

    std::vector<std::vector<uint64_t>> mat;

    switch (TYPEOF(x))
    {

    /******************************
     * IntegerMatrix
     ******************************/
    case INTSXP:
    {
        Rcpp::IntegerMatrix m(x);

        const size_t n = static_cast<size_t>(m.nrow());
        const size_t p = static_cast<size_t>(m.ncol());

        if (byrow) 
        {
          mat.resize(n);
          for (size_t j = 0; j < n; ++j)
              mat[j].reserve(p);  
        }
        else 
        {
          mat.resize(p);
          for (size_t j = 0; j < p; ++j)
              mat[j].reserve(n);        
        }

        std::vector<int> uniq;
        uniq.reserve(n*p);

        for (size_t j = 0; j < p; ++j)
        {
            for (size_t i = 0; i < n; ++i)
            {
                int v = m(i,j);
                if (!Rcpp::IntegerVector::is_na(v))
                    uniq.push_back(v);
            }
        }

        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

        std::unordered_map<int,uint64_t> dict;
        dict.reserve(uniq.size());

        for (uint64_t i = 0; i < uniq.size(); ++i)
            dict[uniq[i]] = i+1;

        if (byrow)
        {
          for (size_t i = 0; i < n; ++i)
          {
              auto &row = mat[i];

              for (size_t j = 0; j < p; ++j)
              {
                  int v = m(i,j);

                  if (Rcpp::IntegerVector::is_na(v))
                      row.push_back({0});
                  else
                      row.push_back(dict[v]);
              }
          }
        }
        else 
        {
          for (size_t j = 0; j < p; ++j)
          {
              auto &col = mat[j];

              for (size_t i = 0; i < n; ++i)
              {
                  int v = m(i,j);

                  if (Rcpp::IntegerVector::is_na(v))
                      col.push_back({0});
                  else
                      col.push_back(dict[v]);
              }
          }
        }

        break;
    }

    /******************************
     * NumericMatrix
     ******************************/
    case REALSXP:
    {
        Rcpp::NumericMatrix m(x);

        const size_t n = static_cast<size_t>(m.nrow());
        const size_t p = static_cast<size_t>(m.ncol());

        if (byrow) 
        {
          mat.resize(n);
          for (size_t j = 0; j < n; ++j)
              mat[j].reserve(p);  
        }
        else 
        {
          mat.resize(p);
          for (size_t j = 0; j < p; ++j)
              mat[j].reserve(n);        
        }

        std::vector<double> uniq;
        uniq.reserve(n*p);

        for (size_t j = 0; j < p; ++j)
        {
            for (size_t i = 0; i < n; ++i)
            {
                double v = m(i,j);
                if (!Rcpp::NumericVector::is_na(v))
                    uniq.push_back(v);
            }
        }

        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

        std::unordered_map<double,uint64_t> dict;
        dict.reserve(uniq.size());

        for (uint64_t i = 0; i<uniq.size(); ++i)
            dict[uniq[i]] = i+1;

        if (byrow)
        {
          for (size_t i = 0; i < n; ++i)
          {
              auto &row = mat[i];

              for (size_t j = 0; j < p; ++j)
              {
                  double v = m(i,j);

                  if (Rcpp::NumericVector::is_na(v))
                      row.push_back({0});
                  else
                      row.push_back(dict[v]);
              }
          }
        }
        else 
        {
          for (size_t j = 0; j < p; ++j)
          {
              auto &col = mat[j];

              for (size_t i = 0; i < n; ++i)
              {
                  double v = m(i,j);

                  if (Rcpp::NumericVector::is_na(v))
                      col.push_back({0});
                  else
                      col.push_back(dict[v]);
              }
          }
        }

        break;
    }

    /******************************
     * CharacterMatrix
     ******************************/
    case STRSXP:
    {
        Rcpp::CharacterMatrix m(x);

        const size_t n = static_cast<size_t>(m.nrow());
        const size_t p = static_cast<size_t>(m.ncol());

        if (byrow) 
        {
          mat.resize(n);
          for (size_t j = 0; j < n; ++j)
              mat[j].reserve(p);  
        }
        else 
        {
          mat.resize(p);
          for (size_t j = 0; j < p; ++j)
              mat[j].reserve(n);        
        }

        std::vector<std::string> uniq;
        uniq.reserve(n*p);

        for (size_t j = 0; j < p; ++j)
        {
            for (size_t i = 0; i < n; ++i)
            {
                if (!Rcpp::CharacterVector::is_na(m(i,j)))
                    uniq.emplace_back(std::string(m(i,j)));
            }
        }

        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

        std::unordered_map<std::string,uint64_t> dict;
        dict.reserve(uniq.size());

        for (uint64_t i = 0; i < uniq.size(); ++i)
            dict[uniq[i]] = i+1;

        if (byrow)
        {
          for (size_t i = 0; i < n; ++i)
          {
              auto &row = mat[i];

              for (size_t j = 0; j < p; ++j)
              {
                  if (Rcpp::CharacterVector::is_na(m(i,j)))
                  {
                      row.push_back({0});
                  }
                  else
                  {
                      std::string key = std::string(m(i,j));
                      row.push_back(dict[key]);
                  }
              }
          }
        }
        else 
        {
          for (size_t j = 0; j < p; ++j)
          {
              auto &col = mat[j];

              for (size_t i = 0; i < n; ++i)
              {
                  if (Rcpp::CharacterVector::is_na(m(i,j)))
                  {
                      col.push_back({0});
                  }
                  else
                  {
                      std::string key = std::string(m(i,j));
                      col.push_back(dict[key]);
                  }
              }
          }
        }

        break;
    }

    default:
        Rcpp::stop("Matrix must be Integer, Numeric, or Character.");
    }

    return mat;
}
