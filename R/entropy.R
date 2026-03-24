#' Entropy Estimation
#'
#' Estimate the entropy of a vector using either category counts
#' (for discrete data) or a k-nearest neighbor estimator (for continuous data).
#'
#' @param vec A vector.
#' @param base (optional) Logarithm base of the entropy.
#'   Defaults to `exp(1)` (nats). Use `2` for bits or `10` for dits.
#' @param type (optional) Estimation method:
#'   `"disc"` for discrete entropy or `"cont"` for continuous entropy (KSG estimator).
#' @param k (optional) Number of nearest neighbors used by the continuous estimator.
#'   Ignored when `type = "disc"`.
#'
#' @return A numerical value.
#' @export
#'
#' @examples
#' set.seed(42)
#' infocaus::entropy(stats::rnorm(100), type = "cont")
#' infocaus::entropy(sample(letters[1:5], 100, TRUE), base = 2, type = "disc")
#'
entropy = \(vec, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  if (type == "disc") {
    return(RcppDiscEntropy(vec, base, TRUE))
  } else {
    return(RcppContEntropy(vec, k, 0, base))
  }
}

#' Joint Entropy
#'
#' Estimate the joint entropy of selected variables.
#'
#' @inheritParams entropy
#' @param data Observation data.
#' @param indices Integer vector of column indices to include in joint entropy calculation.
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
#' infocaus::je(matrix(1:100,ncol=2),1:2)
#'
je = \(data, indices, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscJE(mat, indices, base, TRUE))
  } else {
    return(RcppContJE(mat, indices, k, 0, base))
  }
}

#' Conditional Entropy
#'
#' Estimate the conditional entropy of target variables given conditioning variables.
#'
#' @inheritParams je
#' @param target Integer vector of column indices for the target variables.
#' @param conds Integer vector of column indices for the conditioning variables.
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
#' infocaus::ce(matrix(1:100,ncol=2),1,2)
#'
ce = \(data, target, conds, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscCE(mat, target, conds, base, TRUE))
  } else {
    return(RcppContCE(mat, target, conds, k, 0, base))
  }
}

#' Mutual Information
#'
#' Estimate the mutual information between target and interacting variables.
#'
#' @inheritParams ce
#' @param interact Integer vector of column indices for the interacting variables.
#' @param normalize (optional) Logical; if `TRUE`, return normalized mutual information.
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
#' infocaus::mi(matrix(1:100,ncol=2),1,2)
#'
mi = \(data, target, interact, base = exp(1), type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscMI(mat, target, interact, base, TRUE, normalize))
  } else {
    return(RcppContMI(mat, target, interact, k, 0, base, normalize))
  }
}

#' Conditional Mutual Information
#'
#' Estimate the conditional mutual information between target and interacting
#' variables given conditioning variables.
#'
#' @inheritParams mi
#' @param conds Integer vector of column indices for the conditioning variables.
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
#' set.seed(42)
#' infocaus::cmi(matrix(stats::rnorm(99,1,10),ncol=3),1,2,3)
#'
cmi = \(data, target, interact, conds, base = exp(1), type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscCMI(mat, target, interact, conds, base, TRUE, normalize))
  } else {
    return(RcppContCMI(mat, target, interact, conds, k, 0, base, normalize))
  }
}
