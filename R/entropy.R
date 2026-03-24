#' Entropy Estimation
#'
#' Estimate the entropy of a vector using either discrete counts
#' or a k-nearest neighbor estimator for continuous data.
#'
#' @param vec The input vector.
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
#' entropy(rnorm(100), type = "cont")
#' entropy(sample(letters[1:5], 100, TRUE), base = 2, type = "disc")
#'
entropy = \(vec, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  if (type == "disc") {
    return(RcppDiscEntropy(vec, base, TRUE))
  } else {
    return(RcppContEntropy(vec, k, 0, base))
  }
}

je = \(data, indices, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscJE(mat, indices, base, TRUE))
  } else {
    return(RcppContJE(mat, indices, k, 0, base))
  }
}

ce = \(data, target, conds, base = exp(1), type = c("cont", "disc"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscCE(mat, target, conds, base, TRUE))
  } else {
    return(RcppContCE(mat, target, conds, k, 0, base))
  }
}

mi = \(data, target, interact, base = exp(1), type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscMI(mat, target, interact, base, TRUE, normalize))
  } else {
    return(RcppContMI(mat, target, interact, k, 0, base, normalize))
  }
}

cmi = \(data, target, interact, conds, base = exp(1), type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscCMI(mat, target, interact, conds, base, TRUE, normalize))
  } else {
    return(RcppContCMI(mat, target, interact, conds, k, 0, base, normalize))
  }
}
