#' Estimate entropy
#'
#' Estimate entropy for discrete or continuous vector.
#'
#' @param vec The input vector.
#' @param base (optional) Logarithm base for entropy (`2` for bits, `exp(1)` for nats, `10` for dits).
#' @param type (optional) Estimation method: `"disc"` for discrete, `"cont"` for continuous (KSG).
#' @param k (optional) Number of nearest neighbors for continuous entropy estimation (ignored for discrete input).
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
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
