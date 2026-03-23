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
