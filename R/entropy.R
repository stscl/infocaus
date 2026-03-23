je = \(data, indices, base = 2, type = c("disc", "cont"), k = 3) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscJE(mat, indices, base, TRUE))
  } else {
    return(RcppContJE(mat, indices, k, 0, base))
  }
}
