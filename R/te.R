te = \(data, target, agent, lag_p = 3, lag_q = 3, base = exp(1), 
       type = c("cont", "disc"), k = 3, normalize = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscTE(mat, target, agent, lag_p, lag_q, base, TRUE, normalize))
  } else {
    return(RcppContTE(mat, target, agent, lag_p, lag_q, k, 0, base, normalize))
  }
}
