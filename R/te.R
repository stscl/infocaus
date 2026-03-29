#' Transfer Entropy
#'
#' Estimate the transfer entropy from agent variables to target variables.
#'
#' @inheritParams mi
#' @param agent Integer vector of column indices for the source (agent) variables.
#' @param lag_p (optional) Lag of the target variables.
#' @param lag_q (optional) Lag of the agent variables.
#' @param lag_single (optional) Logical; if `FALSE`, use full lag embedding.
#'
#' @returns A numerical value.
#' @export
#'
#' @examples
#' set.seed(42)
#' infoxtr::te(matrix(stats::rnorm(100,1,10),ncol=2),1,2)
#'
te = \(data, target, agent, lag_p = 3, lag_q = 3, base = exp(1), 
       type = c("cont", "disc"), k = 3, normalize = FALSE, lag_single = FALSE) {
  type = match.arg(type)
  mat = as.matrix(data)
  if (type == "disc") {
    return(RcppDiscTE(mat, target, agent, lag_p, lag_q, base, TRUE, normalize, lag_single))
  } else {
    return(RcppContTE(mat, target, agent, lag_p, lag_q, k, 0, base, normalize, lag_single))
  }
}
