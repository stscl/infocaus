#' Discretization
#'
#' Discretize a numeric vector into categorical classes using several
#' commonly used discretization methods. Missing values (`NA`/`NaN`) 
#' are ignored and returned as class `0`.
#'
#' @note If `x` is not numeric, it will be converted to
#' integer categories via `as.factor()`.
#'
#' @param x A vector.
#' @param n (optional) Number of classes.
#' @param method (optional) Discretization method. One of
#'   `"sd"`, `"equal"`, `"geometric"`, `"quantile"`,
#'   `"manual"`, `"natural("jenks")"`, or `"headtail"("headtails")`.
#' @param large (optional) Threshold sample size for natural breaks sampling.
#' @param prop (optional) Sampling proportion used when `method = "natural"`
#'   and the input size exceeds `large`.
#' @param seed (optional) Random seed used for sampling in natural breaks.
#' @param thr (optional) Threshold used in the head/tail breaks algorithm.
#' @param iter (optional) Maximum number of iterations for head/tail breaks.
#' @param bps (optional) Numeric vector of manual breakpoints used when
#'   `method = "manual"`.
#' @param right_closed (optional) Logical. If `TRUE`, intervals are right-closed
#'   (e.g., `(a, b]`). If `FALSE`, intervals are left-closed `[a, b)`.
#'
#' @returns A discretized integer vector.
#' @export
#'
#' @examples
#' set.seed(42)
#' infoxtr::discretize(stats::rnorm(99,1,10))
#' 
discretize = \(x, n = 5, method = "natural", large = 3000, prop = 0.15,
               seed = 123456789, thr = 0.4, iter = 100, bps = NULL, right_closed = TRUE){
  if (!is.numeric(x)){
    return(as.integer(as.factor(x)))
  }

  return(RcppDisc(x,n,method,large,prop,seed,thr,iter,bps,right_closed))
}
