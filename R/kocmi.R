#' KOCMI
#'
#' Knockoff Conditional Mutual Information
#'
#' @param conds Integer vector of column indices for the conditioning variables.
#' @param knockoff
#' @param null_knockoff
#' @param type (optional) Estimation method:
#'   `"disc"` for discrete entropy or `"cont"` for continuous entropy (KSG estimator).
#' @param nboots
#' @param k
#' @param seed (optional) Random seed used for permutation test.
#' @param contain_null
#'
#' @returns A named numeric vector.
#' @export
#'
#' @examples
kocmi = \(data, target, agent, conds, knockoff, null_knockoff = NULL,
          type = c("cont", "disc"), nboots = 1e4, k = 3, threads = 1,
          seed = 123456789, base = exp(1), method = "equal", contain_null = TRUE) {
  type = match.arg(type)
  mat = .convert2mat(data, type)
  knockoff = .convert2mat(knockoff, contain_type = FALSE)
  null_knockoff = .convert2mat(null_knockoff, contain_type = FALSE)
  return(RcppKOCMI(mat, target, agent, conds, knockoff, null_knockoff, type,
                   nboots, k, 0, threads, seed, base, method, contain_null))
}
