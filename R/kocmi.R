#' KOCMI
#'
#' Knockoff Conditional Mutual Information
#'
#' @param conds Integer vector of column indices for the conditioning variables.
#' @param knockoff A matrix or data frame containing knockoff realizations
#'   constructed for the `agent` variable while keeping the `target`
#'   variable unchanged. Each row corresponds to one Monte Carlo knockoff
#'   sample generated using the remaining variables except the target.
#' @param null_knockoff (optional) A matrix or data frame containing knockoff
#'   realizations generated under the null setting where all variables are
#'   jointly used to construct knockoffs. Each row represents one Monte Carlo
#'   sample. If `contain_null = FALSE`, this argument can be `NULL`.
#' @param type (optional) Estimation method:
#'   `"disc"` for discrete entropy or `"cont"` for continuous entropy (KSG estimator).
#' @param nboots
##' @param k For `type = "cont`, the number of nearest neighbors used by
#'   the continuous conditional mutual information estimator.
#'   For `type = "disc"`, the number of bins used for discretization.
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
