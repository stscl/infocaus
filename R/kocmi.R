#' KOCMI
#'
#' Knockoff Conditional Mutual Information
#'
#' @param conds Integer vector of column indices for the conditioning variables.
#' @param knockoff Knockoff realizations constructed for the `agent` variable
#'   while keeping the `target` variable unchanged. Each column corresponds to
#'   one Monte Carlo knockoff sample generated using the remaining variables
#'   except the target.
#' @param null_knockoff (optional) Knockoff realizations generated under the
#'   null setting where all variables are jointly used to construct knockoffs.
#'   Each column represents one Monte Carlo sample. If `contain_null = FALSE`,
#'   this argument can be `NULL`.
#' @param type (optional) Estimation method: `"disc"` for discrete mutual information
#'   or `"cont"` for continuous mutual information (KSG estimator).
#' @param nboots (optional) Number of permutations used in the sign-flipping permutation
#'   test for evaluating the significance of the mean information difference.
##' @param k (optional) For `type = "cont`, the number of nearest neighbors used by
#'   the continuous conditional mutual information estimator.
#'   For `type = "disc"`, the number of bins used for discretization.
#' @param seed (optional) Random seed used for permutation test.
#' @param contain_null (optional) Logical. If `TRUE`, the test statistic is
#'   computed using knockoffs generated under the null model (provided in
#'   `null_knockoff`). In this case the difference is defined as
#'   \eqn{I(Y; X_null | Z) - I(Y; X_knockoff | Z)}.
#'   If `FALSE`, the original conditional mutual information
#'   \eqn{I(Y; X | Z)} is used instead and compared against the knockoff
#'   estimates \eqn{I(Y; X_knockoff | Z)}.
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
