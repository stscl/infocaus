#' KOCMI
#'
#' Knockoff Conditional Mutual Information
#'
#' @inheritParams surd
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
##' @param k (optional) For `type = "cont"`, the number of nearest neighbors used by
#'   the continuous conditional mutual information estimator.
#'   For `type = "disc"`, the number of bins used for discretization.
#' @param seed (optional) Random seed used for permutation test.
#' @param contain_null (optional) Logical. If `TRUE`, the test statistic is
#'   computed using knockoffs generated under the null model (provided in
#'   `null_knockoff`). In this case the difference is defined as
#'   \eqn{I(Y; X_{null} | Z) - I(Y; X_{knockoff} | Z)}.
#'   If `FALSE`, the original conditional mutual information
#'   \eqn{I(Y; X | Z)} is used instead and compared against the knockoff
#'   estimates \eqn{I(Y; X_{knockoff} | Z)}.
#'
#' @returns A named numeric vector.
#' @export
#' @references 
#' Zhang, X., Chen, L., 2025. Quantifying interventional causality by knockoff operation. Science Advances 11.
#'
#' @examples
#' popd_df = readr::read_csv(system.file("case/popd.csv",package = "spEDM"))
#' popd_df
#' # infoxtr::kocmi(popd_df, 3, 6, c(4,5,7))
#' 
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
