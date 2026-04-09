.surd_ts = \(data, target, agent, lag = 1, bin = 5, method = "equal",
             max.combs = 10, threads = 1, base = 2.0, normalize = TRUE) {
  mat = .convert2mat(data, contain_type = FALSE)
  return(RcppSURD(mat, target, agent, lag, bin, max.combs, 
                  threads, base, normalize, method))
}

.surd_lattice = \(data, target, agent, lag = 1, bin = 5, method = "equal", 
                  max.combs = 10, threads = 1, base = 2.0, normalize = TRUE, nb = NULL) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  mat = .convert2mat(data, contain_type = FALSE)
  return(RcppSURD(mat, target, agent, lag, bin, max.combs, 
                  threads, base, normalize, method, nb))
}

.surd_grid = \(data, target, agent, lag = 1, bin = 5, method = "equal",
               max.combs = 10, threads = 1, base = 2.0, normalize = TRUE) {
  mat = .convert2mat(data, contain_type = FALSE)
  return(RcppSURD(mat, target, agent, lag, bin, max.combs, 
                  threads, base, normalize, method, NULL, terra::nrow(data[[1]])))
}

#' SURD
#' 
#' Synergistic-Unique-Redundant Decomposition
#' 
#' @note SURD only support numeric data.
#'
#' @inheritParams te
#' @param lag (optional) Lag of the agent variables.
#' @param bin (optional) Number of discretization bins.
#' @param method (optional) Discretization method. One of
#'   `"sd"`, `"equal"`, `"geometric"`, `"quantile"`,
#'   `"natural("jenks")"`, or `"headtail"("headtails")`.
#' @param max.combs (optional) Maximum combination order.
#' @param threads (optional) Number of threads used.
#' @param nb (optional) Neighbours list.
#'
#' @return A list.
#' \describe{
#'   \item{vars}{Character vector indicating the variable combination associated with each information component.}
#'   \item{types}{Character vector indicating the information type of each component.}
#'   \item{values}{Numeric vector giving the magnitude of each information component.}
#' }
#'
#' @export
#' @name surd
#' @aliases surd,data.frame-method
#' @references
#' Martinez-Sanchez, A., Arranz, G. & Lozano-Duran, A. Decomposing causality into its synergistic, unique, and redundant components. Nat Commun 15, 9296 (2024).
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' infoxtr::surd(columbus, 1, 2:3)
#'
methods::setMethod("surd", "data.frame", .surd_ts)

#' @rdname surd
methods::setMethod("surd", "sf", .surd_lattice)

#' @rdname surd
methods::setMethod("surd", "SpatRaster", .surd_grid)
