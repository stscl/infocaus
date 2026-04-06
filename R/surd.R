.surd_ts = \(data, target, agent, lag = 1, bin = 5, max.combs = 3, threads = 1){
  if (is.null(bin) || bin <= 0) {
    pfm = RcppDiscMat2PFM(as.matrix(data[,c(target, agents),drop = FALSE]))
  } else {
    pfm = cbind(
      data[,target,drop = TRUE],
      RcppGenTSLagMulti(as.matrix(data[,agents,drop = FALSE]),
                        rep(lag,length.out = length(agents)))
    )
  }
  utils_run_surd(pfm, bin, max.combs, cores)
}

.surd_lattice = \(data, target, agent, lag = 1, bin = 5, max.combs = 3, threads = 1, nb = NULL){
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  data = sf::st_drop_geometry(data)
  if (is.null(bin) || bin <= 0) {
    pfm = RcppDiscMat2PFM(as.matrix(data[,c(target, agents),drop = FALSE]))
  } else {
    pfm = cbind(
      data[,target,drop = TRUE],
      RcppGenLatticeLagMulti(as.matrix(data[,agents,drop = FALSE]),
                             nb,rep(lag,length.out = length(agents)))
    )
  }
  utils_run_surd(pfm, bin, max.combs, cores)
}

.surd_grid = \(data, target, agent, lag = 1, bin = 5, max.combs = 3, threads = 1){
  if (is.null(bin) || bin <= 0) {
    pfm = RcppDiscMat2PFM(terra::values(data[[c(target, agents)]],mat = TRUE,na.rm = FALSE))
  } else {
    pfm = cbind(
      terra::values(data[[target]],mat = TRUE,na.rm = FALSE),
      RcppGenGridLagMulti(terra::values(data[[agents]],mat = TRUE,na.rm = FALSE),
                          rep(lag,length.out = length(agents)),terra::nrow(data))
    )
  }
  utils_run_surd(pfm, bin, max.combs, cores)
}

#' SURD
#' 
#' Synergistic-Unique-Redundant Decomposition of causality
#'
#' @inheritParams te
#' @param lag (optional) Lag of the agent variables.
#' @param bin (optional) Number of discretization bins.
#' @param max.combs (optional) Maximum combination order.
#' @param threads (optional) Number of threads used.
#' @param nb (optional) Neighbours list.
#'
#' @return A list.
#' \describe{
#'   \item{unique}{Unique information contributions per variable.}
#'   \item{synergistic}{Synergistic information components by agent combinations.}
#'   \item{redundant}{Redundant information shared by agent subsets.}
#'   \item{mutual_info}{Mutual information measures for each combination.}
#'   \item{info_leak}{Information leak ratio.}
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
#' \donttest{
#' tryCatch(
#'   infoxtr::surd(columbus, "hoval", c("inc", "crime")),
#'   error = \(e) message("Skipping Python-dependent example: ", e$message)
#' )
#' }
methods::setMethod("surd", "data.frame", .surd_ts)

#' @rdname surd
methods::setMethod("surd", "sf", .surd_lattice)

#' @rdname surd
methods::setMethod("surd", "SpatRaster", .surd_grid)