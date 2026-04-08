.check_vec1d = \(vec, type = "cont", contain_type = TRUE) {
  if (type == "cont") {
    if (!is.numeric(vec)) {
      type_suffix = if (contain_type) " (`type = \"cont\"`)" else ""
      msg = paste0(
        "The input vector must be numeric for continuous variables ",
        type_suffix, ". "
      )
      stop(msg, call. = FALSE)
      stop("The input vector must be numeric for continuous variables (`type = \"cont\"`). ")
    }
  }
  
  return(vec)
}

.convert2mat = \(data, type = "cont", contain_type = TRUE) {
  if (inherits(data, "sf")) {
    mat = as.matrix(sf::st_drop_geometry(data))
  } else if (inherits(data, "SpatRaster")) {
    mat = terra::values(data, mat = TRUE)
  } else {
    mat = as.matrix(data)
  }

  if (type == "cont" && !(typeof(mat) %in% c("integer", "double"))) {
    type_suffix = if (contain_type) " when `type = \"cont\"`" else ""
    msg = paste0(
      "Non-numeric values detected in input data. ",
      "All variables must be numeric",
      type_suffix,
      ". Please remove columns such as dates, characters, or factors."
    )
    stop(msg, call. = FALSE)
  }

  return(mat)
}
