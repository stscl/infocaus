.check_vec1d = \(vec, type = "cont") {
  if (type == "cont") {
    if (!is.numeric(vec)) {
      stop("The input vector must be numeric for continuous variables (`type = \"cont\"`). ")
    }
  }
  
  return(vec)
}

.convert2mat = \(data, type = "cont") {
  if (inherits(data, "sf")) {
    mat = as.matrix(sf::st_drop_geometry(data))
  } else if (inherits(data, "SpatRaster")) {
    mat = terra::values(data, mat = TRUE)
  } else {
    mat = as.matrix(data)
  }

  if (type == "cont" && !(typeof(mat) %in% c("integer", "double"))) {
    stop(
      "Non-numeric values detected in input data. When `type = \"cont\"`, ",
      "all variables must be numeric. Please remove columns such as dates, ",
      "characters, or factors."
    )
  }

  return(mat)
}
