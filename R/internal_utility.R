.check_vec1d = \(vec, type = "cont") {
  if (type == "cont") {
    if (!is.numeric(vec)) {
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
    type_clause = if (contain_type) " When `type = \"cont\"`," else ""
    msg = paste0(
      "Non-numeric values detected in input data.",
      type_clause,
      " All variables must be numeric. ",
      "Please remove columns such as dates, characters, or factors."
    )
    stop(msg, call. = FALSE)
    stop(
      "Non-numeric values detected in input data. When `type = \"cont\"`, ",
      "all variables must be numeric. Please remove columns such as dates, ",
      "characters, or factors."
    )
  }

  return(mat)
}
