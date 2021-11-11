#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

extract_coords <- function(geom) {
  
  # regex to use to extract lat and lon
  re <- "\\(([^()]+)\\)"
  coords <- strsplit(
    gsub(re, "\\1", stringr::str_extract_all(geom, re)),
    "\\s+"
  )

  ## Bind the lon and lat to create the df
  coords <- do.call(rbind,
    lapply(coords, function(x) {
      cbind(as.numeric(x[1]), as.numeric(x[2]))
    })
  )
  coords <- as.data.frame(coords)
  colnames(coords) <- c("lon", "lat")

  return(coords)

}