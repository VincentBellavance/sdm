#' Crop and mask a raster with a polygon
#' 
#' Description
#'
#' @param x Raster*
#' @param y Spatial polygon
#'
#' @return Raster

mask_keep_partial <- function(x, y, qc) {

  if(class(x) == "RasterStack" | class(x) == "RasterBrick") {
    spoly_rast <- raster::rasterize(y, x[[1]], getCover = TRUE)
  } else {
    spoly_rast <- raster::rasterize(y, x, getCover = TRUE)
  }
  spoly_rast[spoly_rast == 0] <- NA
  masked_rast <- mask(x, spoly_rast)
  return(crop(masked_rast, y, snap = "out"))

}