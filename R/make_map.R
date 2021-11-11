#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

make_map <- function(type, mod, rast, sPoly, year, proj) {

  map_type <- ifelse(type == "mean_all", "mean", type)
  map <- mapSpecies::mapSpace(mod,
                              dims = dim(rast)[2:1],
                              type = map_type,
                              sPoly = sPoly)

  raster::crs(map) <- raster::crs(sPoly)

  map <- terra::mask(terra::rast(map), terra::vect(sPoly))

  map <- raster::raster(map)

  raster::crs(map) <- raster::crs(sPoly)

  names(map) <- paste0(type, "_", year)

  return(map)

}