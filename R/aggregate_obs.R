#' Title
#' 
#' Description
#'
#' @param
#'
#' @return

aggregate_obs <- function(obs, years, i, rast) {

  # Rasterize to aggregate by cells
  obs_aggr <- terra::rasterize(terra::vect(obs[obs$year_obs %in% (years+i),]), terra::rast(rast), field = "occurrence", fun = "max")

  # From raster to spatialPoints
  obs_aggr <- as(terra::as.points(obs_aggr), "Spatial")
  names(obs_aggr) <- "occurrence"

  # Not sure it's necessary anymore...
  #obs_aggr <- sp::SpatialPointsDataFrame(obs_aggr[,1:2], as.data.frame(obs_aggr[,3]), proj4string = obs@proj4string)
  
  return(obs_aggr)

}