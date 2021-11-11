#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

aggregate_obs <- function(obs, years, i, rast) {

  # Rasterize to aggregate by cells
  obs_aggr <- raster::rasterize(obs[obs$year_obs %in% (years+i),], rast, field = "occurrence", fun = "max")

  # From raster to spatialPoints
  obs_aggr <- raster::rasterToPoints(obs_aggr, spatial = TRUE)
  names(obs_aggr) <- "occurrence"

  # Not sure it's necessary anymore...
  #obs_aggr <- sp::SpatialPointsDataFrame(obs_aggr[,1:2], as.data.frame(obs_aggr[,3]), proj4string = obs@proj4string)
  
  return(obs_aggr)

}