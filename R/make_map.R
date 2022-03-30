#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

make_map <- function(type, mesh, mod, rast, sPoly, year) {

  mapBasis <- inla.mesh.projector(mesh,
                                  dims = dim(rast)[2:1],
                                  xlim = c(xmin(sPoly), xmax(sPoly)),
                                  ylim = c(ymin(sPoly), ymax(sPoly)),
                                  crs = mesh$crs)

  #mapPred <- inla.mesh.project(mapBasis, 
  #                             mod$summary.random[['i']][['mean']])

  ### Find the mesh edges on which predictions should be made
  ID <- inla.stack.index(Stack, tag="pred")$data
  
  ### Calculate prediction
  mapPred <- inla.mesh.project(mapBasis, 
                               mod$summary.fitted.values[["mean"]][ID])

  mapRaster <- raster::raster(t(mapPred[,ncol(mapPred):1]),
                                xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                                ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                                crs = mesh$crs)

  mapRaster <- terra::mask(terra::rast(mapRaster), terra::vect(sPoly))

  map <- raster::raster(mapRaster)

  names(map) <- paste0(type, "_", year)

  return(map)

}