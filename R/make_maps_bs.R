#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

make_maps_bs <- function(mod, mesh, rast, sPoly, Stack, year) {

  mapBasis <- inla.mesh.projector(mesh,
                                  dims = dim(rast)[2:1],
                                  xlim = c(xmin(sPoly), xmax(sPoly)),
                                  ylim = c(ymin(sPoly), ymax(sPoly)),
                                  crs = mesh$crs)

  ### Find the mesh edges on which predictions should be made
  ID <- inla.stack.index(Stack, tag="pred")$data
  
  pred <- inla.posterior.sample(1000, result = mod, selection = list(APredictor = 0, Predictor = 0))

  maps <- lapply(pred, function(x) {
    fitted <- exp(x$latent[ID,])/(1 + exp(x$latent[ID,]))
    mapPred <- inla.mesh.project(mapBasis, fitted)
    mapRaster <- raster::raster(t(mapPred[,ncol(mapPred):1]),
                                  xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                                  ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                                  crs = mesh$crs)
    mapRaster <- terra::mask(terra::rast(mapRaster), terra::vect(sPoly))
    map <- raster::raster(mapRaster)
  })

  maps <- do.call(raster::stack, maps)

  names(maps) <- paste0(year, "_", 1:1000)

  return(maps)

}
