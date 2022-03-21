#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

make_map_gam <- function(mod, rast) {

  new_coord <- as.data.frame(raster::coordinates(rast))
  colnames(new_coord) <- c("lon", "lat")

  pred <- predict(mod, new_coord, type = "response")
  new_pred <- cbind(new_coord, pred)
  colnames(new_pred)[3] <- "prob_occ"

  map <- rasterFromXYZ(new_pred)

  return(map)

}