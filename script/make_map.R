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
    map <- raster::mask(map, qc)

    raster::crs(map) <- proj

    names(map) <- paste0(type, "_", year)

    return(map)

}