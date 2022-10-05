# Binarize maps with lower CI

#---------- SETUP ----------#

library(raster)
library(INLA)
library(ggplot2)
library(rbdi)
library(rangemap)

org <- raster::raster("data/rast.gri")
qc  <- readRDS("data/qc_spacePoly.rds")
rast_qc <- terra::crop(terra::mask(terra::rast(org),
                                   terra::vect(qc)),
                       terra::vect(qc))
rast_qc <- raster::raster(rast_qc)
values(rast_qc)[values(rast_qc) == 1] <- 0

# List of species to work on
species_list <- list.dirs(path()$maps, 
                     recursive = FALSE,
                     full.names = FALSE)


for(species in species_list[1:length(species_list)]) {

  if(!file.exists(paste0(path_sp()$maps,"/maps_pocc.gri"))) next
  maps <- raster::stack(paste0(path_sp()$maps,"/maps_pocc.gri"))
  maps025 <- raster::stack(paste0(path_sp()$maps,"/maps_025.gri"))
  
  if(nlayers(maps) != 27) next

  names(maps) <- paste0(1992:2018)
  for(i in 1:length(names(maps))) {
    
    map <- maps[[i]]
    map025 <- maps025[[i]]

    map[map025<= 0.001] <- 0
    map[map>0] <- 1

    if(raster::cellStats(map, stat = "max") == 0) next

    xy <- raster::coordinates(map)[map[,,] == 1,]
    xy <- as.data.frame(xy[!is.na(xy[,1]),])
    colnames(xy) <- c("Longitude", "Latitude")
    coordinates(xy) <- c("Longitude", "Latitude")
    proj4string(xy) <- raster::crs(map)@projargs
    xy <- sp::spTransform(xy, sp::CRS("EPSG:4326"))
    xy <- cbind(
      data.frame(Species = species),
      as.data.frame(xy)
    )
    
    if(nrow(xy) >= 4) {
      rangemap <- rangemap_hull(xy,
                                hull_type = "concave", 
                                concave_distance_lim = 350000, 
                                buffer_distance = 10000,
                                final_projection = raster::crs(map)@projargs,
                                extent_of_occurrence = FALSE,
                                area_of_occupancy = FALSE)

      rangemap_qc <- terra::intersect(
                       terra::vect(terra::aggregate(rangemap@species_range)),
                       terra::vect(qc))

      if(terra::geomtype(rangemap_qc) != "none") {
        rangemap_qc <- as(rangemap_qc, "Spatial")

        rangemap_qc <- raster::crop(
          raster::rasterize(rangemap_qc,
                            rast_qc, 
                            mask = TRUE),
          rangemap_qc)
        rangemap_qc[rangemap_qc == 0] <- 1

        rangemap_qc <- raster::merge(rangemap_qc, rast_qc)
        names(rangemap_qc) <- paste0(c(1992:2018)[i])

        if(exists("sdms_range")) {
          sdms_range <- raster::stack(sdms_range, rangemap_qc)
        } else {
          sdms_range <- raster::stack(rangemap_qc)
        }
      }
    
      rm(rangemap_qc,
         rangemap,
         map025)
      gc()    
    }

    map <- raster::crop(map, rast_qc)
    if(exists("sdms_occ")) {
      sdms_occ <- raster::stack(sdms_occ, rangemap_qc)
    } else {
      sdms_occ <- raster::stack(rangemap_qc)
    }

    rm(map)
    gc()
  }

  if(exists("sdms_range")) {
    raster::writeRaster(sdms_range, 
                        paste0("output/maps/",species,"/maps_range"),
                        overwrite = TRUE)
  }
  if(exists("sdms_occ")) {
    raster::writeRaster(sdms_occ, 
                        paste0("output/maps/",species,"/maps_occ"),
                        overwrite = TRUE)
  }

  rm(sdms_range)
  rm(sdms_occ)
  gc()

}
