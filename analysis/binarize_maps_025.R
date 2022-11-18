# Binarize maps with lower CI

#---------- SETUP ----------#

library(raster)
library(INLA)
library(ggplot2)
library(rangemap)

args <- commandArgs(trailingOnly=TRUE)
species <- args[1]
zone <- args[2]
threshold <- as.numeric(args[3])
years <- 1992:2018

source("R/path.R")

org <- raster::raster("data/rast.gri")
qc <- readRDS(paste0("data/",zone,"/qc_spacePoly.rds"))
rast_qc <- terra::crop(terra::mask(terra::rast(org),
                                   terra::vect(qc)),
                       terra::vect(qc))
rast_qc <- raster::raster(rast_qc)
values(rast_qc)[values(rast_qc) == 1] <- 0

# Import spatial opbjects for the species
study_extent <- readRDS(paste0(path_sp(species, zone)$spat, "/study_extent.rds"))
rast <- raster::raster(paste0(path_sp(species, zone)$spat, "/rast.gri"))
qc <- readRDS(paste0("data/",zone,"/qc_spacePoly.rds"))
mesh <- readRDS(paste0(path_sp(species, zone)$spat, "/mesh.rds"))

# maps of 2.5e quantile
maps025 <- stack(paste0(path_sp(species, zone)$maps, "/maps_025.gri"))

#---------- CREATE BINARY MAP FOR EVERY YEAR ----------#

binary_maps <- maps025 >= threshold
qcr <- raster::rasterize(qc, binary_maps, getCover = TRUE)
qcr <- qcr[qcr == 0] <- NA
binary_maps_final <- raster::crop(raster::mask(binary_maps, qcr), qcr)
binary_maps_final <- raster::merge(binary_maps_final, rast_qc, overlap = FALSE)
names(binary_maps_final) <- years

# Save occurrence maps
raster::writeRaster(binary_maps_final, 
                    paste0(path_sp(species, zone)$maps,"/maps_occ_025"),
                    overwrite = TRUE)



#---------- MAKE MAPS OF SPECIES RANGE ----------#

# Crop and mask the rasters with study_extent to make rangemaps 
ser <- raster::rasterize(study_extent, binary_maps, getCover = TRUE)
ser <- ser[ser == 0] <- NA
binary_maps <- raster::crop(raster::mask(binary_maps, ser), ser)

for(i in 1:length(names(binary_maps))) {
    
  map <- binary_maps[[i]]

  # Remove small clumps (3 or less cells together)
  map_clump <- raster::clump(map, directions = 4)
  f <- as.data.frame(freq(map_clump))
  excludeID <- f$value[which(f$count <= 4)]
  map_clump[map_clump %in% excludeID] <- NA
  map_clump[!is.na(values(map_clump))] <- 1 # Transform ID in presence
  values(map)[!is.na(values(map))] <- 0 # Basic map of 0 to merge with map_clump

  map <- raster::merge(map_clump, map)

  # If there is no presence at all
  if(raster::cellStats(map, stat = "max") == 0) {

    rangemap_qc <- rast_qc
    names(rangemap_qc) <- paste0(c(1992:2018)[i])

    if(exists("sdms_range")) {
      sdms_range <- raster::stack(sdms_range, rangemap_qc)
    } else {
      sdms_range <- raster::stack(rangemap_qc)
    }   
    
    rm(rangemap_qc)

    next
  }

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
  
  # Make range maps if there is more than 4 points
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

    # If the range intersect with the Quebec polygon
    if(terra::geomtype(rangemap_qc) != "none") {
      rangemap_qc <- as(rangemap_qc, "Spatial")

      rangemap_qc <- raster::mask(
                       raster::crop(
                         raster::rasterize(rangemap_qc,
                                           rast_qc, 
                                           mask = TRUE),
                         rangemap_qc),
                       qc)
      rangemap_qc[rangemap_qc == 0] <- 1
      rangemap_qc[is.na(rangemap_qc)] <- 0

      rangemap_qc <- raster::merge(rangemap_qc, rast_qc, overlap = FALSE)
      rangemap_qc <- raster::mask(rangemap_qc, qc)
      names(rangemap_qc) <- paste0(c(1992:2018)[i])

      if(exists("sdms_range")) {
        sdms_range <- raster::stack(sdms_range, rangemap_qc)
      } else {
        sdms_range <- raster::stack(rangemap_qc)
      }

    } else {
      
      # If the range DOES NOT intersect with the Quebec polygon
      rangemap_qc <- rast_qc  
      names(rangemap_qc) <- paste0(c(1992:2018)[i])

      if(exists("sdms_range")) {
        sdms_range <- raster::stack(sdms_range, rangemap_qc)
      } else {
        sdms_range <- raster::stack(rangemap_qc)
      }   

    }
    
    rm(rangemap_qc,
       rangemap) 
  
  } else {

    # If there is not enough points to make a range polygon

    rangemap_qc <- raster::crop(raster::mask(map, qc),qc)
    rangemap_qc <- raster::merge(rangemap_qc, rast_qc, overlap = FALSE)
    rangemap_qc <- raster::mask(rangemap_qc, qc)

    if(exists("sdms_range")) {
      sdms_range <- raster::stack(sdms_range, rangemap_qc)
    } else {
      sdms_range <- raster::stack(rangemap_qc)
    }
    rm(rangemap_qc,
       rangemap)
  }
}

if(exists("sdms_range")) {
  raster::writeRaster(sdms_range, 
                      paste0(path_sp(species, zone)$maps,"/maps_range_025"),
                      overwrite = TRUE)
}

rm(sdms_range)
gc()
