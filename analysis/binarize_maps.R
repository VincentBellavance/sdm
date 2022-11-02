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

# List models and stacks
models <- list.files(path_sp(species, zone)$mod)
Stacks <- list.files(path_sp(species, zone)$stack)
years <- as.integer(gsub(".rds", "", models))


#---------- CREATE BINARY MAP FOR EVERY YEAR ----------#

binary_maps <- lapply(1:length(models), function(x) {

  cat(paste0("Sampling for year ",years[x]," in progress : ",x,"/",length(models)), "\n")

  #---------- POSTERIOR SAMPLINGS ----------#
  
  # Import model
  mod <- readRDS(paste0(path_sp(species, zone)$mod, "/", models[x]))
  Stack <- readRDS(paste0(path_sp(species, zone)$stack, "/", Stacks[x]))


  # Make map for entire sPoly to compute AUC
  # mapBasis
  mapBasis <- inla.mesh.projector(mesh,
                                  dims = dim(rast)[2:1],
                                  xlim = c(xmin(rast), 
                                           xmax(rast)),
                                  ylim = c(ymin(rast), 
                                           ymax(rast)),
                                  crs = mesh$crs)

  ### Find the mesh edges on which predictions should be made
  ID <- inla.stack.index(Stack, tag="pred")$data
  pred <- suppressMessages(
            inla.posterior.sample(100, 
                                  result = mod, 
                                  selection = list(APredictor = 0,
                                                   Predictor = 0),
                                  intern = TRUE,
                                  use.improved.mean = FALSE))

  maps_list <- lapply(1:length(pred), function(y) {
    
    fitted <- exp(pred[[y]]$latent[ID,])/(1 + exp(pred[[y]]$latent[ID,]))
    mapPred <- inla.mesh.project(mapBasis, fitted)
    mapRaster <- raster::raster(t(mapPred[,ncol(mapPred):1]),
                                  xmn = min(mapBasis$x), 
                                  xmx = max(mapBasis$x), 
                                  ymn = min(mapBasis$y), 
                                  ymx = max(mapBasis$y),
                                  crs = mesh$crs)
    names(mapRaster) <- paste0(years[x])
    cat("\r", paste0("... Iteration ",y," of ",length(pred)," done!"))
    if(y == length(pred)) cat("\n")
    return(mapRaster)
  })      
    
  maps <- do.call(raster::stack, maps_list)


  #---------- BINARIZE MAPS ----------#

  maps_df <- maps[1:raster::ncell(maps)]

  binary_values <- sapply(1:nrow(maps_df), function(i) {
    distr <- maps_df[i,]
    if(all(is.na(distr))) return(NA)
    quant <- quantile(distr, probs=0.025, na.rm = TRUE)
    if(quant<threshold) {
      return(0)
    } else {
      return(1)
    }
  })

  binary_map <- raster(crs = raster::crs(maps),
                       ext = raster::extent(maps),
                       resolution = round(raster::res(maps)),
                       vals = binary_values)

  return(binary_map)

  cat(paste0("Sampling for year ",years[x]," done! : ",x,"/",length(models)), "\n")
})

# Gather maps in a rasterStack and crop it with Qc sPoly
binary_maps <- do.call(raster::stack, binary_maps)
names(binary_maps) <- years

binary_maps_final <- terra::crop(
                       terra::mask(terra::rast(binary_maps), 
                                   terra::vect(qc)),
                       terra::vect(qc))
binary_maps_final <- raster::stack(binary_maps_final)
binary_maps_final <- raster::merge(binary_maps_final, rast_qc, overlap = FALSE)
names(binary_maps_final) <- years

# Save occurrence maps
raster::writeRaster(binary_maps_final, 
                    paste0(path_sp(species, zone)$maps,"/maps_occ"),
                    overwrite = TRUE)



#---------- MAKE MAPS OF SPECIES RANGE ----------#

for(i in 1:length(names(binary_maps))) {
    
  map <- binary_maps[[i]]
  
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
                      paste0(path_sp(species)$maps,"/maps_range"),
                      overwrite = TRUE)
}

rm(sdms_range)
gc()
