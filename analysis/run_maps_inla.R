# Make maps for all sPoly and for Qc only + compute AUC for different threshold

#--------------------------------------------------------------
# Steps:
# 1. Setup environment
# 2. In a loop:
#   2.1: Make map for entire sPoly
#   2.2: Compute AUC with different threshold
#   2.3: Crop and mask entire sPoly to get map for Qc only
#--------------------------------------------------------------

#--- Setup ---#
# Import arguments
args = commandArgs(trailingOnly=TRUE)

# Import mapSpecies
suppressMessages(library(INLA))
suppressMessages(library(raster))
suppressMessages(library(terra))

# Set variables
proj <- args[1]
species <- gsub("output/maps/inla/", "", args[2])

# Import functions
source("R/plot_map.R")
source("R/make_map.R")
source("R/make_gif.R")

# Import spatial objects
q <- readRDS("data/spacePoly.rds")
qc <- readRDS("data/qc_spacePoly.rds")
rast <- raster::stack("data/rast.gri")
mesh <- readRDS("data/mesh.rds")

# Create directories for species
dir.create(paste0("output/maps/inla/", species))
dir.create(paste0("output/maps/inla/", species,"/qc"))
dir.create(paste0("output/maps/inla/", species,"/region"))
dir.create(paste0("output/maps/inla/", species,"/region_occ"))
 
# List all models from the specific species
models <- list.files(paste0("output/models/inla/", species))

# List years of the models
years <- as.integer(gsub(".rds", "", models))


# Loop on models
for(i in 1:length(models)) {
    
  # Import model
  mod <- readRDS(paste0("output/models/inla/", species, "/", models[i]))

  # Make map for entire sPoly to compute AUC
  map_all <- make_map(type = "mean_all",
                      mesh,
                      mod,
                      rast,
                      sPoly = q,
                      year = years[i])
  
  plot_map(file = paste0("output/maps/inla/", species,"/region/",years[i],".png"),
           map = map_all,
           title = paste0(species,"_",years[i]),
	         region = q)

  png(paste0("output/maps/inla/", species,"/region_occ/",years[i],"_occ.png"), width = 960, height = 960)
  raster::plot(map_all, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = paste0(species,"_",years[i]))
  sp::plot(q, 
           lwd=0.2, 
           add = TRUE)
  dev.off()


  # Make map for Qc only
  map <- terra::crop(terra::rast(map_all), terra::vect(qc))
  map <- terra::mask(map, terra::vect(qc))
  map <- raster::raster(map)

  plot_map(file = paste0("output/maps/inla/", species,"/qc/",years[i],".png"),
           map = map,
           title = paste0(species,"_",years[i]),
           region = qc)

  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack")) {
    map_stack <- raster::stack(map_stack, map)
  } else {
    map_stack <- raster::stack(map)
  }

  # Clean  
  rm(map_all)
  rm(map)
  rm(mod)

  # Save the stack of maps if it's the last year
  if(i == length(models)) {
    raster::writeRaster(map_stack, paste0("output/maps/inla/", species, "/maps"))
  }

  cat(paste0(years[i], " done\n"))
}

make_gif(folder = paste0("output/maps/inla/",species,"/"), region = "region")
make_gif(folder = paste0("output/maps/inla/",species,"/"), region = "qc")

# Clean
rm(map_stack)
