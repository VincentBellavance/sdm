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
suppressMessages(library(mapSpecies))

# Set variables
threshold <- seq(args[1], args[2], by = 0.05)
proj <- args[3]
species <- gsub("output/maps/", "", args[4])

# Import functions
source("R/calc_auc.R")
source("R/make_map.R")

# Import spatial objects
q <- readRDS("data/spacePoly.rds")
qc <- readRDS("data/qc_spacePoly.rds")
rast <- raster::stack("data/rast.gri")

# Create output/maps directory if it doesn't exist
if(!dir.exists("output/maps")) {
  dir.create("output/maps")
}

# Create output/auc directory if it doesn't exist
if(!dir.exists("output/auc")) {
  dir.create("output/auc")
}

# Create directories for species
dir.create(paste0("output/maps/", species))
dir.create(paste0("output/auc/", species))
  
# List all models from the specific species
models <- list.files(paste0("output/models/", species))

# List years of the models
years <- as.integer(gsub(".rds", "", models))

# Make a list to store AUCs
auc <- list()

# Loop on models
for(i in 1:length(models)) {
    
  # Import model
  mod <- readRDS(paste0("output/models/", species, "/", models[i]))

  # Make map for entire sPoly to compute AUC
  map_all <- make_map(type = "mean_all",
                      mod,
                      rast,
                      sPoly = q,
                      year = years[i],
                      proj)

  # Compute threshold
  auc[[i]] <- calc_auc(mod, 
                       map = map_all, 
                       threshold)

  # Make map for Qc only
  map <- make_map(type = "mean",
                  mod,
                  rast,
                  sPoly = qc,
                  year = years[i],
                  proj)

  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack")) {
    map_stack <- raster::stack(map_stack, map)
  } else {
    map_stack <- raster::stack(map)
  }

  # Save the stack of maps if it's the last year
  if(i == length(models)) {
    
    # Save treshold
    names(auc) <- years
    saveRDS(auc, paste0("output/auc/", species, "/auc.rds"))

    # Save maps
    raster::writeRaster(map_stack, paste0("output/maps/", species, "/maps"))
  }

  cat(paste0(years[i], " done\n"))
}