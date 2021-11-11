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

# List directories in '/output/models' that correspond to species names
species <- list.dirs("output/models", recursive = FALSE, full.names = FALSE)


#--- Launch loop ---#
for(i in species) {

  # Create directories for species 'i'
  dir.create(paste0("output/maps", i))
  dir.create(paste0("output/auc", i))
  
  # List all models from the specific species 'i'
  models <- list.files(paste0("output/models/", i))

  # List years of the models
  years <- as.integer(gsub(".rds", "", models))

  # Loop on models
  for(j in 1:length(models)) {
    
    # Import model
    mod <- readRDS("output/models/", i, "/", models[j])

    # Make map for entire sPoly to compute AUC
    map_all <- make_map(type = "mean_all",
                        mod,
                        rast,
                        sPoly = q,
                        year = years[j],
                        proj)

    # Compute threshold
    auc <- calc_auc(mod, 
                    map = map_all, 
                    threshold)

    # Save treshold
    saveRDS(auc, paste0("output/auc/", i, "/", years[j],".rds"))
  
    # Make map for Qc only
    map <- make_map(type = "mean",
                    mod,
                    rast,
                    sPoly = qc,
                    year = years[j],
                    proj)

    # Create stack if it doesn't exist, else stack the map to the existing one
    if(exists("map_stack")) {
      raster::stack(map_stack, map)
    } else {
      raster::stack(map)
    }

    # Save the stack of maps if it's the last year
    if(j == length(models)) {
      raster::writeRaster(map_stack, paste0("output/maps/", i, "/maps"))
    }

    # Clean before the next iteration
    rm(mod, map_all, auc, map)

  }

}