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
species <- args[1]
proj <- args[2]

# Import functions
source("R/plot_map.R")
source("R/make_map.R")
source("R/make_gif.R")
source("R/path.R")
source("R/binarize_maps.R")

# Create directories for species
dir.create(path_sp(species)$maps)
dir.create(path_maps(species)$qc)
dir.create(path_maps(species)$region)
dir.create(path_maps(species)$region_pres)
dir.create(path_maps(species)$region_abs)

# Import spatial objects
study_extent <- readRDS(paste0(path_sp(species)$spat, "/study_extent.rds"))
qc <- readRDS("data/qc_spacePoly.rds")
rast <- raster::stack(paste0(path_sp(species)$spat, "/rast.gri"))
mesh <- readRDS(paste0(path_sp(species)$spat, "/mesh.rds"))
obs_all <- readRDS(paste0(path_sp(species)$spat, "/obs.rds"))

# List all models from the specific species
models <- list.files(path_sp(species)$mod)
Stacks <- list.files(path_sp(species)$stack)

# List years of the models
years <- as.integer(gsub(".rds", "", models))

# Loop on models
for(i in 1:length(models)) {
  
  obs <- obs_all[obs_all$year_obs %in% (years[i]-2):(years[i]+2),]

  # Import model
  mod <- readRDS(paste0(path_sp(species)$mod, "/", models[i]))

  Stack <- readRDS(paste0(path_sp(species)$stack, "/", Stacks[i]))

  # Make map for entire sPoly to compute AUC
  map <- make_map(type = "mean",
                  mesh,
                  mod,
                  rast,
                  sPoly = study_extent,
                  year = years[i],
                  Stack)
  map_0025 <- make_map(type = "0.025quant",
			                 mesh,
			                 mod,
			                 rast,
			                 sPoly = study_extent,
			                 year = years[i],
			                 Stack) 
  map_0975 <- make_map(type = "0.975quant",
                       mesh,
                       mod, 
                       rast,
                       sPoly = study_extent,
                       year = years[i],
                       Stack)  
  
  # make binary maps
  ## Find threshold
  threshold <- find_threshold(map, obs, "sensitivity")
  map_binary <- binarize_pred(map, threshold)

  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack")) {
    map_stack_pocc <- raster::stack(map_stack_pocc, map)
  } else {
    map_stack_pocc <- raster::stack(map)
  }
  
  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack_uncert")) {
    map_stack_uncert <- raster::stack(map_stack_uncert, map_0025, map_0975)
  } else {
    map_stack_uncert <- raster::stack(map_0025, map_0975)
  }

  # Create binary map stack
  if(exists("map_stack_binary")) {
    map_stack_binary <- raster::stack(map_stack_binary, map_binary)
  } else {
    map_stack_binary <- raster::stack(map_binary)
  }

  # Clean  
  rm(map_pocc)
  rm(map_0025)
  rm(map_0975)
  rm(map_binary)
  rm(mod)
  rm(Stack)

  # Save the stack of maps if it's the last year
  if(i == length(models)) {
    raster::writeRaster(map_stack_pocc, paste0(path_sp(species)$maps, "/maps_occ"))
    raster::writeRaster(map_stack_uncert, paste0(path_sp(species)$maps, "/maps_uncert"))
    raster::writeRaster(map_stack_binary, paste0(path_sp(species)$maps, "/maps_binary"))
  }

  cat(paste0(years[i], " done\n"))
}

# Clean
rm(map_stack_pocc)
rm(map_stack_uncert)
rm(map_stack_binary)
