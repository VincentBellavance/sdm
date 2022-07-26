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
suppressMessages(library(dplyr))

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
  thresh_spec_sens <- find_threshold(map, obs, "spec_sens")
  map_spec_sens <- binarize_pred(map, thresh_spec_sens)

  # make binary maps
  ## Find threshold
  thresh_sensitivity <- find_threshold(map, obs, "sensitivity")
  map_sensitivity <- binarize_pred(map, thresh_sensitivity)

  # make binary maps
  ## Find threshold
  thresh_kappa <- find_threshold(map, obs, "kappa")
  map_kappa <- binarize_pred(map, thresh_kappa)

  # make binary maps
  ## Find threshold
  thresh_prevalence <- find_threshold(map, obs, "prevalence")
  map_prevalence <- binarize_pred(map, thresh_prevalence)

  # make binary maps
  ## Find threshold
  thresh_equal_sens_spec <- find_threshold(map, obs, "equal_sens_spec")
  map_equal_sens_spec <- binarize_pred(map, thresh_equal_sens_spec)

  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack_pocc")) {
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
  if(exists("map_stack_spec_sens")) {
    map_stack_spec_sens <- raster::stack(map_stack_spec_sens, map_spec_sens)
  } else {
    map_stack_spec_sens <- raster::stack(map_spec_sens)
  }

  # Create binary map stack
  if(exists("map_stack_sensitivity")) {
    map_stack_sensitivity <- raster::stack(map_stack_sensitivity, map_sensitivity)
  } else {
    map_stack_sensitivity <- raster::stack(map_sensitivity)
  }

  # Create binary map stack
  if(exists("map_stack_kappa")) {
    map_stack_kappa <- raster::stack(map_stack_kappa, map_kappa)
  } else {
    map_stack_kappa <- raster::stack(map_kappa)
  }

  # Create binary map stack
  if(exists("map_stack_prevalence")) {
    map_stack_prevalence <- raster::stack(map_stack_prevalence, map_prevalence)
  } else {
    map_stack_prevalence <- raster::stack(map_prevalence)
  }

  # Create binary map stack
  if(exists("map_stack_equal_sens_spec")) {
    map_stack_equal_sens_spec <- raster::stack(map_stack_equal_sens_spec, map_equal_sens_spec)
  } else {
    map_stack_equal_sens_spec <- raster::stack(map_equal_sens_spec)
  }

  # Clean  
  rm(map_pocc)
  rm(map_0025)
  rm(map_0975)
  rm(map_spec_sens)
  rm(map_sensitivity)
  rm(map_kappa)
  rm(map_prevalence)
  rm(map_equal_sens_spec)
  rm(mod)
  rm(Stack)

  # Save the stack of maps if it's the last year
  if(i == length(models)) {
    raster::writeRaster(map_stack_pocc, paste0(path_sp(species)$maps, "/maps_pocc"))
    raster::writeRaster(map_stack_uncert, paste0(path_sp(species)$maps, "/maps_uncert"))
    raster::writeRaster(map_stack_spec_sens, paste0(path_sp(species)$maps, "/maps_spec_sens"))
    raster::writeRaster(map_stack_sensitivity, paste0(path_sp(species)$maps, "/maps_sensitivity"))  
    raster::writeRaster(map_stack_kappa, paste0(path_sp(species)$maps, "/maps_kappa"))
    raster::writeRaster(map_stack_prevalence, paste0(path_sp(species)$maps, "/maps_prevalence"))
    raster::writeRaster(map_stack_equal_sens_spec, paste0(path_sp(species)$maps, "/maps_equal_sens_spec"))
  }

  cat(paste0(years[i], " done\n"))
}

# Clean
rm(map_stack_pocc)
rm(map_stack_uncert)
rm(map_stack_spec_sens)
rm(map_stack_sensitivity)
rm(map_stack_kappa)
rm(map_stack_prevalence)
rm(map_stack_equal_sens_spec)
