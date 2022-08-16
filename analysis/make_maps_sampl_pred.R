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
source("R/make_maps_bs.R")
source("R/path.R")

# Create directories for species
dir.create(path_sp(species)$maps)

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

maps <- lapply(1:length(years), function(x) {

  # Import model
  mod <- readRDS(paste0(path_sp(species)$mod, "/", models[x]))
  Stack <- readRDS(paste0(path_sp(species)$stack, "/", Stacks[x]))
  
  # Make map for entire sPoly to compute AUC
  maps <- make_maps_bs(mod, 
                       mesh, 
                       rast,
                       qc,
                       Stack,
                       years[x])

  cat(paste0(years[x], " done\n"))
  return(maps)

})

maps <- do.call(raster::stack, maps)

raster::writeRaster(maps, paste0(path_sp(species)$maps, "/maps_bs"))
