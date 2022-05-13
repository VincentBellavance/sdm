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

  Stack <- readRDS(paste0(path_sp(species)$stack, "/", Stack[i]))

  # Make map for entire sPoly to compute AUC
  map_all <- make_map(type = "mean",
                      mesh,
                      mod,
                      rast,
                      sPoly = q,
                      year = years[i],
                      Stack)
  
  map_0025 <- make_map(type = "0.025quant",
			                 mesh,
			                 mod,
			                 rast,
			                 sPoly = q,
			                 year = years[i],
			                 Stack) 
  map_0975 <- make_map(type = "0.975quant",
                       mesh,
                       mod, 
                       rast,
                       sPoly = q,
                       year = years[i],
                       Stack)  

  plot_map(file = paste0(path_maps(species)$region,"/",years[i],".png"),
           map = map_all,
           title = paste0(species,"_",years[i]),
           region = q)

  png(paste0(path_maps(species)$region_pres,"/",years[i],"_pres.png"), width = 1300, height = 1300)
  raster::plot(map_all, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = paste0(species,"_",years[i]))
  sp::plot(obs[obs$occurrence == 1,"occurrence"], 
           lwd=0.2, 
           add = TRUE)
  dev.off()

  png(paste0(path_maps(species)$region_abs,"/",years[i],"_abs.png"), width = 1300, height = 1300)
  raster::plot(map_all, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = paste0(species,"_",years[i]))
  sp::plot(obs[obs$occurrence == 0,"occurrence"], 
           lwd=0.2, 
           add = TRUE)
  dev.off()

  # Make map for Qc only
  map <- terra::crop(terra::rast(map_all), terra::vect(qc))
  map <- terra::mask(map, terra::vect(qc))
  map <- raster::raster(map)

  plot_map(file = paste0(path_maps(species)$qc,"/",years[i],".png"),
           map = map,
           title = paste0(species,"_",years[i]),
           region = qc)

  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack")) {
    map_stack <- raster::stack(map_stack, map)
  } else {
    map_stack <- raster::stack(map)
  }
  
  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("map_stack_uncert")) {
    map_stack_uncert <- raster::stack(map_stack_uncert, map_0025, map_0975)
  } else {
    map_stack_uncert <- raster::stack(map_0025, map_0975)
  }

  # Clean  
  rm(map_all)
  rm(map)
  rm(map_0025)
  rm(map_0975)
  rm(mod)

  # Save the stack of maps if it's the last year
  if(i == length(models)) {
    raster::writeRaster(map_stack, paste0(path_sp(species)$maps, "/maps"))
    raster::writeRaster(map_stack_uncert, paste0(path_sp(species)$maps, "/maps_uncert"))
  }

  cat(paste0(years[i], " done\n"))
}

make_gif(paste0(path_sp(species)$maps, "/"), "region")
make_gif(paste0(path_sp(species)$maps, "/"), "qc")
make_gif(paste0(path_sp(species)$maps, "/"), "region_pres")
make_gif(paste0(path_sp(species)$maps, "/"), "region_abs")

# Clean
rm(map_stack)
