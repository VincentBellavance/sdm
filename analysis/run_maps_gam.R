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
suppressMessages(library(mgcv))
suppressMessages(library(raster))

# Set variables
#threshold <- seq(args[1], args[2], by = 0.05)
#proj <- args[3]
species <- gsub("output/maps/gam/", "", args[1])
obs_all <- readRDS(paste0("occurrences/",species,".rds"))

# Import functions
#source("R/calc_auc.R")
source("R/make_map_gam.R")
source("R/plot_map.R")
source("R/make_gif.R")

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
dir.create(paste0("output/maps/gam/", species))
dir.create(paste0("output/maps/gam/", species,"/qc"))
dir.create(paste0("output/maps/gam/", species,"/region"))
dir.create(paste0("output/maps/gam/", species,"/region_occ/"))
dir.create(paste0("output/maps/gam/", species,"/region_abs/"))
dir.create(paste0("output/maps/gam/", species,"/qc/png/"))
dir.create(paste0("output/maps/gam/", species,"/region/png/"))
#dir.create(paste0("output/auc/", species))
  
# List all models from the specific species
models <- list.files(paste0("output/models/gam/", species))

# List years of the models
years <- as.integer(gsub(".rds", "", models))

# Make a list to store AUCs
#auc <- list()

# Loop on models
for(i in 1:length(models)) {
    
  # Import model
  mod <- readRDS(paste0("output/models/gam/", species, "/", models[i]))
  obs <- obs_all[obs_all$year_obs %in% (years[i]-2):(years[i]+2),]

  # Make map for entire sPoly to compute AUC
  map_all <- make_map_gam(mod, rast)
  names(map_all) <- paste0("region_", years[i])
  raster::crs(map_all) <- raster::crs(rast)
  map_all <- terra::crop(terra::rast(map_all), terra::vect(q)) |>
               {\(.) terra::mask(., terra::vect(q))}() |>
                 {\(.) raster::raster(.)}()
  # Compute threshold
  #auc[[i]] <- calc_auc(mod, 
  #                     map = map_all, 
  #                     threshold)

  # Make map for Qc only
  map <- terra::crop(terra::rast(map_all), terra::vect(qc)) |>
    {\(.) terra::mask(., terra::vect(qc))}() |>
      {\(.) raster::raster(.)}()
  names(map) <- paste0("qc_", years[i])
  raster::crs(map) <- raster::crs(rast)
  
  plot_map(paste0("output/maps/gam/", species, "/qc/png/",years[i],".png"),
           map,
           paste0(species,"_", years[i]),
           qc)
  plot_map(paste0("output/maps/gam/", species, "/region/png/",years[i],".png"),
           map_all,
           paste0(species,"_", years[i]),
           q)
  png(paste0("output/maps/gam/", species,"/region_occ/",years[i],"_pres.png"), width = 1300, height = 1300)
  raster::plot(map_all, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = paste0(species,"_",years[i]))
  sp::plot(obs[obs$occurrence == 1,"occurrence"], 
           lwd=0.2, 
           add = TRUE)
  dev.off()

  png(paste0("output/maps/gam/", species,"/region_abs/",years[i],"_abs.png"), width = 1300, height = 1300)
  raster::plot(map_all, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = paste0(species,"_",years[i]))
  sp::plot(obs[obs$occurrence == 0,"occurrence"], 
           lwd=0.2, 
           add = TRUE)
  dev.off()
  

  # Create stack if it doesn't exist, else stack the map to the existing one
  if(exists("region_stack")) {
    region_stack <- raster::stack(region_stack, map_all)
  } else {
    region_stack <- raster::stack(map_all)
  }

  if(exists("qc_stack")) {
    qc_stack <- raster::stack(qc_stack, map)
  } else {
    qc_stack <- raster::stack(map)
  }

  # Save the stack of maps if it's the last year
  if(i == length(models)) {
    
    # Save treshold
  #  names(auc) <- years
  #  saveRDS(auc, paste0("output/auc/", species, "/auc.rds"))

    # Save maps
    raster::writeRaster(region_stack, paste0("output/maps/gam/", species, "/qc/maps"))
    raster::writeRaster(qc_stack, paste0("output/maps/gam/", species, "/region/maps"))
  }

  cat(paste0(years[i], " done\n"))
}

make_gif(paste0("output/maps/gam/", species, "/qc"))
make_gif(paste0("output/maps/gam/", species, "/region"))
make_gif(paste0("output/maps/gam/", species, "/region_occ"))
make_gif(paste0("output/maps/gam/", species, "/region_abs"))
