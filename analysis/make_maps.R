# Binarize maps with LO_CI

#---------- SETUP ----------#

args <- commandArgs(trailingOnly=TRUE)

suppressMessages(library(raster))
suppressMessages(library(INLA))
suppressMessages(library(terra))
suppressMessages(library(rangemap))

species <- args[1]
proj <- args[2]

source("R/path.R")

# Objet spatial propre à chaque espèce
study_extent <- readRDS(paste0(path_sp(species)$spat, "/study_extent.rds"))
rast <- raster::raster(paste0(path_sp(species)$spat, "/rast.gri")) # Raster pour aggréger les données par cellule
qc <- readRDS("data/qc_spacePoly.rds")
mesh <- readRDS(paste0(path_sp(species)$spat, "/mesh.rds"))
rast_qc <- terra::crop(terra::mask(terra::rast(rast),
                                   terra::vect(qc)),
                       terra::vect(qc))
rast_qc <- raster::raster(rast_qc)
values(rast_qc)[values(rast_qc) == 1] <- 0

# List all models from the specific species
models <- list.files(path_sp(species)$mod)

# List years of the models
years <- as.integer(gsub(".rds", "", models))


for(year in years) {
    
  #---------- Importer modèle et stack ----------#
  
  model <- readRDS(paste0(path_sp(species)$mod, "/",year,".rds"))
  Stack <- readRDS(paste0(path_sp(species)$stack, "/",year,".rds"))


  #---------- Faire les cartes (moyenne et CI) ----------#

  mapBasis <- inla.mesh.projector(mesh,
                                  dims = dim(rast)[2:1],
                                  xlim = c(xmin(study_extent), 
                                           xmax(study_extent)),
                                  ylim = c(ymin(study_extent), 
                                           ymax(study_extent)),
                                  crs = mesh$crs)

  ## Trouver les arêtes du mesh pour y faire les prédiction
  ID <- inla.stack.index(Stack, tag="pred")$data

  ## Prédictions avec 0.5quant
  mapPred <- inla.mesh.project(mapBasis, 
                               model$summary.fitted.values[["0.5quant"]][ID])

  ## Convertir en raster
  map <- raster::raster(t(mapPred[,ncol(mapPred):1]),
                          xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                          ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                          crs = mesh$crs)

  ## Prédictions avec 0.025 & 0.975 quant
  mapPred025 <- inla.mesh.project(mapBasis, 
                              model$summary.fitted.values[["0.025quant"]][ID])
  mapPred975 <- inla.mesh.project(mapBasis, 
                              model$summary.fitted.values[["0.975quant"]][ID])

  ## Convertir en raster
  map025 <- raster::raster(t(mapPred025[,ncol(mapPred025):1]),
                             xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                             ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                             crs = mesh$crs)
  map975 <- raster::raster(t(mapPred975[,ncol(mapPred975):1]),
                             xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                             ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                             crs = mesh$crs)


  #---------- Crop les maps avec study_extent ----------#

  map <- terra::crop(terra::mask(terra::rast(map), 
                                 terra::vect(study_extent)), 
                     terra::vect(study_extent))
  map025 <- terra::crop(terra::mask(terra::rast(map025), 
                                    terra::vect(study_extent)), 
                        terra::vect(study_extent))
  map975 <- terra::crop(terra::mask(terra::rast(map975), 
                                    terra::vect(study_extent)), 
                        terra::vect(study_extent))

  ## Reconvertir en raster
  map_pocc <- raster::raster(map)
  names(map_pocc) <- paste0(year)
  map025 <- raster::raster(map025)
  names(map025) <- paste0(year)
  map975 <- raster::raster(map975)
  names(map975) <- paste0(year)


  #---------- Binariser les maps ----------#

  # Est-ce que la limite inférieure inclut 0?
  
  map <- map_pocc
  map[map025<= 0.01] <- 0
  map[map>0] <- 1
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
  rangemap <- rangemap_hull(xy,
                            hull_type = "concave", 
                            concave_distance_lim = 350000, 
                            buffer_distance = 10000,
                            final_projection = proj,
                            extent_of_occurrence = TRUE,
                            area_of_occupancy = TRUE)
  rangemap_qc <- terra::intersect(
                   terra::vect(rangemap@species_range),
                   terra::vect(qc))
  rangemap_qc <- as(rangemap_qc, "Spatial")
  map <- raster::crop(
    raster::rasterize(raster::aggregate(rangemap_qc),
                      rast, 
                      mask = TRUE),
    rangemap_qc)
  map <- raster::merge(map, rast_qc)
  names(map) <- paste0(year)

  #---------- Stack et sauvegarder les maps ----------#

  ## Stack des cartes de chaque année ensemble
  if(exists("sdms")) {
    sdms <- raster::stack(sdms, map)
  } else {
    sdms <- raster::stack(map)
  }

  ## Stack des cartes avec pocc de chaque année ensemble
  if(exists("sdms_pocc")) {
    sdms_pocc <- raster::stack(sdms_pocc, map_pocc)
  } else {
    sdms_pocc <- raster::stack(map_pocc)
  }

  ## Stack des cartes de chaque année ensemble
  if(exists("sdms025")) {
    sdms025 <- raster::stack(sdms025, map025)
  } else {
    sdms025 <- raster::stack(map025)
  }
  if(exists("sdms975")) {
    sdms975 <- raster::stack(sdms975, map975)
  } else {
    sdms975 <- raster::stack(map975)
  }

  ## Si c'est la dernière année, sauvegarder le stack
  if(year == 2018) {
    raster::writeRaster(sdms, 
                        paste0(path_sp(species)$maps, "/maps_binary"))
    raster::writeRaster(sdms_pocc, 
                        paste0(path_sp(species)$maps, "maps_pocc"))     
    raster::writeRaster(sdms025, 
                        paste0(path_sp(species)$maps, "/maps_025"))
    raster::writeRaster(sdms975, 
                        paste0(path_sp(species)$maps, "/maps_975"))
  }
}
