# Binarize maps with LO_CI

#---------- SETUP ----------#

args <- commandArgs(trailingOnly=TRUE)

suppressMessages(library(raster))
suppressMessages(library(INLA))
suppressMessages(library(terra))

species <- args[1]
source("R/path.R")
zone <- args[2]
dir.create(path_sp(species, zone)$maps)

# Objet spatial propre à chaque espèce
study_extent <- readRDS(paste0(path_sp(species, zone)$spat, "/study_extent.rds"))
rast <- raster::raster(paste0(path_sp(species, zone)$spat, "/rast.gri")) # Raster pour aggréger les données par cellule
qc <- readRDS(paste0("data/",zone,"/qc_spacePoly.rds")
mesh <- readRDS(paste0(path_sp(species, zone)$spat, "/mesh.rds"))
org <- raster::raster("data/rast.gri")

# List all models from the specific species
models <- list.files(path_sp(species, zone)$mod)

# List years of the models
years <- as.integer(gsub(".rds", "", models))

if(length(years) == 27) {
  for(year in years) {

    #---------- Importer modèle et stack ----------#

    model <- readRDS(paste0(path_sp(species, zone)$mod, "/",year,".rds"))
    Stack <- readRDS(paste0(path_sp(species, zone)$stack, "/",year,".rds"))


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


    #---------- Project Raster pour être sur le même origin ----------#

    map_pocc <- raster::projectRaster(map, org)
    map025 <- raster::projectRaster(map025, org)
    map975 <- raster::projectRaster(map975, org)

    #---------- Stack et sauvegarder les maps ----------#

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
      raster::writeRaster(sdms_pocc, 
                          paste0(path_sp(species,zone)$maps, "/maps_pocc"))     
      raster::writeRaster(sdms025, 
                          paste0(path_sp(species,zone)$maps, "/maps_025"))
      raster::writeRaster(sdms975, 
                          paste0(path_sp(species,zone)$maps, "/maps_975"))
    }
  }
}
