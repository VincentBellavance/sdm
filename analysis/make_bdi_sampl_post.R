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
iter <- as.integer(args[1])

# Import functions
# Import rbdi library
source("../rbdi/R/calc_bdi.R")
source("../rbdi/R/calc_lgr.R")
source("../rbdi/R/utils.R")
source("R/path.R")

qc <- readRDS("data/qc_spacePoly.rds")
species <- as.character(read.table("data/species_vect.txt", sep = " ")[1,])

sdms <- lapply(species, function(i) {

  # List all models from the specific species
  models <- list.files(path_sp(i)$mod)

  if(length(models) == 27) {
  
    Stacks <- list.files(path_sp(i)$stack)

    # Create directories for species
    dir.create(path_sp(i)$maps)

    # Import spatial objects
    study_extent <- readRDS(paste0(path_sp(i)$spat, "/study_extent.rds"))
    rast <- raster::stack(paste0(path_sp(i)$spat, "/rast.gri"))
    mesh <- readRDS(paste0(path_sp(i)$spat, "/mesh.rds"))
    obs_all <- readRDS(paste0(path_sp(i)$spat, "/obs.rds"))

    # List years of the models
    years <- as.integer(gsub(".rds", "", models))

    maps_all_years <- lapply(1:length(years), function(x) {

      # Import model
      mod <- readRDS(paste0(path_sp(i)$mod, "/", models[x]))
      Stack <- readRDS(paste0(path_sp(i)$stack, "/", Stacks[x]))

      # Make map for entire sPoly to compute AUC
      # mapBasis
      mapBasis <- inla.mesh.projector(mesh,
                                      dims = dim(rast)[2:1],
                                      xlim = c(xmin(study_extent), xmax(study_extent)),
                                      ylim = c(ymin(study_extent), ymax(study_extent)),
                                      crs = mesh$crs)

      ### Find the mesh edges on which predictions should be made
      ID <- inla.stack.index(Stack, tag="pred")$data

      pred <- inla.posterior.sample(1, 
                                    result = mod, 
                                    selection = list(APredictor = 0, Predictor = 0),
                                    intern = TRUE,
                                    use.improved.mean = FALSE)

      fitted <- exp(pred[[1]]$latent[ID,])/(1 + exp(pred[[1]]$latent[ID,]))
      mapPred <- inla.mesh.project(mapBasis, fitted)
      mapRaster <- raster::raster(t(mapPred[,ncol(mapPred):1]),
                                    xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                                    ymn = min(mapBasis$y), ymx = max(mapBasis$y),
                                    crs = mesh$crs)
      mapRaster <- terra::mask(terra::rast(mapRaster), terra::vect(qc))
      map <- raster::raster(mapRaster)
      names(map) <- paste0(years[x])
      cat(paste0(years[x], " done\n"))
      return(map)  
    })

    maps_all_years <- do.call(raster::stack, maps_all_years)

  }

})

filter_sp <- !sapply(sdms, is.null) & sapply(sdms, function(x) length(names(x)) == 27)
species <- species[filter_sp]
sdms <- sdms[filter_sp]
names(sdms) <- species

# Calculate BDI
bdi <- calc_bdi(sdms, first_year = 1992)

# Save bdi object
write.table(bdi,
            paste0("~/scratch/output/bdi_improved_mean/bdi",iter,".txt"),
            sep = ",",
            dec = ".")
