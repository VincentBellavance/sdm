
# Make Species Distribution Models that will be used in the BDI calculation

#--------------------------------------------------------------
# Steps:
# 1. Setup environment
# 2. Get observations
# 3. For every time windows:
#   3.1: Aggregate observations
#   3.2: Make model and save it as RDS object
#--------------------------------------------------------------

#--- Setup ---#
# Import arguments
args = commandArgs(trailingOnly=TRUE)

# Import mapSpecies
suppressMessages(library(INLA))
suppressMessages(library(raster))

# Set variables
species <- args[1]
source("R/path.R")
window_width <- as.integer(args[2])
half_wind <- (window_width-1)/2
num_threads <- as.integer(args[3])
zone <- args[4]
output_dir <- args[5]

if(length(list.files(path_sp(species, output_dir, zone = zone)$spat)) != 0) {

  obs_all <- readRDS(paste0(path_sp(species, output_dir, zone = zone)$spat, "/obs.rds"))

  #--- Continue setup ---#
  # Import spatial objects to make models
  q <- readRDS(paste0(path_sp(species, output_dir, zone = zone)$spat,"/study_extent.rds"))
  qc <- readRDS(paste0("data/",zone,"/qc_spacePoly.rds"))
  mesh <- readRDS(paste0(path_sp(species, output_dir, zone = zone)$spat,"/mesh.rds"))
  rast <- raster::stack(paste0(path_sp(species, output_dir, zone = zone)$spat,"/rast.gri"))

  source("R/make_spde.R")
  source("R/make_stack.R")
  source("R/success_trials.R")

  # Species name as a folders
  dir.create(path_sp(species, output_dir, zone = zone)$mod)
  dir.create(path_sp(species, output_dir, zone = zone)$log)
  dir.create(path_sp(species, output_dir, zone = zone)$stack)

  #--- Make models ---#
  # First time window
  years <- attributes(obs_all)$years_mod

  # Make models for every time windows
  for(j in years) {
    
    ## Filter observations
    obs <- obs_all[obs_all$year_obs %in% (j-half_wind):(j+half_wind),]
    obs <- convert_in_success_trials(obs, rast)
    obs <- obs[obs$observations > 0,]

    ## Make model and pray that it makes sense
    tryCatch(
      expr = {

        spde <- make_spde(mesh)
        Stack <- make_stack(mesh, obs, spde)

        # Step 8 - Building the model
        if(exists("previous_model")) {

          model <- inla(presences ~ 0 + f(field, model = spde),
                      data = inla.stack.data(Stack),
                      family="binomial",
                      Ntrials=observations,
                      control.family =list(link="logit"),
                      control.compute=list(waic=TRUE,
                                           openmp.strategy = "huge",
  				                                 config = TRUE),
                      control.predictor=list(A=inla.stack.A(Stack),
                                             compute=TRUE, 
                                             link = 1),
                      control.mode = list(theta = previous_model$mode$theta,
                                          restart = TRUE),
                      control.inla=list(int.strategy = "ccd"),
                      verbose = TRUE,
                      debug = TRUE)
        } else {
          model <- inla(presences ~ 0 + f(field, model = spde),
                      data = inla.stack.data(Stack),
                      family="binomial",
                      Ntrials=observations,
                      control.family =list(link="logit"),
                      control.compute=list(waic=TRUE,
                                           openmp.strategy = "huge",
  				                                 config = TRUE),
                      control.predictor=list(A=inla.stack.A(Stack),
                                             compute=TRUE,
                                             link = 1),
                      control.inla=list(int.strategy = "ccd"),
                      verbose = TRUE,
                      debug = TRUE)
        }

        saveRDS(model, paste0(path_sp(species, output_dir, zone = zone)$mod,"/",j,".rds"))
        saveRDS(Stack, paste0(path_sp(species, output_dir, zone = zone)$stack,"/",j,".rds"))
        previous_model <- model
        rm(model)
        rm(obs)
        rm(spde)
        rm(Stack)
      },
      error=function(cond) {
        cat(paste0("Error ", species, " for year ", j, ": ", cond), 
            sep = "\n\n", 
            file = paste0(path_sp(species, output_dir, zone = zone)$log, "/log"), 
            append = file.exists(paste0(path_sp(species, output_dir, zone = zone)$log, "/log")))
      }
    )
  }

}
