
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
year_start <- as.integer(args[2])
year_end <- as.integer(args[3])
window_width <- as.integer(args[4])
num_threads <- as.integer(args[5])
zone <- args[6]
obs_all <- readRDS(paste0(path_sp(species, zone)$spat, "/obs.rds"))

#--- Continue setup ---#
# Import spatial objects to make models
q <- readRDS(paste0(path_sp(species, zone)$spat,"/study_extent.rds"))
qc <- readRDS(paste0("data/",zone,"qc_spacePoly.rds"))
mesh <- readRDS(paste0(path_sp(species, zone)$spat,"/mesh.rds"))
rast <- raster::stack(paste0(path_sp(species, zone)$spat,"/rast.gri"))

source("R/make_spde.R")
source("R/make_stack.R")
source("R/success_trials.R")

# Species name as a folders
dir.create(path_sp(species, zone)$mod)
dir.create(path_sp(species, zone)$log)
dir.create(path_sp(species, zone)$stack)

#--- Make models ---#
# First time window
years <- year_start:(year_start+window_width-1)

# Make models for every time windows
for(j in 0:(year_end-years[length(years)])) {
 ## Mean year of the time window
  year <- mean((years+j))

   ## Filter observations
  obs <- obs_all[obs_all$year_obs %in% (years+j),]
  obs <- convert_in_success_trials(obs, rast)
  obs <- obs[obs$observations > 0,]

  if(sum(obs$presences) > 25) {
    
    ## Make model and pray that it makes sense
    tryCatch(
      expr = {

        spde <- make_spde(mesh)
        Stack <- make_stack(mesh, obs, spde)

        # Step 8 - Building the model
        if(j != 0) {

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
                      control.mode = list(theta = previous_model$mode$theta, restart = TRUE),
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

        saveRDS(model, paste0(path_sp(species, zone)$mod,"/",year,".rds"))
        saveRDS(Stack, paste0(path_sp(species, zone)$stack,"/",year,".rds"))
        previous_model <- model
        rm(model)
        rm(obs)
        rm(spde)
        rm(Stack)
      },
      error=function(cond) {
        cat(paste0("Error ", species, " for year ", year, ": ", cond), 
            sep = "\n\n", 
            file = paste0(path_sp(species, zone)$log, "/log"), 
            append = file.exists(paste0(path_sp(species, zone)$log, "/log")))
      }
    )
  }
}
