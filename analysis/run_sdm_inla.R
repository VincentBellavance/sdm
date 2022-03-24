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
species <- gsub("output/models/inla/", "", args[1])
obs_all <- readRDS(paste0("occurrences/", species, ".rds"))
species <- gsub("_", " ", species)
species <- stringr::str_to_sentence(species)
year_start <- as.integer(args[2]) # 1990
year_end <- as.integer(args[3]) # 2020
window_width <- as.integer(args[4]) # 5
num_threads <- as.integer(args[5])

#--- Continue setup ---#
# Import spatial objects to make models
q <- readRDS("data/spacePoly.rds")
qc <- readRDS("data/qc_spacePoly.rds")
mesh <- readRDS("data/mesh.rds")

# List of species
species_list <- readRDS("data/species.rds")
which_species <- lapply(species_list, "[[", 1) %in% species
species <- species_list[[which(which_species)]]

#--- Launch loop ---#
# Models species
# Species name as a folders
folder <- tolower(gsub(" ", "_", species$accepted))
dir.create(paste0("output/models/inla/", folder))
dir.create(paste0("output/log/inla/", folder))

#--- Make models ---#
# First time window
years <- year_start:(year_start+window_width-1)


# Make models for every time windows
for(j in 0:(year_end-years[length(years)])) {
 
 ## Mean year of the time window
  year <- mean((years+j))

   ## Filter observations
  obs <- obs_all[obs_all$year_obs %in% (years+j),]

  if(nrow(obs[obs$occurrence,]) > 25) {
    
    ## Make model and pray that it makes sense
    tryCatch(
      expr = {

        # Step 2 - Define the stochastic partial differential equation object
        spde <- inla.spde2.pcmatern(mesh=mesh,
                                    alpha=2,
                                    prior.range=c(100, 0.01),
                                    prior.sigma=c(0.5, 0.01))


        # Step 3 - Index matrix
        field <- inla.spde.make.index("field", n.spde=mesh$n)


        # Step 4 - A matrix
        xy <- as.matrix(cbind(obs$lon, obs$lat))
        colnames(xy) <- c("x", "y")
        ## For estimation
        AEst <- inla.spde.make.A(mesh, loc=xy)
        ## For prediction
        APred <- inla.spde.make.A(mesh)


        # Step 5 - Organise the A matrix into a list
        ## For estimation
        AEstlist <- list(AEst)
        ## For prediction
        APredlist <- list(APred)


        # Step 6 - Organise the effects (field variables)
        ## Estimation
        effectEst <- list(i = 1:spde$n.spde)
        ## Prediction
        effectPred <- list(i = 1:spde$n.spde)


        # Step 7 - Build stack
        ## Stack for estimation
        StackEst <- inla.stack(data=list(occurrence = as.data.frame(obs)[,"occurrence"]),
                               A = AEstlist,
                               effects = effectEst,
                               tag="est")
        ## Stack for prediction
        StackPred <- inla.stack(data=list(occurrence = NA),
                                A=APredlist,
                                effects=effectPred,
                                tag="pred")
        ## Organise StackEst and StackPred into a single stack
        Stack <- inla.stack(StackEst,StackPred)


        # Step 8 - Building the model
        model <- inla(occurrence ~ 0 + f(i, model = SPDE),
                      data = inla.stack.data(Stack),
                      family="binomial",
                      control.family =list(link="logit"),
                      control.compute=list(waic=TRUE,
                                           openmp.strategy = "huge"),
                      num.threads = num_threads,
                      control.inla=list(int.strategy = "ccd"),
                      verbose = TRUE)

        saveRDS(model, paste0("output/models/inla/",folder,"/",year,".rds"))
        rm(model)
        rm(obs)
      },
      error=function(cond) {
        cat(paste0("Error ", folder, " for year ", year, ": ", cond), 
            sep = "\n\n", 
            file = paste0("output/log/inla/", folder, "/log"), 
            append = file.exists(paste0("output/log/", folder, "/log")))
      }
    )
  }
}
