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
suppressMessages(library(mapSpecies))

# Set variables
species <- gsub("output/models/inla/", "", args[1])
obs <- readRDS(paste0("occurrences/", species, ".rds"))
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
rast <- raster::stack("data/rast.gri")
explana <- readRDS("data/explana.rds")

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

# Setup priors
#prange <- 150
#psigma <- 5

# Make models for every time windows
for(j in 0:(year_end-years[length(years)])) {
 
 ## Mean year of the time window
  year <- mean((years+j))

   ## Filter observations
  obs_filtered <- obs[obs$year_obs %in% (years+j),]

  if(nrow(obs_filtered[obs_filtered$occurrence,]) > 25) {
    
    ## Make model and pray that it makes sense
    tryCatch(
      expr = {
        mod <- mapSpecies::uniSpace(occurrence ~ none,
                                    sPoints = obs_filtered,
                                    explanaMesh = explana, 
                                    family = "binomial",
                                    link = "logit", 
                                    control.compute = list(waic = TRUE,
							   openmp.strategy="huge"),
                                    num.threads = num_threads,
                                    control.inla=list(int.strategy="ccd"), # with ccd or eb, seems less costly in RAM...
                                    prior.range = c(100, 0.01),
                                    prior.sigma = c(0.1, 0.01),
                                    verbose = TRUE)
        
#       prange <- unlist(summary(mod)$spaceCoefficients["Range", c("mean"), drop = T], use.names=F)
#	psigma <- unlist(summary(mod)$spaceCoefficients["Stdev", c("mean"), drop = T], use.names=F)

        saveRDS(mod, paste0("output/models/inla/",folder,"/",year,".rds"))
        rm(mod)
	rm(obs_filtered)
      },
      error=function(cond) {
        cat(paste0("Error ", folder, " for year ", year, ": ", cond), 
            sep = "\n\n", 
            file = paste0("output/log/inla/", folder, "/log"), 
            append = file.exists(paste0("output/log/inla/", folder, "/log")))
      }
    )
  }
}
