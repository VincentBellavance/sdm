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
suppressMessages(library(mgcv))

# Set variables
species <- gsub("output/models/", "", args[1])
obs <- readRDS(paste0("occurrences/", species, ".rds"))
species <- gsub("_", " ", species)
species <- stringr::str_to_sentence(species)
year_start <- as.integer(args[2]) # 1990
year_end <- as.integer(args[3]) # 2020
window_width <- as.integer(args[4]) # 5

# List of species
species_list <- readRDS("data/species.rds")
which_species <- lapply(species_list, "[[", 1) %in% species
species <- species_list[[which(which_species)]]

#--- Launch loop ---#
# Models species
# Species name as a folders
folder <- tolower(gsub(" ", "_", species$accepted))
dir.create(paste0("output/models/", folder))
dir.create(paste0("output/log/", folder))

#--- Make models ---#
# First time window
years <- year_start:(year_start+window_width-1)

# Make models for every time windows
for(j in 0:(year_end-years[length(years)])) {

  ## Mean year of the time window
  year <- mean((years+j))

  ## Aggregate observations
  obs_filtered <- obs[obs$year_obs %in% (years+j),]

  if(nrow(obs_filtered[obs_filtered$occurrence,]) > 30) {
    
    ## Make model and pray that it makes sense
    tryCatch(
      expr = {
        mod <- mgcv::gam(occurrence ~ s(lon, lat, k = 2), data = obs_filtered, family = "binomial", method = "REML")

        saveRDS(mod, paste0("output/models/",folder,"/",year,".rds"))
      },
      error=function(cond) {
        cat(paste0("Error ", folder, " for year ", year, ": ", cond), 
            sep = "\n\n", 
            file = paste0("output/log/", folder, "/log"), 
            append = file.exists(paste0("output/log/", folder, "/log")))
      }
    )
  }
}
