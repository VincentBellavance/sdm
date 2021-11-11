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
year_start <- as.integer(args[1]) # 1990
year_end <- as.integer(args[2]) # 2020
window_width <- as.integer(args[3]) # 5
buffer <- as.integer(args[4]) # 21
proj <- args[5] #"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
num_threads <- as.integer(args[6])
cat("variables imported\n")

# Import functions
for(j in list.files("R", full.names = T)) {
  source(j)
}
cat("functions imported\n")

# Import spatial objects to make models
q <- readRDS("data/spacePoly.rds")
qc <- readRDS("data/qc_spacePoly.rds")
mesh <- readRDS("data/mesh.rds")
rast <- raster::stack("data/rast.gri")
explana <- readRDS("data/explana.rds")
cat("spatObject imported\n")

# Create output/models directory if it doesn't exist
if(!dir.exists("output/models")) {
  dir.create("output/models")
}

# List of species
species <- readRDS("data/species.rds")


#--- Launch loop ---#
# Models every species in a loop
for(i in 1:length(species)) {
  cat(paste0("species ",i,"\n"))
  # Species name as a folders
  folder <- tolower(gsub(" ", "_", species[[i]]$accepted))
  dir.create(paste0("output/models/", folder))

  #--- Get observations ---#
  obs <- get_obs(species[[i]]$accepted)
  cat("observations imported\n")

  # Filter data with dates of observations
  obs <- obs[obs$year_obs >= year_start & obs$year_obs <= year_end, ]
  obs <- filter_dates(obs, species[[i]], buffer)

  # Transform TRUE of FALSE for 1 or 0
  obs$occurrence <- ifelse(obs$occurrence == FALSE, 0, 1)

  # Extract coordinates from geom column
  coords <- extract_coords(obs$geom)

  # Bind coords and obs together
  obs <- cbind(coords, obs)

  # Remove geom column
  obs <- obs[,-which(colnames(obs) == "geom")]

  # Make spatialPointsDataFrame for the observations
  obs <- sp::SpatialPointsDataFrame(coords, obs[,-c(1:2)], proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  obs <- sp::spTransform(obs, sp::CRS(proj))
  cat("observations cleaned\n")

  #--- Make models ---#
  # First time window
  years <- year_start:(year_start+window_width-1)

  # Make models for every time windows
  for(j in 0:(year_end-years[length(years)])) {

    ## Mean year of the time window
    year <- mean((years+j))

    ## Aggregate observations
    obs_aggr <- aggregate_obs(obs, years, j, rast)
    cat(paste0("observations aggregated for year ", year, "\n"))
    if(nrow(obs_aggr[obs_aggr$occurrence,]) > 25) {
      ## Make model and pray that it makes sense
      mod <- mapSpecies::uniSpace(occurrence ~ none,
                                  sPoints = obs_aggr,
                                  explanaMesh = explana, 
                                  family = "binomial",
                                  link = "logit", 
                                  control.compute = list(waic = TRUE),
                                  num.threads = num_threads,
                                  control.inla=list(int.strategy="ccd"), # with ccd or eb, seems less costly in RAM...
                                  prior.range = c(150, 0.05),
                                  prior.sigma = c(50, 0.05))

      tmp <- mod["summary.fitted.values"]
      attributes(tmp) <- attributes(mod)[c("meshSpace", "Stack", "sPoints")]
      names(tmp) <- "summary.fitted.values"

      mod <- tmp

      saveRDS(mod, paste0("output/models/",folder,"/",year,".rds"))
    }
  }
}
