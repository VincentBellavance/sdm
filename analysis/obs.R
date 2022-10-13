#--- Setup ---#
# Import arguments
args = commandArgs(trailingOnly=TRUE)

# Set variables
species <- as.vector(read.table("data/species_vect.txt", sep = " "))
species <- gsub("_", " ", species)
species <- stringr::str_to_sentence(species)
year_start <- as.integer(args[1]) # 1990
year_end <- as.integer(args[2]) # 2020
buffer <- as.integer(args[3]) # 21
proj <- args[4] #"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

# Import functions
source("R/get_obs.R")
source("R/filter_dates.R")
source("R/extract_coords.R")

# Create occurrence folder
dir.create("data/occurrences")

# Connection to DB
con <- atlasBE::conn(user=Sys.getenv("user"), pwd=Sys.getenv("pwd"), host=Sys.getenv("host"), dbname=Sys.getenv("dbname"))


for (i in species) {
  
  # Get observations
  obs <- get_obs(con, i, c(year_start, year_end))

  # Extract coordinates from geom column
  coords <- extract_coords(obs$geom)

  # Bind coords and obs together
  obs <- cbind(coords, obs)

  # Filter data with dates of observations
  obs <- filter_dates(con, obs, i, buffer)

  # Transform TRUE of FALSE for 1 or 0
  obs$occurrence <- ifelse(obs$occurrence == FALSE, 0, 1)

  # Remove geom column
  obs <- obs[,-which(colnames(obs) == "geom")]

  # Make spatialPointsDataFrame for the observations
  obs <- sp::SpatialPointsDataFrame(obs[,c(1:2)], 
                                    obs[,-c(1:2)], 
                                    proj4string = sp::CRS("EPSG:4326"))
  obs <- sp::spTransform(obs, sp::CRS(proj))

  # Save spdf
  lower_sp <- gsub(" ", "_", tolower(species$accepted))
  saveRDS(obs, paste0("data/occurrences/", lower_sp, ".rds"))

  rm(obs, coords, lower_sp)
  
}
