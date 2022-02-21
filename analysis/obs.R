#--- Setup ---#
# Import arguments
args = commandArgs(trailingOnly=TRUE)

# Set variables
species <- gsub("occurrences/", "", args[1])
species <- gsub(".rds", "", species)
species <- gsub("_", " ", species)
species <- stringr::str_to_sentence(species)
year_start <- as.integer(args[2]) # 1990
year_end <- as.integer(args[3]) # 2020
window_width <- as.integer(args[4]) # 5
buffer <- as.integer(args[5]) # 21
proj <- args[6] #"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

# Import functions
for(j in list.files("R", full.names = T)) {
  source(j)
}

# List of species
species_list <- readRDS("data/species.rds")
which_species <- lapply(species_list, "[[", 1) %in% species
species <- species_list[[which(which_species)]]

#--- Get observations ---#
obs <- get_obs(species$accepted, c(year_start, year_end))

# Extract coordinates from geom column
coords <- extract_coords(obs$geom)

# Bind coords and obs together
obs <- cbind(coords, obs)

# Keep species or not
keep <- keep_species(obs)

if(!keep) {
  readLines(con = "data/species_vect.txt") |>
    stringr::str_replace(
      pattern = gsub(" ", "_", tolower(species$accepted)),
      replace = ""
    ) |>
      stringr::str_replace(
        pattern = "  ",
        replace = " "
      ) |>
        writeLines(con = "data/species_vect.txt")

} else {

  # Filter data with dates of observations
  obs <- filter_dates(obs, species, buffer)

  # If at least 5 observations by year for more than 5 consecutive years
  #obs_by_year <- table(obs[obs$occurrence, "year_obs"]) > 5
  #test <- rle(as.logical(obs_by_year))

  # Transform TRUE of FALSE for 1 or 0
  obs$occurrence <- ifelse(obs$occurrence == FALSE, 0, 1)

  # Remove geom column
  obs <- obs[,-which(colnames(obs) == "geom")]

  # Make spatialPointsDataFrame for the observations
  obs <- sp::SpatialPointsDataFrame(obs[,c(1:2)], obs[,-c(1:2)], proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  obs <- sp::spTransform(obs, sp::CRS(proj))

  saveRDS(obs, paste0("occurrences/", gsub(" ", "_", tolower(species$accepted)), ".rds"))

}
