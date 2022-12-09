# This script create, for every species,: 
#   - study extent (a polygon around the presence of a species)
#   - mesh
#   - raster
#   - observations filtered with the study extent

#--- Setup ---#
# Import arguments
args = commandArgs(trailingOnly=TRUE)

# Import mapSpecies
suppressMessages(library(raster))

# Set variables
species <- args[1]
zone <- args[2]
obs_dir <- args[5]
output_dir <- args[6]
year_start <- as.integer(args[7])
year_end <- as.integer(args[8])
time_window <- as.integer(args[9])
source("R/path.R")
source("R/mask_keep_partial.R")
q <- readRDS(paste0("data/",zone,"/spacePoly.rds"))
obs <- readRDS(path_sp(species, obs_dir = obs_dir)$obs)
dist_buffer <- as.integer(args[3])

# Make directory to store new spatial object
if(!exists(path_sp(species, output_dir, zone = zone)$spat)) {
  dir.create(path_sp(species, output_dir, zone = zone)$spat)
}

# Crop observations with new polygon
obs <- terra::intersect(terra::vect(obs), terra::vect(q))
obs <- as(obs, "Spatial")

# Filter obs with obs/models (<30, remove the first year of the time window)
filt_years <- sapply((1:length((year_start+2):(year_end-2)))-1, function(x) {
  years <- year_start:(year_start+time_window-1)+x
  pres <- obs[obs$year_obs %in% years & obs$occurrence == 1, ]
  return(nrow(pres) >= 30)
})

pres_thresh <- FALSE
j = 1
while(pres_thresh == FALSE) {
  pres_thresh <- filt_years[j]
  if(!pres_thresh) {
    obs <- obs[obs$year_obs != ((year_start):(year_end))[j],]
  }
  j <- j+1
  if(j > length(filt_years)) break
}

#TODO: Pas certain du fonctionnement ici, revoir, quoi que pas urgent puisqu'il ne devrait pas y avoir beaucoup d'espèce qui n'ont pas 30 observations en 2018....
pres_thresh <- FALSE
j = length((year_start+2):(year_end-2))
while(pres_thresh == FALSE) {
  pres_thresh <- filt_years[j]
  if(!pres_thresh) {
    obs <- obs[obs$year_obs != ((year_start):(year_end))[j],]
  }
  j <- j-1
  if(j < 0) break
}

if(nrow(obs) > 0 & length(unique(obs$year_obs)) >= 5) {

  # Filter for presence only to define the polygon
  obs_pres <- obs[obs$occurrence == 1,]
  
  # Make the study extent around presence of the species
  xy <- as.data.frame(obs_pres@coords)
  coords.t <- chull(xy[,1], xy[,2])
  xy.bord <- xy[coords.t,]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  study_extent <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(as.matrix(xy.bord))), 1)))
  study_extent <- sf::st_as_sf(study_extent)
  sf::st_crs(study_extent) <- raster::crs(obs_pres@proj4string)
  
  if (!is.null(dist_buffer)) {
    study_extent <-  sf::st_buffer(study_extent, dist =  dist_buffer)
  }
  
  # Filter observations with study extent
  study_extent <- terra::intersect(terra::vect(study_extent), terra::vect(q))
  obs <- terra::intersect(terra::vect(obs), study_extent)
  obs <- as(obs, "Spatial")
  study_extent <- as(study_extent, "Spatial")
  
  # Make mesh
  pedge <- as.numeric(args[4])
  edge <- min(c(diff(raster::bbox(study_extent)[1,])*pedge,diff(raster::bbox(study_extent)[2,])*pedge))
  mesh <- INLA::inla.mesh.2d(boundary = study_extent,
                             max.edge = c(edge, edge*5), 
                             min.angle = 21,
                             offset = c(edge, edge*2), 
                             cutoff = edge/2, 
                             crs = raster::crs(study_extent))
  
  # Crop raster with study extent
  rast <- raster::raster("data/rast.gri")
  rast <- mask_keep_partial(rast, study_extent)
  
  # Save all four objects
  raster::writeRaster(rast, 
                      paste0(path_sp(species,output_dir,zone=zone)$spat,"/rast"))
  saveRDS(mesh, paste0(path_sp(species,output_dir,zone=zone)$spat,"/mesh.rds"))
  # Put years to model as attributes of obs
  range_years <- range(obs$year_obs)
  attr(obs, "years_mod") <- (range_years[1]+2):(range_years[2]-2)
  saveRDS(obs, paste0(path_sp(species,output_dir,zone=zone)$spat,"/obs.rds"))
  saveRDS(study_extent, 
          paste0(path_sp(species,output_dir,zone=zone)$spat,"/study_extent.rds"))

}
