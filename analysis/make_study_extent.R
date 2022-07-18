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
source("R/path.R")
q <- readRDS("data/spacePoly.rds")
obs <- readRDS(path_sp(species)$occ)
dist_buffer <- as.integer(args[2])

# Make directory to store new spatial object
if(!exists(path_sp(species)$spat)) {
  dir.create(path_sp(species)$spat)
}

# Crop observations with new polygon
obs <- terra::intersect(terra::vect(obs), terra::vect(q))
obs <- as(obs, "Spatial")
raster::crs(obs) <- raster::crs(q)

# Filter for presence only to define the polygon
obs_pres <- obs[obs$occurrence == 1,]

# Make the study extent around presence of the species
xy <- as.data.frame(obs_pres@coords)
coords.t <- chull(xy[,1], xy[,2])
xy.bord <- xy[coords.t,]
xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
study_extent <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(as.matrix(xy.bord))), 1)))
study_extent <- sf::st_as_sf(study_extent, crs = raster::crs(obs_pres@proj4string))

if (!is.null(dist_buffer)) {
  study_extent <-  sf::st_buffer(study_extent, dist =  dist_buffer)
}

# Filter observations with study extent
study_extent <- terra::intersect(study_extent, terra::vect(q))
obs <- terra::intersect(terra::vect(obs), study_extent)
obs <- as(obs, "Spatial")
raster::crs(obs) <- raster::crs(obs_pres@proj4string)

# Make mesh
pedge <- as.numeric(args[3])
edge <- min(c(diff(raster::bbox(study_extent)[1,])*pedge,diff(raster::bbox(study_extent)[2,])*pedge))
mesh <- INLA::inla.mesh.2d(boundary = study_extent,
                           max.edge = c(edge, edge*5), 
                           min.angle = 21,
                           offset = c(edge, edge*2), 
                           cutoff = edge/2, 
                           crs = raster::crs(study_extent))

# Crop raster with study extent
rast <- raster::raster("data/rast.gri")
rast <- terra::mask(terra::crop(terra::rast(rast), study_extent)
study_extent <- as(study_extent, "Spatial")
raster::crs(study_extent) <- raster::crs(obs_pres@proj4string)
rast <- raster::raster(rast)
raster::crs(rast) <- raster::crs(obs_pres@proj4string)


# Save all four objects
raster::writeRaster(raster, paste0(path_sp(species)$spat,"/rast"))
saveRDS(mesh, paste0(path_sp(species)$spat,"/mesh.rds"))
saveRDS(obs, paste0(path_sp(species)$spat,"/obs.rds"))
saveRDS(study_extent, paste0(path_sp(species)$spat,"/study_extent.rds"))

