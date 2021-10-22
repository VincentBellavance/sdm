# This script will serve as a test to make species distribution models for the Antigone Canadensis using mapSpecies package.
# The goal here is to make a SDM using only space and time

sdm_bdi <- function(species, year_start = 1990, year_end = 2020, window_width = 5, buffer = 30, proj = "+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0", threshold = seq(0.05, 0.95, by = 0.05)) {

  # Get observations
  obs <- get_obs(species)

  # Filter data with dates of observations
  obs <- obs[obs$year_obs >= year_start & obs$year_obs <= year_end, ]
  obs <- filter_dates(obs, species, buffer)
  obs$occurrence <- ifelse(obs$occurrence == FALSE, 0, 1)


  # Deal with geom (extract lon and lat from geom character string)
  ## Extract coordinates from geom column
  coords <- extract_coords(obs$geom)
  ## Bind coords and obs together
  obs <- cbind(coords, obs)
  ## Remove geom column
  obs <- obs[,-which(colnames(obs) == "geom")]

  # Make spatialPointsDataFrame for the observations
  obs <- sp::SpatialPointsDataFrame(coords, obs[,-c(1:2)], proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  obs <- sp::spTransform(obs, sp::CRS(proj))

  # Make spatialPolygons to use it as boundary for the mesh
  q <- prep_spat_poly(proj)

  # Make spatialPolygons of Qc province to make map
  qc <- prep_qc_spat_poly(proj)

  # Make the mesh and convince yourself it's good enough (or save the mesh you did before!)
  mesh <- prep_mesh(q, obs)

  # Import raster template (will be used as a "predictors" with 1 everywhere)
  rast <- prep_rast_pred(q, obs)

  # Make the explanaMesh and pretend you know what's going on
  explana <- prep_explana(q, mesh, rast)

  # Creates folders for every output
  folder <- prep_folders(year_start, window_width, species)

  # First years
  years <- year_start:(year_start+window_width-1)

  # Make models for every time windows
  for(i in 0:(year_end-years[length(years)])) {
    
    year <- mean((years+i))

    obs_aggr <- aggregate_obs(obs, years, i, rast)

    # Make model and pray that it makes sense
    mod <- mapSpecies::uniSpace(occurrence ~ none,
                                sPoints = obs_aggr,
                                explanaMesh = explana, 
                                family = "binomial",
                                link = "logit", 
                                control.compute = list(waic = TRUE),
                                num.threads = 30,
                                control.inla=list(int.strategy="eb")) # with ccd or eb, seems less costly in RAM...
    saveRDS(mod, paste0(folder, "/mod/mod_",year,".rds"))

    
    ## Map the model to make compute the AUC
    map <- make_map("mean_all", mod, rast, q, year, proj)
    #calc_auc(mod, map, threshold, folder, year)

    ## Stack the map to save it at the end of the loop
    if(exists("map_stack")){
      map_stack <- raster::stack(map_stack, map)
    } else {
      map_stack <- raster::stack(map)
    }

    if(i == year_end-years[length(years)]) {
      raster::writeRaster(map_stack, file.path(folder, "/map/mean"))
    }

  }
}

#
#
#
#
#
#

### Plot the map
#pdf(paste0(folder, "/qc/map_",year,".pdf"), width = 18)
#if(plot_ci == TRUE) {
#  par(mfrow = c(1,3), oma = c(0,0,0,2))
#  # Plot mean map
#  plot_map(map,"Mean")
#  # Plot CI
#  plot_map(map_025, "CI 0.025")
#  plot_map(map_975, "CI 0.975")
#} else {
#  # If plot_ci == FALSE, plot mean only
#  plot_map(map,"Mean")
#}
#dev.off()

### Plot the map for entire zone
#pdf(paste0(folder, "/all/map",mean((years+i)),".pdf"))
#raster::plot(map_all, zlim = c(0, 1),
#             axes = FALSE, box = FALSE, main = "Mean")
#sp::plot(q, lwd = 0.2, add = TRUE)
#dev.off()
