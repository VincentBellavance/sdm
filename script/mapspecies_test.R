# This script will serve as a test to make species distribution models for the Antigone Canadensis using mapSpecies package.
# The goal here is to make a SDM using only space and time

# Set environment
library(atlasBE)
library(mapSpecies)
library(RPostgres)
library(ratlas)

# Connection to DB
con <- atlasBE::conn(Sys.getenv("user"), Sys.getenv("pwd"), Sys.getenv("host"), Sys.getenv("dbname"))

# Species
sciname = "Antigone canadensis"

# Import data
## Get id of species
id <- RPostgres::dbGetQuery(con, paste0("SELECT id FROM taxa WHERE scientific_name = '",sciname,"';"))
obs <- RPostgres::dbGetQuery(con, paste0("SELECT * FROM api.get_bird_presence_absence(",id,");"))
RPostgres::dbDisconnect(con)
# Filter data
years <- c(1990, 2020)
months <- c(5:9)
#hours <- c(hms::as_hms("05:00:00"), hms::as_hms("09:00:00"))
obs <- obs[!is.na(obs$month_obs) & !is.na(obs$time_obs), ]
obs <- obs[obs$year_obs >= years[1] & obs$year_obs <= years[2], ]
obs <- obs[obs$month_obs >= months[1] & obs$month_obs <= months[2], ]
#obs <- obs[obs$time_obs >= hours[1] & obs$time_obs <= hours[2], ]
obs$occurrence <- ifelse(obs$occurrence == FALSE, 0, 1)

# Deal with geom (extract lon and lat from geom character string)
## Split lon and lat from the rest of the character
re <- "\\(([^()]+)\\)"
geom <- strsplit(
  gsub(re, "\\1", stringr::str_extract_all(obs$geom, re)),
  "\\s+"
)

## Bind the lon and lat to create the df
geom <- do.call(rbind,
  lapply(geom, function(x) {
    cbind(as.numeric(x[1]), as.numeric(x[2]))
  })
)
geom <- as.data.frame(geom)
colnames(geom) <- c("lon", "lat")

## Bind geom and obs together
obs <- cbind(geom, obs)
obs <- obs[,-which(colnames(obs) == "geom")]


# Make spatialPointsDataFrame for the observations
proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
obs <- sp::SpatialPointsDataFrame(geom, obs[,-c(1:2)], proj4string = sp::CRS(proj))


# Make spatialPolygons
if(file.exists("data/spacePoly.rds")) {
  q <- readRDS("data/spacePoly.rds")
} else {
  ## Visualize observations on the map
  states <- c("CA.ON", "CA.MB", "CA.NB", "CA.NF", "CA.NS", "CA.NU", "CA.PE", "CA.QC", "US.VT", "US.ME", "US.NY", "US.NH", "US.MA", "US.RI", "US.CT", "US.PA", "US.NJ")

  ### Get states data
  q1<-raster::getData('GADM', country='CAN', level=1, download = TRUE) # Canada
  q2<-raster::getData('GADM', country='USA', level=1, download = TRUE) # US

  ### Keep states where there is ebird data
  q1 <- q1[q1$HASC_1 %in% states,]
  q2 <- q2[q2$HASC_1 %in% states,]

  ### Bind US and Canada together
  q <- rbind(q1, q2) # This is the region of interest we want to model
  q <- raster::aggregate(q, dissolve = T)
  #pdf("test.pdf")
  #sp::plot(q)
  #sp::plot(obs, add = TRUE, pch = 1, col = "red", cex = 0.5)
  #dev.off()
}

# Make the mesh and convince yourself it's good enough (or save the mesh you did before!)
if(file.exists("data/mesh.rds")) {
  mesh <- readRDS("data/mesh.rds")
} else {
  pedge <- 0.02
  edge <- min(c(diff(raster::bbox(q)[1,])*pedge,diff(raster::bbox(q)[2,])*pedge))
  mesh <- INLA::inla.mesh.2d(boundary = q, 
                              max.edge = c(edge, edge*5), 
                              min.angle = 21,
                              offset = c(edge, edge*2), 
                              cutoff = edge/2, 
                              crs = obs@proj4string)
}

# Import raster template (will be used as a "predictors" with 1 everywhere)
if(file.exists("data/rast.rds")) {
  rast <- readRDS("data/rast.rds")
} else {
  ext <- raster::extent(q)
  rast <- raster::raster(crs = obs@proj4string, ext = ext, resolution = 0.05, vals = 1)
  rast <- raster::mask(rast, q)
  names(rast) <- "none"
}


# Make the explanaMesh and pretend you know what's going on
if(file.exists("data/explana.rds")) {
  explana <- readRDS("data/explana.rds")
} else {
  explana <- mapSpecies::explanaMesh(sPoly = q,
                                     meshSpace = mesh,
                                     meshTime = NULL,
                                     X = rast)
}


# Creates folders for every output
years <- c(1990:1994) # TODO: set the width of the temporal window as an argument / Set the years outside the loop
folder_name <- tolower(gsub(" ", "_", sciname))
folder <- paste0("models/", folder_name)
dir.create(folder)
dir.create(paste0(folder,"/all"))
dir.create(paste0(folder,"/map"))
dir.create(paste0(folder,"/mod"))
dir.create(paste0(folder,"/qc"))
dir.create(paste0(folder,"/auc"))


# Make models
for(i in 0:(2020-years[length(years)])) {

  obs_aggr <- raster::rasterize(obs[obs$year_obs %in% (years+i),], rast, field = "occurrence", fun = "max")
  obs_aggr <- raster::rasterToPoints(obs_aggr)
  obs_aggr <- sp::SpatialPointsDataFrame(obs_aggr[,1:2], as.data.frame(obs_aggr[,3]), proj4string = obs@proj4string)
  names(obs_aggr) <- "occurrence"

  ### Test to make sure observation points are not the same as the centroid of the cells obtained in the obs_aggr
  #pdf("test4.pdf")
  #sp::plot(test[which(test@data == 1),])
  #sp::plot(obs[obs$year_obs %in% (years+i) & obs$occurrence == 1,], pch = 1, col = "red", add = TRUE, cex = 0.5)
  #dev.off()


  # Make model and pray that it makes sense
  mod <- mapSpecies::uniSpace(occurrence ~ none,
                              sPoints = obs_aggr,
                              explanaMesh = explana, 
                              family = "binomial",
                              link = "logit", 
                              control.compute = list(waic = TRUE),
                              verbose = TRUE,
                              num.threads = 30,
                              control.inla=list(int.strategy="ccd")) # with ccd, seems less costly in RAM...
  saveRDS(mod, paste0(folder, "/mod/mod_",mean((years+i)),".rds"))

  # Make map and cry
  ## Take Qc spatialPolygon
  qc <- raster::getData('GADM', country='CAN', level=1, download = TRUE) # Canada
  qc <- qc[qc$NAME_1 == "Québec",]

  ## Map the model
  map <- mapSpace(mod,
                  dims = dim(rast)[2:1],
                  type = "mean",
                  sPoly = qc)
  map <- raster::mask(map, qc)

  map_sd <- mapSpace(mod,
                     dims = dim(rast)[2:1],
                     type = "sd",
                     sPoly = qc)
  map_sd <- raster::mask(map_sd, qc)

  map_025 <- mapSpace(mod,
                      dims = dim(rast)[2:1],
                      type = "0.025quant",
                      sPoly = qc)
  map_025 <- raster::mask(map_025, qc)

  map_975 <- mapSpace(mod,
                      dims = dim(rast)[2:1],
                      type = "0.975quant",
                      sPoly = qc)
  map_975 <- raster::mask(map_975, qc)

  ## Map the model on the entire zone
  map_all <- mapSpace(mod,
                      dims = dim(rast)[2:1],
                      type = "mean",
                      sPoly = q)
  map_all <- raster::mask(map_all, q)

  ## Set crs to raster
  raster::crs(map) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(map_sd) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(map_025) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(map_975) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(map_all) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

  ## Set names of rasters
  names(map) <- paste0("mean_", mean((years+i)))
  names(map_sd) <- paste0("sd_", mean((years+i)))
  names(map_025) <- paste0("ci025_", mean((years+i)))
  names(map_975) <- paste0("ci975_", mean((years+i)))
  names(map_all) <- paste0("mean_all_", mean((years+i)))

  ## Compute AUC with map_all and observation used to make the model
  ### Vector of threshold to get 1/0 from model prediction 
  threshold <- seq(0.05, 0.95, by = 0.05)
  df <- as.data.frame(attributes(mod)$sPoints)
  df[,"pred"] <- raster::extract(map_all, df[,c("x", "y")])
  auc <- data.frame(threshold = threshold, auc = NA)
  for(j in 1:length(threshold)) {
    pred_tmp <- ifelse(df$pred >= threshold[j], 1, 0)
    auc[j, "auc"] <- pROC::auc(pROC::roc(df$occurrence, pred_tmp))
  }
  saveRDS(auc, file.path(folder, paste0("auc/auc_",mean((years+i)),".rds")))
  
  ## Save maps from QC
  if(exists("map_stack")){
    map_stack <- raster::stack(map_stack, map, map_sd, map_025, map_975)
  } else {
    map_stack <- raster::stack(map, map_sd, map_025, map_975)
  }

  ## Save maps for entire zone
  if(exists("map_all_stack")){
    map_all_stack <- raster::stack(map_all_stack, map_all)
  } else {
    map_all_stack <- raster::stack(map_all)
  }

}


# Save rasterStack map_stack and map_all_stack
raster::writeRaster(map_stack, file.path(folder, "/map/maps"))
raster::writeRaster(map_all_stack, file.path(folder, "/map/maps_all"))


# Plot the maps
for(i in 0:(2020-years[length(years)])) {
  
  obs_aggr <- raster::rasterize(obs[obs$year_obs %in% (years+i),], rast, field = "occurrence", fun = "max")
  obs_aggr <- raster::rasterToPoints(obs_aggr)
  obs_aggr <- sp::SpatialPointsDataFrame(obs_aggr[,1:2], as.data.frame(obs_aggr[,3]), proj4string = obs@proj4string)
  names(obs_aggr) <- "occurrence"

  map <- map_stack[[paste0("mean_", mean((years+i)))]]
  map_sd <- map_stack[[paste0("sd_", mean((years+i)))]]
  map_025 <- map_stack[[paste0("ci025_", mean((years+i)))]]
  map_975 <- map_stack[[paste0("ci975_", mean((years+i)))]]
  map_all <- map_all_stack[[paste0("mean_all_", mean((years+i)))]]

  ## Plot the map
  pdf(paste0(folder, "/qc/map",mean((years+i)),".pdf"), width = 12)
  par(mfrow = c(1,2), oma = c(0,0,0,2))
  raster::plot(map, zlim = c(0, 1),
               axes = FALSE, box = FALSE, main = "Mean")
  sp::plot(qc, lwd=0.2, add = TRUE)
  raster::plot(map_sd, zlim = c(0, ceiling(max(map_stack[[paste0("sd_", 1992:2018)]][,,], na.rm = TRUE)*100)/100),
               axes = FALSE, box = FALSE, main = "SD")
  sp::plot(qc, lwd=0.2, add = TRUE)
  dev.off()

  ## Plot the ci map
  pdf(paste0(folder, "/qc/ci_map",mean((years+i)),".pdf"), width = 12)
  par(mfrow = c(1,2), oma = c(0,0,0,2))
  raster::plot(map_025, zlim = c(0, 1),
               axes = FALSE, box = FALSE, main = "CI 0.025")
  sp::plot(qc, lwd=0.2, add = TRUE)
  raster::plot(map_975, zlim = c(0, 1),
               axes = FALSE, box = FALSE, main = "CI 0.975")
  sp::plot(qc, lwd=0.2, add = TRUE)
  dev.off()

  ## Plot the map for entire zone
  pdf(paste0(folder, "/all/map",mean((years+i)),".pdf"))
  raster::plot(map_all, zlim = c(0, 1),
               axes = FALSE, box = FALSE, main = "Mean")
  sp::plot(q, lwd = 0.2, add = TRUE)
  sp::plot(obs_aggr[obs_aggr$occurrence == 1,], add = TRUE, pch = 1, col = "blue", cex = 0.1)
  dev.off()
}
  

# Make a GIF for the QC map (both mean and CI)
imgs <- list.files(paste0(folder, "/qc"), pattern = "^map", full.names = TRUE)
img_list <- lapply(imgs, magick::image_read_pdf, density = 500)
img_joined <- magick::image_join(img_list)
img_animated <- magick::image_animate(img_joined, fps = 5)
magick::image_write(img_animated, path = paste0(folder,"/map.gif"))

imgs <- list.files(paste0(folder, "/qc"), pattern = "ci_map", full.names = TRUE)
img_list <- lapply(imgs, magick::image_read_pdf, density = 500)
img_joined <- magick::image_join(img_list)
img_animated <- magick::image_animate(img_joined, fps = 5)
magick::image_write(img_animated, path = paste0(folder,"/ci_map.gif"))
