# This script will serve as a test to make species distribution models for the Antigone Canadensis using mapSpecies package.
# The goal here is to make a SDM using only space and time

# Set environment
library(atlasBE)
library(mapSpecies)
library(RPostgres)
library(ratlas)

con <- atlasBE::conn(Sys.getenv("user"), Sys.getenv("pwd"), Sys.getenv("host"), Sys.getenv("dbname"))


# Import data
## Get id of Antigone Canadensis
sciname = "Falco peregrinus"
id <- RPostgres::dbGetQuery(con, paste0("SELECT id FROM taxa WHERE scientific_name = '",sciname,"';"))
obs <- RPostgres::dbGetQuery(con, paste0("SELECT * FROM api.get_bird_presence_absence(",id,");"))


# Filter data
years <- c(1990, 2020)
months <- c(5, 8)
hours <- c(hms::as_hms("05:00:00"), hms::as_hms("09:00:00"))
obs <- obs[!is.na(obs$month_obs) & !is.na(obs$time_obs), ]
obs <- obs[obs$year_obs >= years[1] & obs$year_obs <= years[2], ]
obs <- obs[obs$month_obs >= months[1] & obs$month_obs <= months[2], ]
obs <- obs[obs$time_obs >= hours[1] & obs$time_obs <= hours[2], ]
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

# Import raster template (will be used as a "predictors" with 1 everywhere)
ext <- raster::extent(q)
rast <- raster::raster(crs = obs@proj4string, ext = ext, resolution = 0.05, vals = 1)
rast <- raster::mask(rast, q)
names(rast) <- "none"

pedge<-0.02
edge<-min(c(diff(raster::bbox(q)[1,])*pedge,diff(raster::bbox(q)[2,])*pedge))
edge

# Make the mesh and convince yourself it's good enough
mesh <- INLA::inla.mesh.2d(boundary = q, 
                            max.edge = c(edge, edge*5), 
                            min.angle = 21,
                            offset = c(edge,edge*2), 
                            cutoff = edge/2, 
                            crs = obs@proj4string)

## Visualize the shit you just made
#pdf("mesh2.pdf")
#par(mfrow = c(2,1))
#plot(mesh, main = "", asp = 1)
#plot(mesh2, main = "", asp = 1)
#dev.off()


#######################
# Testing mesh
#mesh1 <- INLA::inla.mesh.2d(boundary = q, max.edge = c(1.5, 6), offset = c(0.05, 0.5), cutoff = 0.05, crs = obs@proj4string)
#mesh2 <- INLA::inla.mesh.2d(boundary = q, max.edge = c(1.5, 6), offset = c(0.05, 0.5), cutoff = 0.5, crs = obs@proj4string)
#mesh3 <- INLA::inla.mesh.2d(boundary = q, max.edge = c(1.5, 6), offset = c(0.05, 0.5), cutoff = 1.5, crs = obs@proj4string)
#mesh4 <- INLA::inla.mesh.2d(boundary = q, max.edge = c(1.5, 6), offset = c(0.05, 0.5), cutoff = 5, crs = obs@proj4string)
#
#pdf("mesh_test.pdf")
#par(mfrow = c(2,2))
#plot(mesh1, main = "1", asp = 1)
#plot(mesh2, main = "2", asp = 1)
#plot(mesh3, main = "3", asp = 1)
#plot(mesh4, main = "4", asp = 1)
#dev.off()
#######################


# Make the explanaMesh and pretend you know what's going on
explana <- mapSpecies::explanaMesh(sPoly = q,
                                   meshSpace = mesh,
                                   meshTime = NULL,
                                   X = rast)

# Aggregate observations
years <- c(1990:1994)
folder_name <- tolower(gsub(" ", "_", sciname))
folder <- paste0("models/", folder_name)
dir.create(folder)
dir.create(paste0(folder,"/all"))
dir.create(paste0(folder,"/map"))
dir.create(paste0(folder,"/mod"))
dir.create(paste0(folder,"/qc"))


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
                              num.threads = 30)
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


  ## Set crs to raster
  raster::crs(map) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  raster::crs(map_sd) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

  ## Save maps
  map_stack <- raster::stack(map, map_sd)
  raster::writeRaster(map_stack, file.path(folder, "/map/", paste0("map_",mean((years+i)), ".grd")))

  ## Plot the map
  pdf(paste0(folder, "/qc/map",mean((years+i)),".pdf"), width = 12)
  par(mfrow = c(1,2), oma = c(0,0,0,2))
  raster::plot(map, zlim = c(0, 1),
               axes = FALSE, box = FALSE, main = "Mean")
  sp::plot(qc, lwd=0.2, add = TRUE)
  raster::plot(map_sd, zlim = c(0, 0.5),
               axes = FALSE, box = FALSE, main = "SD")
  sp::plot(qc, lwd=0.2, add = TRUE)
  dev.off()


  ## Map the model on the entire zone
  map_all <- mapSpace(mod,
                      dims = dim(rast)[2:1],
                      type = "mean",
                      sPoly = q)
  map_all <- raster::mask(map_all, q)

  ## Set crs to raster
  raster::crs(map_all) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

  ## Plot the map
  pdf(paste0(folder, "/all/map",mean((years+i)),".pdf"))
  raster::plot(map_all, zlim = c(0, 1),
               axes = FALSE, box = FALSE, main = "Mean")
  sp::plot(q, lwd = 0.2, add = TRUE)
  sp::plot(obs_aggr[obs_aggr$occurrence == 1,], add = TRUE, pch = 1, col = "blue", cex = 0.1)
  dev.off()
}


# Make a GIF for the QC map
imgs <- list.files(paste0(folder, "/qc"), full.names = TRUE)
img_list <- lapply(imgs, magick::image_read_pdf, density = 500)

img_joined <- magick::image_join(img_list)

img_animated <- magick::image_animate(img_joined, fps = 5)

magick::image_write(img_animated, path = paste0(folder,"/map.gif"))
