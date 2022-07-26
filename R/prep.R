#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_spat_poly <- function(proj) {

  # Make spatialPolygons
  if(file.exists("data/spacePoly.rds")) {
    return(readRDS("data/spacePoly.rds"))
  } else {
    ## Visualize observations on the map
    us_states <- c("US.VT", "US.ME", "US.NY", "US.NH", "US.MA", "US.RI", "US.CT", "US.PA", "US.NJ")
    can_provinces <- c("CA.NB", "CA.NS", "CA.PE", "CA.QC")
    lab_divisions <- c("CA.NF.TE")
    on_districts_rm <- c("CA.ON.KR", "CA.ON.RR", "CA.ON.TB", "CA.ON.WB")

    ### Get states data
    can1<-raster::getData('GADM', country='CAN', level=1, download = TRUE)
    can2 <- raster::getData('GADM', country='CAN', level=2, download = TRUE)
    us<-raster::getData('GADM', country='USA', level=1, download = TRUE)

    ### Keep states where there is ebird data
    can1 <- can1[can1$HASC_1 %in% can_provinces,1]
    on <- can2[can2$NAME_1 == "Ontario" & !can2$HASC_2 %in% on_districts_rm,1]
    lab <- can2[can2$HASC_2 %in% lab_divisions,1]
    us <- us[us$HASC_1 %in% us_states,1]

    ### Bind US and Canada together
    q <- rbind(can1,on,lab,us) # This is the region of interest we want to model
    q <- terra::vect(q)
    q <- terra::aggregate(q, dissolve = T)
    q <- terra::project(q, proj)
    q <- as(sf::st_as_sf(q), "Spatial")
    raster::crs(q) <- proj

    saveRDS(q, "data/spacePoly.rds")
    return(q)
  }
}


#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_mesh <- function(q, proj) {

  # Make the mesh and convince yourself it's good enough (or read the mesh you did before!)
  if(file.exists("data/mesh.rds")) {
    return(readRDS("data/mesh.rds"))
  } else {
    pedge <- 0.035
    edge <- min(c(diff(raster::bbox(q)[1,])*pedge,diff(raster::bbox(q)[2,])*pedge))
    mesh <- INLA::inla.mesh.2d(boundary = q, 
                                max.edge = c(edge, edge*5), 
                                min.angle = 21,
                                offset = c(edge, edge*2), 
                                cutoff = edge/2, 
                                crs = sp::CRS(proj))
    saveRDS(mesh, "data/mesh.rds")
    return(mesh)
  }
}


#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_rast_pred <- function(q, proj, res) {

  if(file.exists("data/rast.gri")) {
    return(raster::stack("data/rast.gri"))
  } else {
    ext <- raster::extent(q)
    rast <- raster::raster(crs = sp::CRS(proj), ext = ext, resolution = res, vals = 1)
    rast <- terra::crop(terra::rast(rast), terra::vect(q))
    rast <- terra::mask(rast, terra::vect(q))
    rast <- raster::raster(rast)
    names(rast) <- "none"
    raster::writeRaster(rast, "data/rast")
    return(rast)
  }
}


#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_explana <- function(q, mesh, rast) {

  if(file.exists("data/explana.rds")) {
    return(readRDS("data/explana.rds"))
  } else {
    explana <- mapSpecies::explanaMesh(sPoly = q,
                                       meshSpace = mesh,
                                       meshTime = NULL,
                                       X = rast)
    saveRDS(explana, "data/explana.rds")
    return(explana)
  }
}


#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_qc_spat_poly <- function(proj) {

  # Make spatialPolygons
  if(file.exists("data/qc_spacePoly.rds")) {
    return(readRDS("data/qc_spacePoly.rds"))
  } else {
    qc <- raster::getData('GADM', country='CAN', level=1, download = TRUE) # Canada
    qc <- qc[qc$NAME_1 == "Québec",] 
    qc <- sp::spTransform(qc, sp::CRS(proj))

    # Save it to not make it everytime
    saveRDS(qc, "data/qc_spacePoly.rds")
    return(qc)
  }
}
