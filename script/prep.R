#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_spat_poly <- function() {

  # Make spatialPolygons
  if(file.exists("data/spacePoly.rds")) {
    return(readRDS("data/spacePoly.rds"))
  } else {
    ## Visualize observations on the map
    states <- c("CA.ON", "CA.MB", "CA.NB", "CA.NF", "CA.NS", "CA.NU", "CA.PE", "CA.QC", "US.VT", "US.ME", "US.NY", "US.NH", "US.MA", "US.RI", "US.CT", "US. PA", "US.NJ")

    ### Get states data
    q1<-raster::getData('GADM', country='CAN', level=1, download = TRUE) # Canada
    q2<-raster::getData('GADM', country='USA', level=1, download = TRUE) # US

    ### Keep states where there is ebird data
    q1 <- q1[q1$HASC_1 %in% states,]
    q2 <- q2[q2$HASC_1 %in% states,]

    ### Bind US and Canada together
    q <- rbind(q1, q2) # This is the region of interest we want to model
    q <- raster::aggregate(q, dissolve = T)
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

prep_mesh <- function(q, obs) {

  # Make the mesh and convince yourself it's good enough (or read the mesh you did before!)
  if(file.exists("data/mesh.rds")) {
    return(readRDS("data/mesh.rds"))
  } else {
    pedge <- 0.02
    edge <- min(c(diff(raster::bbox(q)[1,])*pedge,diff(raster::bbox(q)[2,])*pedge))
    mesh <- INLA::inla.mesh.2d(boundary = q, 
                                max.edge = c(edge, edge*5), 
                                min.angle = 21,
                                offset = c(edge, edge*2), 
                                cutoff = edge/2, 
                                crs = obs@proj4string)
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

prep_rast_pred <- function(q, obs) {

  if(file.exists("data/rast.rds")) {
    return(readRDS("data/rast.rds"))
  } else {
    ext <- raster::extent(q)
    rast <- raster::raster(crs = obs@proj4string, ext = ext, resolution = 0.05, vals = 1)
    rast <- raster::mask(rast, q)
    names(rast) <- "none"
    saveRDS(rast, "data/rast.rds")
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
    saveRDS("data/explana.rds")
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

prep_folders <- function(year_start, window_width, species) {

  years <- seq(from = year_start, by = 1, length.out = window_width)
  folder_name <- tolower(gsub(" ", "_", species))
  folder <- paste0("models/", folder_name)
  dir.create(folder)
  dir.create(paste0(folder,"/all"))
  dir.create(paste0(folder,"/map"))
  dir.create(paste0(folder,"/mod"))
  dir.create(paste0(folder,"/qc"))
  dir.create(paste0(folder,"/auc"))

}


#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

prep_qc_spat_poly <- function() {

  # Make spatialPolygons
  if(file.exists("data/qc_spacePoly.rds")) {
    return(readRDS("data/qc_spacePoly.rds"))
  } else {
    qc <- raster::getData('GADM', country='CAN', level=1, download = TRUE) # Canada
    qc <- qc[qc$NAME_1 == "Québec",] # Save it to not make it everytime
    saveRDS("data/qc_spacePoly.rds")
    return(qc)
  }
}
