#' Title
#' 
#' Description
#'
#' @param
#'
#' @return

path <- function() {
  list(
    spat="~/scratch/output/spatial/",
    mod="~/scratch/output/models/",
    log="~/scratch/output/log/",
    stack="~/scratch/output/stack/",
    check="~/scratch/output/checks/",
    maps="~/scratch/output/maps/",
    occ="~/projects/def-dgravel/belv1601/sdm/data/occurrences/"
  )
}


#' Title
#' 
#' Description
#'
#' @param
#'
#' @return

path_sp <- function(species) {
  list(
    occ=paste0(path()$occ, species, ".rds"),
    spat=paste0(path()$spat, species),
    mod=paste0(path()$mod, species),
    log=paste0(path()$log, species),
    stack=paste0(path()$stack, species),
    check=paste0(path()$check, species),
    maps=paste0(path()$maps, species)
  )
}

#' Title
#' 
#' Description
#'
#' @param
#'
#' @return 

path_maps <- function(species, maps_dir = path()$maps) {
  list(
    qc = paste0(path_sp(species)$maps, "/qc"),
    region = paste0(path_sp(species)$maps, "/region"),
    region_pres = paste0(path_sp(species)$maps, "/region_pres"),
    region_abs = paste0(path_sp(species)$maps, "/region_abs")
  ) 
}
