#' Title
#' 
#' Description
#'
#' @param
#'
#' @return

path <- function(zone) {
  list(
    spat=paste0("~/scratch/output/spatial/",zone,"/"),
    mod=paste0("~/scratch/output/models/",zone,"/"),
    log=paste0("~/scratch/output/log/",zone,"/"),
    stack=paste0("~/scratch/output/stack/",zone,"/"),
    check=paste0("~/scratch/output/checks/",zone,"/"),
    maps=paste0("~/scratch/output/maps/",zone,"/"),
    occ="~/projects/def-dgravel/belv1601/sdm/data/occurrences/"
  )
}


path_sp <- function(species, zone="quebec") {
  list(
    occ=paste0(path(zone)$occ, species, ".rds"),
    spat=paste0(path(zone)$spat, species),
    mod=paste0(path(zone)$mod, species),
    log=paste0(path(zone)$log, species),
    stack=paste0(path(zone)$stack, species),
    check=paste0(path(zone)$check, species),
    maps=paste0(path(zone)$maps, species)
  )
}
