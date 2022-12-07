#' Title
#' 
#' Description
#'
#' @param
#'
#' @return

path <- function(zone, output_dir, obs_dir) {
  list(
    spat=paste0(output_dir,"/spatial/",zone,"/"),
    mod=paste0(output_dir,"/models/",zone,"/"),
    log=paste0(output_dir,"/log/",zone,"/"),
    stack=paste0(output_dir,"/stack/",zone,"/"),
    maps=paste0(output_dir,"/maps/",zone,"/"),
    obs=obs_dir
  )
}


path_sp <- function(species,
                    output_dir = "/home/belv1601/scratch/output", 
                    obs_dir = "data/occurrences",
                    zone="quebec") {
  list(
    obs=paste0(path(zone, output_dir, obs_dir)$obs, "/",species, ".rds"),
    spat=paste0(path(zone, output_dir, obs_dir)$spat, species),
    mod=paste0(path(zone, output_dir, obs_dir)$mod, species),
    log=paste0(path(zone, output_dir, obs_dir)$log, species),
    stack=paste0(path(zone, output_dir, obs_dir)$stack, species),
    check=paste0(path(zone, output_dir, obs_dir)$check, species),
    maps=paste0(path(zone, output_dir, obs_dir)$maps, species)
  )
}
