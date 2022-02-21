#'
#' 
#'
#'
#'
#'
#'

keep_species <- function(obs) {

  species_info <- ebird_scraping(species$accepted)
  
  if("Breeding season" %in% species_info) {
    qc <- obs_qc(obs)
    qc_mo <- table(qc$month_obs)
    return(which(max(qc_mo) == qc_mo) %in% c(5:9))
  } else {
    return(TRUE)
  }
  
}


#'
#' 
#'
#'
#'
#'

obs_qc <- function(obs, lat1=44.99, lat2=62.59, lon1=-79.76, lon2=-57.10) {

  # Keep observation inside qc lat and lon and occurrence = 1
  return(obs[obs$lat > lat1 & 
             obs$lat < lat2 & 
             obs$lon > lon1 &
             obs$lon < lon2 &
             obs$occurrence,])

}