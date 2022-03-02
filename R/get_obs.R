#'
#' 
#' 
#' 
#' 
#' 
#' 

get_obs <- function(con, species, years) {

  # Import data
  ## Get id of species
  obs <- RPostgres::dbGetQuery(con, paste0("SELECT * FROM api.get_bird_presence_absence('",species,"')"))
  ## Filter dates
  obs <- obs[obs$year_obs %in% years[1]:years[2],]
  
  return(obs)
}
