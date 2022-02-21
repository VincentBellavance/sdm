#'
#' 
#' 
#' 
#' 
#' 
#' 

get_obs <- function(species, years) {

  # Import data
  ## Get id of species
  obs <- ratlas::get_bird_presence_absence(taxa_name=species)
  print("yay!")
  ## Filter dates
  obs <- obs[obs$year_obs %in% years[1]:years[2],]
  
  return(obs)
}
