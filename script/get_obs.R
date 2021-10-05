#'
#' 
#' 
#' 
#' 
#' 
#' 

get_obs <- function(species) {
  # Connection to DB
  con <- atlasBE::conn(Sys.getenv("user"), Sys.getenv("pwd"), Sys.getenv("host"), Sys.getenv("dbname"))

  # Import data
  ## Get id of species
  id <- RPostgres::dbGetQuery(con, paste0("SELECT id FROM taxa WHERE scientific_name = '",species,"';"))
  obs <- RPostgres::dbGetQuery(con, paste0("SELECT * FROM api.get_bird_presence_absence(",id,");"))
  RPostgres::dbDisconnect(con)

  return(obs)
}
