#'
#' 
#' 
#' 
#' 
#' 
#' 

get_obs <- function(species, years) {
  # Connection to DB
  con <- atlasBE::conn(Sys.getenv("user"), Sys.getenv("pwd"), Sys.getenv("host"), Sys.getenv("dbname"))

  # Import data
  ## Get id of species
  id <- RPostgres::dbGetQuery(con, paste0("SELECT id FROM api.bird_quebec_taxa_ref WHERE scientific_name = '",species,"';"))[1,"id"]
  obs <- RPostgres::dbGetQuery(con, paste0("SELECT * FROM api.get_bird_presence_absence(",id,") WHERE year_obs BETWEEN ",years[1]," AND ",years[2],";"))
  RPostgres::dbDisconnect(con)

  return(obs)
}
