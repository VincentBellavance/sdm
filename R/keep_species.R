#'
#' 
#'
#'
#'
#'
#'

keep_species <- function(obs, species) {

  sp_code <- get_sp_code(species)

  if(!is.null(sp_code)) {
    species_info <- ebird_scraping(sp_code)

    qc <- obs_qc(obs)
    qc_mo <- table(qc$month_obs)
    summer <- which(max(qc_mo) == qc_mo) %in% c(5:9)

    if("Breeding season" %in% species_info & summer) {
      if(species_info[which(species_info %in% "Breeding season")+1] == "Not shown") {
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else if(length(species_info) == 0 & summer) {
      return(TRUE)
    } else if("Year-round" %in% species_info) {
      return(TRUE)
    } else {
      line <- gsub("_", " ", species) |>
                stringr::str_to_sentence() |>
                  paste0(";",paste0(qc_mo, collapse=";"))
      write(line, file = "rm_sp_month.csv", append = TRUE, sep = "\n")
      write(paste0(species,";Max obs is not in summer"), file = "removed_species.csv", append = TRUE, sep = "\n")
      return(FALSE)
    }
  } else {
    line <- gsub("_", " ", species) |>
          stringr::str_to_sentence() |>
            paste0(";",paste0(rep("", 12), collapse=";"))
    write(line, file = "rm_sp_month.csv", append = TRUE, sep = "\n")
    write(paste0(species,";Could not find info on species"), file = "removed_species.csv", append = TRUE, sep = "\n")
    return(FALSE)
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