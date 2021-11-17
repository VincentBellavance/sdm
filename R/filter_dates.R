#'
#' 
#' 
#' 
#' 
#' 

filter_dates <- function(obs, species, buffer) {

  mig_dates <- suppressMessages(time_interval(species, buffer))
  if(!is.null(mig_dates)) {
    obs <- obs[!is.na(obs$month_obs) & !is.na(obs$day_obs), ]
    tmp_md <- paste0("2020-",obs$month_obs, "-", obs$day_obs)
    obs[,"tmp_md"] <- format(as.Date(tmp_md , format = "%Y-%m-%d") - buffer, "%m-%d")
    obs <- obs[!is.na(obs$tmp_md),]
    obs <- obs[obs$tmp_md >= mig_dates[1] & obs$tmp_md <= mig_dates[2], ]
    obs <- obs[,-grep("tmp_md", colnames(obs))]
  }

  return(obs)

}


#'
#' 
#' 
#' 
#' 
#' 

time_interval <- function(species, buffer) {

  ebird_taxa <- c(rebird::ebirdtaxonomy("species")[,"sciName"])[[1]]
  if(species$accepted %in% ebird_taxa) {
  
    species_info <- ebird_scraping(species$accepted)

    if ("Breeding season" %in% species_info) {
      return(get_migration_dates(species_info, buffer))
    } else {
      return(NULL)
    }
  
  } else {

    syn_in_ebird <- species$synonym %in% ebird_taxa
    
    if(any(syn_in_ebird)) {
    
      synonym <- species$synonym[syn_in_ebird]
    
      species_info <- ebird_scraping(synonym)
    
      if ("Breeding season" %in% species_info) {
        return(get_migration_dates(species_info, buffer))
      } else {
        return(NULL)
      }
    } else {
      return(c("06-01", "09-01"))
    }
  }
}


#'
#' 
#' 
#' 
#' 
#' 

ebird_scraping <- function(species) {
  
  # Get species code to paste in the url
  sp_code <- rebird::species_code(species)
  
  url <- paste0("https://ebird.org/science/status-and-trends/",sp_code,"/abundance-map")
  tmp <- rvest::read_html(url)
  tmp <- rvest::html_nodes(tmp, ".VisProduct-meta-main")
  tmp <- rvest::html_text(tmp)

  cleaned <- strsplit(tmp, split = "\t")[[1]]
  cleaned <- gsub("\n", "", cleaned)
  cleaned <- cleaned[-which(nchar(cleaned) %in% c(0,1))]

  return(cleaned)

}


#'
#' 
#' 
#' 
#' 
#' 

get_migration_dates <- function(species_info, buffer) {
    
  predate <- strsplit(species_info[grep("Pre-breeding migratory season", species_info)+1], split = " - ")[[1]][2]
  postdate <- strsplit(species_info[grep("Post-breeding migratory season", species_info)+1], split = " - ")[[1]][1]

  predate <- format(as.Date(predate , format = "%b %d") - buffer, "%m-%d")
  postdate <- format(as.Date(postdate , format = "%b %d") + buffer, "%m-%d")

  return(c(predate, postdate))

}
