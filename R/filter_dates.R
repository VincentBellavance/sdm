#'
#' 
#' 
#' 
#' 
#' 

filter_dates <- function(con, obs, species, buffer) {

  mig_dates <- suppressMessages(time_interval(con, species, buffer))
  
  if(!is.null(mig_dates)) {
    obs <- obs[!is.na(obs$month_obs) & !is.na(obs$day_obs), ]
    obs[,"tmp_md"] <- format(as.Date(paste0(obs$month_obs,"-",obs$day_obs) , format = "%m-%d") - buffer, "%m-%d")
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

time_interval <- function(con, species, buffer) {

  sp_code <- get_sp_code(species$accepted, con)
  
  if(length(sp_code) > 0) {
  
    species_info <- ebird_scraping(sp_code)

    if(length(species_info) == 0) {
      return(c("06-01", "09-01"))
    } else if("Year-round" %in% species_info) {
      return(NULL)
    } else if ("Breeding season" %in% species_info) {
      return(get_migration_dates(species_info, buffer))
    } else {
      return(NULL)
    }
  }
}


#'
#' 
#' 
#' 
#' 
#' 

ebird_scraping <- function(sp_code) {
  
  # Get URL and content of the page  
  url <- paste0("https://ebird.org/science/status-and-trends/",sp_code,"/abundance-map")
  tmp <- rvest::read_html(url)
  tmp <- rvest::html_nodes(tmp, ".VisProduct-meta-main")
  tmp <- rvest::html_text(tmp)

  if(length(tmp) != 0) {
    # Clean the content
    cleaned <- strsplit(tmp, split = "\t")[[1]]
    cleaned <- gsub("\n", "", cleaned)
    cleaned <- cleaned[-which(nchar(cleaned) %in% c(0,1))]
  
    return(cleaned)
  } else {
    return(NULL)
  }

  
}


#'
#' 
#' 
#' 
#' 
#' 

get_migration_dates <- function(species_info, buffer) {
  
  if(species_info[grep("Pre-breeding migratory season", species_info)+1] == "Not shown") {
    breed <- strsplit(species_info[grep("Breeding season", species_info)+1], split = " - ")[[1]][1]
    predate <- format(as.Date(breed, format = "%b %d") -7, "%b %d")
  } else {
    predate <- strsplit(species_info[grep("Pre-breeding migratory season", species_info)+1], split = " - ")[[1]][2]
  }

  if(species_info[grep("Post-breeding migratory season", species_info)+1] == "Not shown") {
    breed <- strsplit(species_info[grep("Breeding season", species_info)+1], split = " - ")[[1]][2]
    postdate <- format(as.Date(breed, format = "%b %d") + 7, "%b %d")
  } else {
    postdate <- strsplit(species_info[grep("Post-breeding migratory season", species_info)+1], split = " - ")[[1]][1]
  }
  
  Sys.setlocale("LC_ALL","en_US.UTF-8")
  predate <- format(as.Date(predate , format = "%b %d") - buffer, "%m-%d")
  postdate <- format(as.Date(postdate , format = "%b %d") + buffer, "%m-%d")

  return(c(predate, postdate))

}


#'
#' 
#' 
#' 
#' 
#' 

get_sp_code <- function(species, con) {
  
  ebird_taxa <- as.data.frame(rebird::ebirdtaxonomy("species"))

  sp_code <- ebird_taxa[ebird_taxa$sciName %in% species, "speciesCode"]

  if(length(sp_code) == 0) {
    com <- as.vector(unlist(RPostgres::dbGetQuery(con, paste0("SELECT vernacular_en FROM api.taxa WHERE valid_scientific_name = '",species,"';"))))
    if(com == "fandhill crane") com <- "sandhill crane"
    sp_code <- ebird_taxa[tolower(ebird_taxa$comName) %in% com, "speciesCode"]
  }

  if(length(sp_code) == 0) {
    return(NULL)
  } else {
    return(sp_code)
  }

}