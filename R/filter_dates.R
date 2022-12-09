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

  sp_code <- get_sp_code(species, con)
  
  if(length(sp_code) > 0) {
    ebird_scraping(sp_code)
  } else {
    return(c("05-15", "10-01"))
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
  tmp_primary <- rvest::html_nodes(tmp, ".StVizLegendItem-label-primary")
  tmp_primary <- rvest::html_text(tmp_primary)

  if("Year-round" %in% tmp_primary & 
     !"Breeding season" %in% tmp_primary) return(NULL)

  if(!"Year-round" %in% tmp_primary & 
     !"Breeding season" %in% tmp_primary) return(c("05-01", "10-01"))
  
  if("Breeding season" %in% tmp_primary) {
    
    tmp_primary <- tmp_primary[grep("season", 
                                    tmp_primary)]
    tmp_secondary <- rvest::html_nodes(tmp, ".StVizLegendItem-label-secondary")
    tmp_secondary <- rvest::html_text(tmp_secondary)
    tmp_secondary <- tmp_secondary[grep("abundance", 
                                        tmp_secondary, 
                                        invert = TRUE)]

    dates <- data.frame(period = tmp_primary,
                        dates = tmp_secondary)

    range <- strsplit(dates[dates$period == "Breeding season", "dates"],
                      split = " - ")[[1]]

    return(c(format(as.Date(range[1], format = "%d %b") - buffer, "%m-%d"),
             format(as.Date(range[2], format = "%d %b") + buffer, "%m-%d")))
  }

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