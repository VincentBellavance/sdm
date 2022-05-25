#' Title
#' 
#' Description
#'
#' @param
#'
#' @return

make_ts <- function(sdms) {
  
  pocc <- data.frame(years = names(sdms), pocc = NA)
  for(i in names(sdms)) {
   pocc[pocc$years == i, "pocc"] <- sum(sdms[[i]][,,], na.rm = T)
  }
  pocc$years <- as.integer(gsub("[^0-9.-]", "", pocc$years))

  pocc <- pocc[order(pocc$years),]

  return(pocc)

}


#' Title
#' 
#' Description
#'
#' @param
#'
#' @return 

compare_to_mean <- function(pocc) {
  
  comp_mean <- sapply(1:nrow(pocc), function(x) {
    filter <- c(x-2,x-1,x+1,x+2)
    filter <- filter[filter > 0 & filter < nrow(pocc)]
    pocc[x,"pocc"]/mean(pocc[filter,"pocc"], na.rm = T)
  })
  comp_mean <- comp_mean/sum(comp_mean)
  return(comp_mean)
  
}