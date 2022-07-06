#' Function to convert observations in observations (trials) and success (presence)
#' 
#' Description
#'
#' @param
#'
#' @return 

convert_in_success_trials <- function(obs, rast) {

  # Number of trials
  observations <- pointcount(obs, rast)

  # Number of success (presence)
  presences <- pointcount(obs[obs$occurrence == 1,], rast)

  # Make a stack with both
  obs <- raster::stack(observations, presences)
  names(obs) <- c("observations", "presences")

  # Transform in spdf
  return(as(obs, "SpatialPointsDataFrame"))
}



#' Calculate number of points in every cell 
#' 
#' Description
#'
#' @param
#'
#' @return

pointcount <- function(obs, rast){
  # make a raster of zeroes like the input
  r = rast
  r[] = 0
  # get the cell index for each point and make a table:
  counts = table(raster::cellFromXY(rast,obs))
  # fill in the raster with the counts from the cell index:
  r[as.numeric(names(counts))] = counts
  return(r)
}