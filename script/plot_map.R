#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

plot_map <- function(map, title) {
  raster::plot(map, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = title)
  sp::plot(qc, 
           lwd=0.2, 
           add = TRUE)
}