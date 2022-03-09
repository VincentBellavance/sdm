#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

plot_map <- function(file, map, title, region) {

  png(file, width = 960, height = 960)
  raster::plot(map, 
               zlim = c(0, 1),
               axes = FALSE, 
               box = FALSE, 
               main = title)
  sp::plot(qc, 
           lwd=0.2, 
           add = TRUE)
  dev.off()
}



### Plot the map
#pdf(paste0(folder, "/qc/map_",year,".pdf"), width = 18)
#if(plot_ci == TRUE) {
#  par(mfrow = c(1,3), oma = c(0,0,0,2))
#  # Plot mean map
#  plot_map(map,"Mean")
#  # Plot CI
#  plot_map(map_025, "CI 0.025")
#  plot_map(map_975, "CI 0.975")
#} else {
#  # If plot_ci == FALSE, plot mean only
#  plot_map(map,"Mean")
#}
#dev.off()

### Plot the map for entire zone
#pdf(paste0(folder, "/all/map",mean((years+i)),".pdf"))
#raster::plot(map_all, zlim = c(0, 1),
#             axes = FALSE, box = FALSE, main = "Mean")
#sp::plot(q, lwd = 0.2, add = TRUE)
#dev.off()
