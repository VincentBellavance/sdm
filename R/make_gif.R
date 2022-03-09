#'
#' 
#'
#' 
#' 
#' 
#' 
#' 

make_gif <- function(folder, region) {

  # Make a GIF for the QC map (both mean and CI)
  imgs <- list.files(paste0(folder, region), full.names = TRUE)
  img_list <- lapply(imgs, magick::image_read)
  img_joined <- magick::image_join(img_list)
  img_animated <- magick::image_animate(img_joined, fps = 2)
  magick::image_write(img_animated, path = paste0(folder,"/map_",region,".gif"))

}