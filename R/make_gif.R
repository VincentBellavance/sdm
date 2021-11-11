#'
#' 
#'
#' 
#' 
#' 
#' 
#' 

#make_gif() {
#
#  # Make a GIF for the QC map (both mean and CI)
#  imgs <- list.files(paste0(folder, "/qc"), pattern = "^map", full.names = TRUE)
#  img_list <- lapply(imgs, magick::image_read_pdf, density = 500)
#  img_joined <- magick::image_join(img_list)
#  img_animated <- magick::image_animate(img_joined, fps = 5)
#  magick::image_write(img_animated, path = paste0(folder,"/map.gif"))
#
#}