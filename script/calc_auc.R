#'
#' 
#' 
#' 
#' 
#' 
#' 
#' 

calc_auc <- function(mod, map, threshold, folder, year) {

  ### Using vector of threshold to get 1/0 from model prediction 
  df <- as.data.frame(attributes(mod)$sPoints)
  df[,"pred"] <- raster::extract(map_all, df[,c("x", "y")])
  auc <- data.frame(threshold = threshold, auc = NA)
  for(j in 1:length(threshold)) {
    pred_tmp <- ifelse(df$pred >= threshold[j], 1, 0)
    auc[j, "auc"] <- pROC::auc(pROC::roc(df$occurrence, pred_tmp))
  }
  saveRDS(auc, file.path(folder, paste0("auc/auc_",year,".rds")))

}