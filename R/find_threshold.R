#' @title Find threshold to binarize an sdm output
#' 
#' @name find_threshold
#' @param sdm raster, sdm prediction
#' @param occs dataframe, presence points
#' @param bg dataframe, background points
#' @param type string, method to choose the threshold, one of "mtp", "kappa", "spec_sens", "no_omission", "prevalence", "equal_sens_spec", "sensitivity"
#' @return sa float
#' @import dplyr dismo raster
#' @export

find_threshold <- function(sdm, obs, type = "sensitivity", sensitivity = 0.99) {

  if (type == "spse") {
    type <- "spec_sens"
  }
  # extract model estimated suitability for occurrence localities
  # occs_vals <- terra::extract(sdm, occs)
  occs_vals <- raster::extract(sdm, obs[obs$occurrence == 1,])

  # extract model estimated suitability for background
  bg_vals <- raster::extract(sdm, obs[obs$occurrence == 0,])

  # evaluate predictive ability of model
  ev <- dismo::evaluate(occs_vals, bg_vals)
  # detect possible thresholds
  thr.table <- dismo::threshold(ev, sensitivity = sensitivity)
  threshold <- thr.table[, type]

  return(threshold)
}
