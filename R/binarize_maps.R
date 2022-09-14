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

find_threshold <- function(sdm, obs, type = "sensitivity") {
  if (!inherits(sdm, "SpatRaster")) {
    sdm <- terra::rast(sdm)
  }
  if (!inherits(obs, "SpatVect")) {
    obs <- terra::vect(obs)
  }

  if (type == "spse") {
    type <- "spec_sens"
  }
  # extract model estimated suitability for occurrence localities
  # occs_vals <- terra::extract(sdm, occs)
  occs_vals <- terra::extract(sdm, obs[obs$occurrence == 1,]) %>%
    dplyr::select(-ID) %>%
    dplyr::pull(1)

  # extract model estimated suitability for background
  bg_vals <- terra::extract(sdm, obs[obs$occurrence == 0,]) %>%
    dplyr::select(-ID) %>%
    dplyr::pull(1)

  # evaluate predictive ability of model
  ev <- dismo::evaluate(occs_vals, bg_vals)
  # detect possible thresholds
  thr.table <- dismo::threshold(ev)
  thresh <- thr.table[, type]

  return(thresh)
}

#' @title Binarize a prediction raster
#' 
#' @name binarize_pred
#' @param obs data frame, containing the coordinates to reproject
#' @param predictors, raster
#' @param lon string, name of the longitude column (same projection as predictor raster)
#' @param lat string, name of the latitude column (same projection as predictor raster)
#' @param proj character, initial projection of the xy coordinates
#' @return spatial points
#' @import dplyr dismo raster
#' @export

binarize_pred <- function(sdm, threshold, ci025) {
  sdm[sdm >= threshold] <- 1
  sdm[sdm < threshold & ci025 <= 0.001] <- 0
  return(sdm)
}