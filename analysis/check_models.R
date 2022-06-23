# Script to check if any models did not converged / if changes are smooth throught time
# Decisions are taken with ...

args = commandArgs(trailingOnly=TRUE)
library(raster)
source("R/make_ts.R")
source("R/path.R")

# The folder that contains maps
species <- args[1]
stack_exists <- "maps.gri" %in% list.files(path_sp(species)$maps)

if(!stack_exists) {
  dir.create(path_sp(species)$check)
  saveRDS("", paste0(path_sp(species)$check, "/not_completed.rds"))
} else {

  # Read raster stack
  sdms <- raster::stack(paste0(path_sp(species)$maps, "/maps.gri"))

  # Calculate sum of p(occ) for every year
  pocc <- make_ts(sdms)

  # Check for missing years
  missing_years <- c(1992:2018)[!1992:2018 %in% pocc$years]

  # Compare the first value of a 2 (or 3,4,5) value subset to the standard deviation of this subset
  # Raise a flag if the value represent 0.1 of the total standard deviation or more (explain this better)
  ## Threshold for decision
  mean_threshold <- 0.075
  sd_threshold1 <- 0.075
  sd_threshold2 <- 0.075
  sd_threshold3 <- 0.075

  ## Tests
  sd1 <- zoo::rollapply(pocc$pocc, 2, sd)/sum(zoo::rollapply(pocc$pocc, 2, sd))
  sd2 <- zoo::rollapply(pocc$pocc, 3, sd)/sum(zoo::rollapply(pocc$pocc, 3, sd))
  sd3 <- zoo::rollapply(pocc$pocc, 4, sd)/sum(zoo::rollapply(pocc$pocc, 4, sd)) 
  comp_mean <- compare_to_mean(pocc)

  ## Decision
  sus_mean <- which(comp_mean>mean_threshold)
  sus_sd1 <- which(sd1>=sd_threshold1)
  sus_sd2 <- which(sd2>=sd_threshold2)
  sus_sd3 <- which(sd3>=sd_threshold3)

  ## Compare results with mean and sd
  sus_mean_sd1 <- pocc[sus_sd1[sus_sd1 %in% sus_mean], "years"]
  sus_mean_sd2 <- pocc[sus_sd2[sus_sd2 %in% sus_mean], "years"]
  sus_mean_sd3 <- pocc[sus_sd3[sus_sd3 %in% sus_mean], "years"]
  sus_mean_only <- pocc[sus_mean[!sus_mean %in% unique(c(sus_sd1,sus_sd2))], "years"]
  sus_all <- pocc[sus_mean[sus_mean %in% unique(c(sus_sd1,sus_sd2, sus_sd3))], "years"]

  # Save the missing years and the failed years, if there is any
  if(length(sus_all) != 0 & length(missing_years) > 0) {
    check_folder <- path_sp(species)$check
    dir.create(check_folder)

    if(nrow(sus_all) > 0) {
      saveRDS(list(sus_all = sus_all,
                   sus_mean_sd1 = sus_mean_sd1,
                   sus_mean_sd2 = sus_mean_sd2,
                   sus_mean_sd3 = sus_mean_sd3,
                   sus_mean_only = sus_mean_only),
              paste0(check_folder, "/failed.rds"))
    }

    if(length(missing_years) > 0) {
      saveRDS(missing_years, paste0(check_folder, "/missing.rds"))
    }
  }
}
