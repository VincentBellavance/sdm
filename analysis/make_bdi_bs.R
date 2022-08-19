# Calculate the BDI with the generated maps by resampling the posteriors of the models' fitted values

# Import index number
args <- commandArgs(trailingOnly=TRUE)
i <- as.integer(args[1])

# Import rbdi library
library(rbdi)
library(raster)

source("R/path.R")
# Get vector of species name
output_folder <- path()$maps
species <- list.dirs(output_folder, recursive = FALSE, full.names = FALSE)

# Import species SDMs
sdms <- lapply(species, function(x) {
  maps_file <- paste0(path_sp(x)$maps, "/maps_bs.gri")
  if(file.exists(maps_file)){
    print("ok")
    tmp <- raster::stack(maps_file, 
                         bands = seq(from = 0, by = 1000, length.out = 27)+i)
    names(tmp) <- gsub("_[0-9]+$", "", names(tmp))
    return(tmp)
  }
})
species <- species[!sapply(sdms, is.null) & sapply(sdms, function(x) length(names(x))) == 27]
sdms <- sdms[!sapply(sdms, is.null)]
sdms <- sdms[sapply(sdms, function(x) length(names(x))) == 27]
names(sdms) <- species

# Calculate BDI
bdi <- calc_bdi(sdms, first_year = 1992)

# Save bdi object
write.table(bdi, 
            paste0("~/scratch/output/bdi/bdi",i,".txt"), 
            sep = ",",
            dec = ".")
