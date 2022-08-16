# Calculate the BDI with the generated maps by resampling the posteriors of the models' fitted values

# Import index number
args <- commandArgs(trailingOnly=TRUE)
i <- args[1]

# Import rbdi library
library(rbdi)

# Get vector of species name
output_folder <- path()$maps
species <- list.dir(output_folder, recursive = FALSE, full.names = FALSE)

# Import species SDMs
sdms <- lapply(species, function(x) {
  maps_file <- paste0(path_sp(x), "maps_bs.gri")
  if(file.exists(maps_file)){
      raster::stack(maps_file, 
                    bands = seq(from = 0, by = 1000, length.out = 27)+i)
  }
})
names(sdms) <- species

# Calculate BDI
bdi <- rbdi::calc_bdi(sdms, first_year = 1992)

# Save bdi object
write.table(bdi, 
            paste0("~/scratch/output/bdi/bdi",i,".txt"), 
            sep = ",",
            dec = ".")
