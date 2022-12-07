# Permute the pair species-occurrence in the long format table of observations for 25 random species to make a null model of the Biodiversity Distribution Index

#--- SETUP ---#
library(raster)
args = commandArgs(trailingOnly=TRUE)
i = as.integer(args[1])
subset_sp <- read.table("data/species_null.txt", sep = " ")

# Import observations table of all species from the subset
obs_sp <- lapply(subset_sp, function(x) {
  tmp <- readRDS(paste0("data/occurrences/",x,".rds"))
  tmp$taxa_scientific_name <- stringr::str_to_sentence(gsub("_", " ", x))
  return(as.data.frame(tmp))
})

# Build one big table of all observations in dataframe class with coords
obs_df <- do.call(rbind, obs_sp)

# Make the permutations
cols_perm <- c("taxa_scientific_name", "occurrence")
sampl <- sample(1:nrow(obs_df))
obs_df[,cols_perm] <- obs_df[sampl,cols_perm]

# Make sp objects for every species and save object in occurrence_null directory
perm_obs_sp <- lapply(subset_sp, function(x) {
  tmp <- obs_df[obs_df$taxa_scientific_name == stringr::str_to_sentence(gsub("_", " ", x)),]
  tmp_sp <- sp::SpatialPointsDataFrame(tmp[,c("lon", "lat")], 
                             tmp[,-which(colnames(tmp) %in% c("lon", "lat"))])
  saveRDS(tmp_sp, paste0("data/occurrences_null/iter_",i,"/",x,".rds"))
  return(tmp_sp)
})

