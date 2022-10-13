# Get species names

#--------------------------------------------------------------
# Steps:
#   1. Import list of Quebec species (from QuébecOiseaux)
#   2. Keep Quebec breeding species
#   3. List family to remove (to remove marine species) and remove them from the list fo Quebec species
#   4. Remove northern nesting species
#   5. Validate taxonomy with Atlas taxonomy
#--------------------------------------------------------------

# Import list of Quebec species
breeding <- read.csv2("data/list_sp_qc.csv")
breeding <- breeding[breeding$status == "Nicheur", c("species", "family")]


# Remove family with marine species
marine_family <- c("anatidae",
                   "podicipedidae",
                   "stercorariidae",
                   "alcidae",
                   "laridae",
                   "gaviidae",
                   "hydrobatidae",
                   "sulidae",
                   "phalacrocoracidae")
keep_family <- c("rallidae",
                 "charadriidae",
                 "ardeidae")

breeding <- breeding[!breeding$family %in% marine_family, "species"]


# Remove northern nesting species
north_nest <- read.csv2("data/list_sp_qc_winter.csv")
north_nest <- north_nest[grep("\\(N\\)", north_nest$status), "species"]


# Validate names in breeding and north_west to compare with the same names
valid_breeding <- lapply(breeding, function(x) {
  tmp <- ratlas::get_taxa(scientific_name = x)
  return(unique(tmp[tmp$rank == "species", "valid_scientific_name"]))
})
valid_north_nest <- lapply(north_nest, function(x) {
  tmp <- ratlas::get_taxa(scientific_name = x)
  return(unique(tmp[tmp$rank == "species", "valid_scientific_name"]))
})
breeding <- as.vector(unlist(valid_breeding))
north_nest <- as.vector(unlist(valid_north_nest))


# Final list of species
final_list <- breeding[!breeding %in% north_nest]

sapply(final_list, function(x) {
  if(x == final_list[length(final_list)]) {
    cat(tolower(gsub(" ", "_", x)), 
        file = "data/species_vect.txt", 
        append = TRUE)
  } else {
    cat(paste0(tolower(gsub(" ", "_", x))," "), 
        file = "data/species_vect.txt", 
        append = TRUE)
  }
})

saveRDS(final_list, "data/species.rds")
