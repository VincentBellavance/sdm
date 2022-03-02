# Get species names (accepted and synonym)

#--------------------------------------------------------------
# Steps:
#--------------------------------------------------------------

# Query species for three different sources
gbif <- ratlas::list_bird_taxa(rank = "species", source_name = "GBIF Backbone Taxonomy") |>
          (\(.) as.data.frame(.))()
col <- ratlas::list_bird_taxa(rank = "species", source_name = "Catalogue of Life") |>
          (\(.) as.data.frame(.))()

# Remove species that are not breeding
breeding <- read.csv2("data/list_sp_qc.csv") |>
              {\(x) x[x$status == "Nicheur", "species"]}()
col <- col[col$valid_srid %in% col[col$scientific_name %in% breeding, "valid_srid"],]
gbif <- gbif[gbif$valid_srid %in% gbif[gbif$scientific_name %in% breeding, "valid_srid"],]

# Remove marine species (for now)
marine <- read.csv2("data/list_marine_sp.csv")
gbif <- gbif[!gbif$valid_srid %in% gbif[gbif$scientific_name %in% marine$species, "valid_srid"],]
col <- col[!col$valid_srid %in% col[col$scientific_name %in% marine$species, "valid_srid"],]

# Remove species that do not have enough observations
avg_limit <- 30
tot_limit <- 300
col <- col[col$avg_yearly_count >= avg_limit & col$total_count >= tot_limit, ]
gbif <- gbif[gbif$avg_yearly_count >= avg_limit & gbif$total_count >= tot_limit, ]


list_ref <- lapply(unique(col$valid_srid), function(x) {
  
  accepted <- col[col$valid_srid == x &
                  col$valid, "scientific_name"]
  
  tmp <- list(
    accepted = accepted
  )

  return(tmp)

})


# Synonym in GBIF that are accepted in COL
accepted_col <- unlist(lapply(list_ref, "[[", 1))

for(i in unique(gbif[!gbif$valid, "scientific_name"])) {
  if(i %in% accepted_col) {
    valid_srid <- gbif[gbif$scientific_name == i, "valid_srid"]
    syn_gbif <- gbif[gbif$valid_srid == valid_srid & gbif$valid, "scientific_name"]
    list_ref[[which(accepted_col %in% i)]][["synonym"]] <- unique(c(syn_gbif, list_ref[[which(accepted_col %in% i)]][["synonym"]]))
  }
}

# Accepted GBIF that are not in COL
add_gbif <- gbif[!gbif$scientific_name %in% unlist(list_ref) & gbif$valid, 
                        c("scientific_name", "valid_srid", "valid", "source_record_id")]

gbif_list <- lapply(unique(add_gbif$valid_srid), function(x) {
  
  accepted <- gbif[gbif$valid_srid == x &
                  gbif$valid, "scientific_name"]
  synonym <- gbif[gbif$valid_srid == x & 
                 !gbif$valid, "scientific_name"]
  
  tmp <- list(
    accepted = accepted
  )

  if(length(synonym) > 0) {
    tmp[["synonym"]] <- synonym
  }

  return(tmp)

})

final_list <- c(list_ref, gbif_list)

#final_list <- list(list(accepted = "Grus canadensis", synonym = "Antigone canadensis"), list(accepted = "Setophaga petechia"))

lapply(final_list, function(x) {
  cat(paste0(tolower(gsub(" ", "_", x$accepted))," "), file = "data/species_vect.txt", append = TRUE)
})

saveRDS(final_list, "data/species.rds")
