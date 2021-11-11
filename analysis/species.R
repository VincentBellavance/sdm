# Get species names (accepted and synonym)

#--------------------------------------------------------------
# Steps:
#--------------------------------------------------------------

sp <- list(
        list(
          accepted = "Antigone canadensis"
        ),
        list(
          accepted = "Cardinalis cardinalis"
        )
      )

saveRDS(sp, "data/species.rds")

lapply(sp, function(x) {
  cat(paste0(tolower(gsub(" ", "_", x$accepted)),"  "), file = "data/species_vect.txt", append = TRUE)
})
