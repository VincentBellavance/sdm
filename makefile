# SDMs workflow
## Make spatial object necessary for the models
run_spat_data=analysis/run_spat_data.R
## Get species names
get_species=analysis/species.R
## Analysis
run_sdm=analysis/run_sdm.R
## Make map(entire zone + qc)
run_maps=analysis/run_maps.R

# Folders
## Make spatial object necessary for the models
spacePoly=data/spacePoly.rds
qc=data/qc_spacePoly.rds
explana=data/explana.rds
mesh=data/mesh.rds
rast=data/rast*
## SDMs targets (multiple targets - one by species)
FILE:=data/species_vect.txt
species_vect:=$(file < $(FILE))
sdms=$(addprefix output/models/, $(species_vect))
## Make map(entire zone + qc)
maps=$(addprefix output/maps/, $(species_vect))
## Compute AUC
auc=$(addprefix output/auc/, $(species_vect))

# Arguments
res=2
proj='+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0'
year_start=1990
year_end=2020
window=5
buffer=21
num_threads=30
t1=0.05
t2=0.55


## Make map(entire zone + qc) and compute AUC
$(maps) $(auc): $(run_maps) $(sdms) 
	@Rscript $< $(t1) $(t2) $(proj) $@

# Run SDMs for every species
$(sdms): $(run_sdm) $(spacePoly) $(rast) $(mesh) $(explana)
	@Rscript $< $@ $(year_start) $(year_end) $(window) $(buffer) $(proj) $(num_threads)

# Make spatial object necessary for the models
$(spacePoly) $(qc) $(explana) $(mesh) $(rast): $(run_spat_data)
	@Rscript --vanilla $(run_spat_data) $(res) $(proj)

# Declare all
all: $(maps)

# Make species objects
species: $(get_species)
	@Rscript -e "source('$(get_species)')"

# install dependencies
install:
	Rscript -e 'if (!require(raster)) install.packages("raster");if (!require(terra)) install.packages("terra");if (!require(sp)) devtools::install_github("sp");if (!require(pROC)) install.packages("pROC");if (!require(stringr)) install.packages("stringr");if (!require(rebird)) install.packages("rebird");if (!require(rvest)) install.packages("rvest");if (!require(RPostgres)) install.packages("RPostgres");if (!require(atlasBE)) devtools::install_github("VincentBellavance/atlasBE");if (!require(mapSpecies)) devtools::install_github("ReseauBiodiversiteQuebec/mapSpecies");if (!require(rgdal)) install.packages("rgdal")'

.PHONY: install species
