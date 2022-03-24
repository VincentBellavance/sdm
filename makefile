# SDMs workflow
## Make spatial object necessary for the models
run_spat_data=analysis/run_spat_data.R
## Get species names
get_species=analysis/species.R
## Get species occurrences
get_occ=analysis/obs.R
## Analysis
run_sdm=analysis/run_sdm_inla.R
run_sdm_gam=analysis/run_sdm_gam.R
## Make map(entire zone + qc)
run_maps=analysis/run_maps.R
run_maps_gam=analysis/run_maps_gam.R

# Folders
## Make spatial object necessary for the models
spacePoly=data/spacePoly.rds
qc=data/qc_spacePoly.rds
explana=data/explana.rds
mesh=data/mesh.rds
rast=data/rast*
## SDMs targets (multiple targets - one by species)
sdms=$(addprefix output/models/inla/, $(species))
## Make map(entire zone + qc)
maps=$(addprefix output/maps/inla/, $(species))
## Compute AUC
auc=$(addprefix output/auc/, $(species))
## Occurrences
occ=$(addsuffix .rds, $(addprefix occurrences/, $(species)))

# Arguments
res=10
proj='+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0'
year_start=2013
year_end=2017
window=5
buffer=21
num_threads=$(cpu_task)
t1=0.05
t2=0.55


# Make map(entire zone + qc) and compute AUC
$(maps): $(run_maps) $(sdms) 
	@Rscript $< $(proj) $@

maps: $(maps)

# Run SDMs for every species
$(sdms): $(run_sdm)
	@Rscript $< $@ $(year_start) $(year_end) $(window) $(num_threads)

models: $(sdms)

# Make map(entire zone + qc) and compute AUC
#$(maps): $(run_maps_gam)
#	@Rscript $< $@

#maps_gam: $(maps)

# Run SDMs for every species
#$(sdms): $(run_sdm_gam)
#	@Rscript $< $@ $(year_start) $(year_end) $(window)

#models_gam: $(sdms)

# Make spatial object necessary for the models
$(spacePoly) $(qc) $(explana) $(mesh) $(rast): $(run_spat_data)
	@Rscript $< $(res) $(proj)

spatial: $(spacePoly) $(qc) $(explana) $(mesh) $(rast)

# Get species occurrences
$(occ): $(get_occ)
	@Rscript $< $@ $(year_start) $(year_end) $(window) $(buffer) $(proj)

occurrences: $(occ)

# Make output folder
out_dir: 
	mkdir output
	mkdir output/auc
	mkdir output/log
	mkdir output/maps
	mkdir output/models
	mkdir output/models/gam
	mkdir output/models/inla
	mkdir occurrences

# Make species objects
species: $(get_species)
	@Rscript -e "source('$(get_species)')"

# install dependencies
install:
	Rscript -e 'if (!require(raster)) install.packages("raster", repos="https://cloud.r-project.org/");if (!require(terra)) install.packages("terra", repos="https://cloud.r-project.org/");if (!require(sp)) devtools::install_github("sp", repos="https://cloud.r-project.org/");if (!require(pROC)) install.packages("pROC", repos="https://cloud.r-project.org/");if (!require(stringr)) install.packages("stringr", repos="https://cloud.r-project.org/");if (!require(rebird)) install.packages("rebird", repos="https://cloud.r-project.org/");if (!require(rvest)) install.packages("rvest", repos="https://cloud.r-project.org/");if (!require(RPostgres)) install.packages("RPostgres", repos="https://cloud.r-project.org/");if (!require(atlasBE)) devtools::install_github("VincentBellavance/atlasBE");if (!require(INLA)) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE);if (!require(mapSpecies)) devtools::install_github("ReseauBiodiversiteQuebec/mapSpecies");if (!require(rgdal)) install.packages("rgdal", repos="https://cloud.r-project.org/")'

# Clean
clean:
	rm -r output/models/* output/log/*

.PHONY: install species occurrences out_dir spatial models_gam maps_gam clean
