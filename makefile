# SDMs workflow
## Make species specific spatial objects necessaries
run_study_extent=analysis/make_study_extent.R
## Get species names
get_species=analysis/species.R
## Get species occurrences
get_occ=analysis/obs.R
## Analysis
run_sdm_inla=analysis/run_sdm_inla.R
run_sdm_gam=analysis/run_sdm_gam.R
## Make map(entire zone + qc)
run_maps_inla=analysis/run_maps_inla.R
run_maps_gam=analysis/run_maps_gam.R

# Folders
## Species specific spatial objects
study_extent=$(addsuffix /study_extent.rds, $(addprefix output/spatial/, $(species)))
mesh=$(addsuffix /mesh.rds, $(addprefix output/spatial/, $(species)))
rast=$(addsuffix /rast.gri, $(addprefix output/spatial/, $(species)))
filtered_obs=$(addsuffix /obs.rds, $(addprefix output/spatial/, $(species)))
## SDMs targets (multiple targets - one by species)
sdms_inla=$(addprefix output/models/inla/, $(species))
sdms_gam=$(addprefix output/models/gam/, $(species))
## Make map(entire zone + qc)
maps_inla=$(addprefix output/maps/inla/, $(species))
maps_gam=$(addprefix output/maps/gam/, $(species))
## Compute AUC
auc=$(addprefix output/auc/, $(species))
## Occurrences
occ=$(addsuffix .rds, $(addprefix occurrences/, $(species)))

# Arguments
res=10
proj='+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0'
year_start=1990
year_end=2020
window=5
time_buffer=21
spat_buffer=300
pedge=0.025
num_threads=$(cpu_task)
t1=0.05
t2=0.55


# Make map(entire zone + qc) and compute AUC
$(maps_inla): $(run_maps_inla) $(sdms_inla) 
	@Rscript $< $(proj) $@

maps_inla: $(maps_inla)

# Run SDMs for every species
$(sdms_inla): $(run_sdm_inla)
	@Rscript $< $@ $(year_start) $(year_end) $(window) $(num_threads)

models_inla: $(sdms_inla)

# Make map(entire zone + qc) and compute AUC
$(maps_gam): $(run_maps_gam)
	@Rscript $< $@

maps_gam: $(maps_gam)

# Run SDMs for every species
$(sdms_gam): $(run_sdm_gam)
	@Rscript $< $@ $(year_start) $(year_end) $(window)

models_gam: $(sdms_gam)

# Make spatial object necessary for the models
$(study_extent) $(mesh) $(rast) $(filtered_obs): $(run_study_extent)
	@Rscript $< $(species) $(spat_buffer) $(pedge)

spatial: $(study_extent) $(mesh) $(rast) $(filtered_obs)

# Make output folder
out_dir: 
	mkdir output
	mkdir output/check
	mkdir output/log
	mkdir output/maps
	mkdir output/models
	mkdir output/spatial
	mkdir occurrences

# Get species occurrences
$(occ): $(get_occ)
	@Rscript $< $@ $(year_start) $(year_end) $(window) $(buffer) $(proj)

occurrences: $(occ)

# Make species objects
species: $(get_species)
	@Rscript -e "source('$(get_species)')"

# install dependencies
install:
	Rscript -e 'if (!require(raster)) install.packages("raster", repos="https://cloud.r-project.org/");if (!require(terra)) install.packages("terra", repos="https://cloud.r-project.org/");if (!require(sp)) devtools::install_github("sp", repos="https://cloud.r-project.org/");if (!require(pROC)) install.packages("pROC", repos="https://cloud.r-project.org/");if (!require(stringr)) install.packages("stringr", repos="https://cloud.r-project.org/");if (!require(rebird)) install.packages("rebird", repos="https://cloud.r-project.org/");if (!require(rvest)) install.packages("rvest", repos="https://cloud.r-project.org/");if (!require(RPostgres)) install.packages("RPostgres", repos="https://cloud.r-project.org/");if (!require(atlasBE)) devtools::install_github("VincentBellavance/atlasBE");if (!require(INLA)) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE);if (!require(mapSpecies)) devtools::install_github("ReseauBiodiversiteQuebec/mapSpecies");if (!require(rgdal)) install.packages("rgdal", repos="https://cloud.r-project.org/")'

# Clean
clean:
	rm -r output/models/* output/log/* output/maps/*

.PHONY: install species occurrences out_dir spatial models_gam maps_gam clean
