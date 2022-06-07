# SDMs workflow
## Make species specific spatial objects necessaries
make_study_extent=analysis/make_study_extent.R
## Get species names
get_species=analysis/species.R
## Get species occurrences
get_occ=analysis/occ.R
## Analysis
make_sdms=analysis/make_models.R
## Make map(entire zone + qc)
make_maps=analysis/make_maps.R
## Run check
check_models=analysis/check_models.R

# Data folders
data_folder=$(data_folder)
## Occurrences
obs_folder=$(addprefix $(data_folder), /occurrences/)
obs=$(addsuffix .rds, $(addprefix $(obs_folder), $(species)))

# Output folders
output=$(output_folder)
output_spatial=	$(addprefix $(output_folder), /spatial)
output_models=$(addprefix $(output_folder), /models)
output_maps=$(addprefix $(output_folder), /maps)
output_stack=$(addprefix $(output_folder), /stack)

## Species specific spatial objects
study_extent=$(addsuffix /study_extent.rds, $(addprefix $(output_spatial), $(species)))
mesh=$(addsuffix /mesh.rds, $(addprefix $(output_spatial), $(species)))
rast=$(addsuffix /rast.gri, $(addprefix $(output_spatial), $(species)))
filtered_obs=$(addsuffix /obs.rds, $(addprefix $(output_spatial), $(species)))
## SDMs targets (multiple targets - one by species)
sdms=$(addprefix $(output_models), $(species))
## Make map(entire zone + qc)
maps=$(addprefix $(output_maps), $(species))
## Check
checks=$(addprefix $(output_check), $(species))

# Arguments
res=10
proj='+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0'
year_start=1990
year_end=2020
window=5
time_buffer=21
spat_buffer=300
pedge=0.03
num_threads=$(cpu_task)
t1=0.05
t2=0.55

# Run checks on models
$(checks): $(check_models) $(maps)
	@Rscript $< $(species)

checks: $(checks)

# Make map(entire zone + qc) and compute AUC
$(maps): $(make_maps) $(sdms) 
	@Rscript $< $(species) $(proj)

maps: $(maps)

# Run SDMs for every species
$(sdms): $(make_sdms)
	@Rscript $< $(species) $(year_start) $(year_end) $(window) $(num_threads)

models: $(sdms)

# Make spatial object necessary for the models
$(study_extent) $(mesh) $(rast) $(filtered_obs): $(make_study_extent)
	@Rscript $< $(species) $(spat_buffer) $(pedge)

spatial: $(study_extent) $(mesh) $(rast) $(filtered_obs)

# Make output folder
out_dir: 
	mkdir $(output_folder)
	mkdir $(addprefix $(output_folder), /checks)
	mkdir $(addprefix $(output_folder), /log)
	mkdir $(output_spatial)
	mkdir $(output_models)
	mkdir $(output_maps)
	mkdir $(output_stack)
	mkdir $(obs_folder)

# Get species occurrences
$(obs_folder): $(get_obs)
	@Rscript $< $@ $(year_start) $(year_end) $(window) $(time_buffer) $(proj)

occurrences: $(obs)

# Make species objects
species: $(get_species)
	@Rscript -e "source('$(get_species)')"

# install dependencies
install:
	Rscript -e 'if (!require(raster)) install.packages("raster", repos="https://cloud.r-project.org/");if (!require(terra)) install.packages("terra", repos="https://cloud.r-project.org/");if (!require(sp)) devtools::install_github("sp", repos="https://cloud.r-project.org/");if (!require(pROC)) install.packages("pROC", repos="https://cloud.r-project.org/");if (!require(stringr)) install.packages("stringr", repos="https://cloud.r-project.org/");if (!require(rebird)) install.packages("rebird", repos="https://cloud.r-project.org/");if (!require(rvest)) install.packages("rvest", repos="https://cloud.r-project.org/");if (!require(RPostgres)) install.packages("RPostgres", repos="https://cloud.r-project.org/");if (!require(atlasBE)) devtools::install_github("VincentBellavance/atlasBE");if (!require(INLA)) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE);if (!require(mapSpecies)) devtools::install_github("ReseauBiodiversiteQuebec/mapSpecies");if (!require(rgdal)) install.packages("rgdal", repos="https://cloud.r-project.org/")'

# Clean
clean:
	rm -r output/models/* output/log/* output/maps/*

.PHONY: models maps install species occurrences out_dir spatial clean
