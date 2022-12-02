# SDMs workflow
## Get species names
get_species=analysis/species.R
## Get species occurrences
get_occ=analysis/obs.R
## Make species specific spatial objects necessaries
make_study_extent=analysis/make_study_extent.R
## Analysis
make_sdms=analysis/make_models.R
## Make map(entire zone + qc)
make_maps=analysis/make_maps.R
## Make binary maps
make_binary_maps=analysis/binarize_maps.R

## Occurrences
obs_dir=$(obs_dir)

# Output dirs
output_spatial=$(addsuffix $(zone), $(addprefix $(output_dir), /spatial/))
output_models=$(addsuffix $(zone), $(addprefix $(output_dir), /models/))
output_maps=$(addsuffix $(zone), $(addprefix $(output_dir), /maps/))
output_stack=$(addsuffix $(zone), $(addprefix $(output_dir), /stack/))
output_log=$(addsuffix $(zone), $(addprefix $(output_dir), /log/))
output_out=$(addsuffix $(zone), $(addprefix $(output_dir), /out/))
output_bdi=$(addsuffix $(zone), $(addprefix $(output_dir), /bdi/))

## Species specific spatial objects
study_extent=$(addsuffix /study_extent.rds, $(addprefix $(output_spatial), $(species)))
mesh=$(addsuffix /mesh.rds, $(addprefix $(output_spatial), $(species)))
rast=$(addsuffix /rast.gri, $(addprefix $(output_spatial), $(species)))
filtered_obs=$(addsuffix /obs.rds, $(addprefix $(output_spatial), $(species)))
## SDMs targets (multiple targets - one by species)
sdms=$(addprefix $(output_models), $(species))
## Make maps
maps=$(addprefix $(addprefix $(output_maps), $(species)), /maps_pocc.gri)
binary_maps_range=$(addprefix $(addprefix $(output_maps), $(species)), /maps_range.gri)
binary_maps_occ=$(addprefix $(addprefix $(output_maps), $(species)), /maps_occ.gri)

# Arguments
res=10
proj='+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=62 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0'
year_start=1990
year_end=2020
window=5
time_buffer=28
spat_buffer=300
pedge=0.03
num_threads=$(cpu_task)
sensitivity=0.99

# Make binary maps
$(binary_maps_range) $(binary_maps_occ): $(make_binary_maps)
	@Rscript $< $(species) $(zone) $(sensitivity) $(output_dir)

binary_maps: $(binary_maps_range) $(binary_maps_occ)

# Make map(entire zone + qc) and compute AUC
$(maps): $(make_maps) 
	@Rscript $< $(species) $(zone) $(output_dir)

maps: $(maps)

# Run SDMs for every species
$(sdms): $(make_sdms)
	@Rscript $< $(species) $(year_start) $(year_end) $(window) $(num_threads) $(zone) $(output_dir)

models: $(sdms)

# Make spatial object necessary for the models
$(study_extent) $(mesh) $(rast) $(filtered_obs): $(make_study_extent)
	@Rscript $< $(species) $(zone) $(spat_buffer) $(pedge) $(obs_dir) $(output_dir)

spatial: $(study_extent) $(mesh) $(rast) $(filtered_obs)

# Make output dir
out_dir: 
	if ! [ -d $(output_dir) ]; then \
	  mkdir $(output_dir); \
	  mkdir $(addprefix $(output_dir), /log); \
	  mkdir $(addprefix $(output_dir), /out); \
	  mkdir $(addprefix $(output_dir), /models); \
	  mkdir $(addprefix $(output_dir), /stack); \
	  mkdir $(addprefix $(output_dir), /maps); \
	  mkdir $(addprefix $(output_dir), /spatial); \
	  mkdir $(addprefix $(output_dir), /bdi); \
	fi
	mkdir $(output_spatial)
	mkdir $(output_models)
	mkdir $(output_maps)
	mkdir $(output_stack)
	mkdir $(output_log)
	mkdir $(output_out)
	mkdir $(output_bdi)

# Get species occurrences
$(obs_dir): $(get_obs)
	@Rscript $< $(year_start) $(year_end) $(time_buffer) $(proj) $(obs_dir)

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
