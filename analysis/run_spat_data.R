#!/usr/bin/env Rscript
# Make spatial objects that will be used in the models

# Import mapSpecies
suppressMessages(library(mapSpecies))

#------------------------------------------------
# 1. Make spatial polygon of the modeled area
# 2. Make spatial polygon of Quebec
# 3. Make the mesh
# 4. Make raster that will be used for prediction
# 5. Make explana Mesh
#------------------------------------------------

# Set variables arguments
args = commandArgs(trailingOnly=TRUE)
res <- as.integer(args[1]) # 10km
proj <- args[2]

# Import functions
source("R/prep.R")

# Make spatial polygon of the modeled area
q <- prep_spat_poly(proj)
print("q done")
# Make spatial polygon of Quebec
qc <- prep_qc_spat_poly(proj)
print("qc done")
# Make the mesh
mesh <- prep_mesh(q, proj)
print("mesh done")
# Make raster that will be used for prediction
rast <- prep_rast_pred(q, proj, res)
print("rast done")
# Make explana Mesh
explana <- prep_explana(q, mesh, rast)
print("explana done")
