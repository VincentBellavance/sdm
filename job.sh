#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH --mem=40GB
#SBATCH --time=04:00:00
#SBATCH --mail-user=vincent.bellavance@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=GAM_all_species
#SBATCH --array=1-196

module load gcc/9.3.0 gdal/3.2.3 r/4.1.2
export R_LIBS=~/R/x86_64-pc-linux-gnu-library/4.1
sp=$(cut -d' ' "-f${SLURM_ARRAY_TASK_ID}" <<<cat data/species_vect.txt)
make models_gam species=$sp
make maps_gam species=$sp
