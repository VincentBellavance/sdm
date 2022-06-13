#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH --mem=20GB
#SBATCH --time=02:00:00
#SBATCH --mail-user=vincent.bellavance@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=run_sdm
#SBATCH --output=/home/belv1601/scratch/output/out/%a.out
#SBATCH --array=1

module use /home/belv1601/.local/easybuild/modules/2020/avx2/Compiler/gcc9/
module load StdEnv/2020  gcc/9.3.0 r-inla/21.05.02 geos/3.9.1 gdal/3.0.4 proj/7.0.1 udunits
sp=$(cut -d' ' "-f${SLURM_ARRAY_TASK_ID}" <<<cat data/species_vect.txt)
make spatial species=$sp cpu_task=1 output_folder=/home/belv1601/scratch/output data_folder=/home/belv1601/projects/def-dgravel/belv1601/sdm/data
make models species=$sp cpu_task=1
make maps species=$sp cpu_task=1
make checks species=$sp
