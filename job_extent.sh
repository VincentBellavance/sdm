#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH --time=1-12:00:00
#SBATCH --mail-user=vincent.bellavance@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name=study_extent
#SBATCH --output=study_extent.out

module use /home/belv1601/.local/easybuild/modules/2020/avx2/Compiler/gcc9/
module load StdEnv/2020  gcc/9.3.0 r-inla/21.05.02 geos/3.9.1 gdal/3.0.4 proj/7.0.1 udunits
#sp=$(cut -d' ' "-f${SLURM_ARRAY_TASK_ID}" <<<cat data/species_vect.txt)
make spatial species="pinicola_enucleator" cpu_task=20
