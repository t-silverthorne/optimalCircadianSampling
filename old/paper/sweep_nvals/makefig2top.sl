#!/bin/bash -l
#SBATCH --job-name=makefig2_top_p3
#SBATCH --account=def-stinch   
#SBATCH --time=8:00:00          # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1      # currently set for parallel commands
#SBATCH --mem-per-cpu=12000      # currently set for parallel
#SBATCH --mail-user=turner.silverthorne@utoronto.ca 
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2022b.2
## matlab -nodisplay -r "clear; partname=1; Nmeasvals=8:15; sweep_makefig2_top"
## matlab -nodisplay -r "clear; partname=2; Nmeasvals=16:23;sweep_makefig2_top"
matlab -nodisplay -r "clear; partname=3; Nmeasvals=24:31;sweep_makefig2_top"
