#!/bin/bash -l
#SBATCH --job-name=optimization_test
#SBATCH --account=def-stinch   
#SBATCH --time=7:00:00          # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30      # currently set for parallel commands
#SBATCH --mem-per-cpu=4500      # currently set for parallel
#SBATCH --mail-user=turner.silverthorne@utoronto.ca 
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2022b.2
matlab -nodisplay -r "clear; optimization_test"

