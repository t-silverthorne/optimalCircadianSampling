#!/bin/bash -l
#SBATCH --job-name=scalarization
#SBATCH --account=def-stinch   
#SBATCH --time=8:00:00          # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30      # currently set for parallel commands
#SBATCH --mem-per-cpu=1600      # currently set for parallel
#SBATCH --mail-user=turner.silverthorne@utoronto.ca 
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2022b.2
matlab -nodisplay -r "run_sandbox"
