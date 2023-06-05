#!/bin/bash -l
#SBATCH --job-name=testperms
#SBATCH --account=def-stinch   
#SBATCH --time=0:10:00          # adjust this to match the walltime of your job

#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12      # adjust this if you are using parallel commands
#SBATCH --mem=30000             # adjust this according to the memory requirement per node you need

#SBATCH --mail-user=turner.silverthorne@utoronto.ca 
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2022b.2
matlab -nodisplay -r "clear; d=10;N=1e4;sandbox;"
matlab -nodisplay -r "clear; d=10;N=1e5;sandbox;"
matlab -nodisplay -r "clear; d=10;N=1e6;sandbox;"
matlab -nodisplay -r "clear; d=10;N=1e7;sandbox;"
