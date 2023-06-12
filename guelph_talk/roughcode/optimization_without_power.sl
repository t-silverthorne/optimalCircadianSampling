#!/bin/bash -l
#SBATCH --job-name=optimize_no_power
#SBATCH --account=def-stinch   
#SBATCH --time=1:00:00          # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25      # currently set for parallel commands
#SBATCH --mem-per-cpu=2500      # currently set for parallel
#SBATCH --mail-user=turner.silverthorne@utoronto.ca 
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2022b.2
matlab -nodisplay -r "clear; p.Nmeas=8;p.freq=3.8;p.Amp=1.5;optimize_without_power"

