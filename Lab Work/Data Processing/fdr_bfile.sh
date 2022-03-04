#!/bin/bash

#SBATCH --job-name=fdr_check   # job name
#SBATCH --output=fdr_check.txt  # output log file
#SBATCH --error=fdr_check.err  # error file
#SBATCH --time=01:00:00  # 1hr of wall time
#SBATCH --nodes=1        # 1 GPU node
#SBATCH --partition=gpu2 # GPU2 partition
#SBATCH --ntasks=1       # 1 CPU core to drive GPU
#SBATCH --gres=gpu:1     # Request 1 GPU

# Load all required modules below. 
module load python/anaconda-2020.02

# Add lines here to run your GPU-based computations.
python3 FDR_value_check.py