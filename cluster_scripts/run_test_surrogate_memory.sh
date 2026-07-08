#!/bin/bash

#SBATCH --job-name=test_surrogate_memory
#SBATCH --output=logs/log-%j-%x.out

#initialize module command
source /etc/profile

#load anaconda
module load anaconda/Python-ML-2025a
module load julia/1.11.3


PROJECT_DIR="$HOME/End2EndThermalImg.jl"
export JULIA_CONDAPKG_BACKEND="Null"
/usr/bin/time -v julia --project=${PROJECT_DIR} scripts/test_surrogate_memory.jl
