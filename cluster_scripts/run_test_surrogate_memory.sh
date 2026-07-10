#!/bin/bash

#SBATCH --job-name=test_surrogate_memory
#SBATCH --output=logs/log-%j-%x.out
#SBATCH --cpus-per-task=8

#initialize module command
source /etc/profile

#load anaconda
module load anaconda/Python-ML-2025a
module load julia/1.11.3


PROJECT_DIR="$HOME/End2EndThermalImg.jl"
export JULIA_CONDAPKG_BACKEND="Null"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1      # add this too — some numpy builds use MKL rather than OpenBLAS
/usr/bin/time -v julia --project=${PROJECT_DIR} scripts/test_surrogate_memory.jl
