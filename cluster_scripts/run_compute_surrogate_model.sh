#!/bin/bash

#SBATCH --job-name=compute_surrogate_model
#SBATCH --output=logs/log-%j-%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

#initialize module command
source /etc/profile

#load anaconda
module load anaconda/Python-ML-2025a
module load julia/1.11.3

export OPENBLAS_NUM_THREADS=3
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export JULIA_CONDAPKG_BACKEND="Null"
PROJECT_DIR="$HOME/End2EndThermalImg.jl"

julia --project=${PROJECT_DIR} -p 15 scripts/compute_surrogate_model.jl