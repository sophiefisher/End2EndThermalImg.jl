#!/bin/bash

#SBATCH --job-name=test_surrogate_memory
#SBATCH --output=logs/log-%j-%x.out

#initialize module command
source /etc/profile

#load anaconda
module load anaconda/Python-ML-2025a
module load julia/1.11.3


PROJECT_DIR="~/End2EndThermalImg.jl/Project.toml"
export JULIA_CONDAPKG_BACKEND="Null"
julia --project=${PROJECT_DIR} -p ${SLURM_NTASKS} scripts/test_surrogate_memory.jl
