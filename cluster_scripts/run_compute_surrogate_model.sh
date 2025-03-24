#!/bin/bash

#SBATCH --job-name=compute_surrogate_model
#SBATCH --output=logs/log-%j-%x.out
#SBATCH -c 10
#SBATCH --no-requeue   
#initialize module command
source /etc/profile

#load anaconda
module load anaconda/2023a
module load julia/1.10.1  

PROJECT_DIR="~/End2EndThermalImg.jl/Project.toml"
export JULIA_CONDAPKG_BACKEND="Null"
julia --project=${PROJECT_DIR} --threads 10 scripts/compute_surrogate_model.jl
