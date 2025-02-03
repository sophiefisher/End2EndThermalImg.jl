#!/bin/bash

#SBATCH --job-name=time_compute_surrogate_model
#SBATCH --output=logs/log-%j-%x.out
#SBATCH -c 40
#initialize module command
source /etc/profile

#load anaconda
module load anaconda/2023a
module load julia/1.10.1  

PROJECT_DIR="~/End2EndThermalImg.jl/Project.toml"
julia --project=${PROJECT_DIR} --threads 40 scripts/time_compute_surrogate_model.jl
