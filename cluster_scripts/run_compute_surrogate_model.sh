#!/bin/bash

#SBATCH --job-name=compute_surrogate_model
#SBATCH --output=logs/log-%j-%x.out
#SBATCH -n 15
#SBATCH -c 24
#initialize module command
source /etc/profile

#load anaconda
module load anaconda/2023a
module load julia/1.10.1  

echo "Number of tasks: $SLURM_NTASKS"
echo "Cores per task: $SLURM_CPUS_PER_TASK"
TOTAL_CORES=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))
echo "Total logical cores used: $TOTAL_CORES"

PROJECT_DIR="~/End2EndThermalImg.jl/Project.toml"
export JULIA_CONDAPKG_BACKEND="Null"
julia --project=${PROJECT_DIR} -p ${SLURM_NTASKS} scripts/compute_surrogate_model.jl
