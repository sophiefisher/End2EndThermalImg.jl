#!/bin/bash

#SBATCH --job-name=test_incident_fields_speed
#SBATCH --output=logs/log-%j-%x.out
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

#initialize module command
source /etc/profile

#load anaconda
module load anaconda/Python-ML-2025a
module load julia/1.11.3

export JULIA_CONDAPKG_BACKEND="Null"
PROJECT_DIR="$HOME/End2EndThermalImg.jl"

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} scripts/test_incident_fields_speed.jl
