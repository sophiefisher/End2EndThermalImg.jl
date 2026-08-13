#!/bin/bash

#SBATCH --job-name=test_multinode_setup
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

julia --project=${PROJECT_DIR} scripts/test_multinode_setup.jl
