#!/bin/bash

#SBATCH --job-name=add_psfs_initial
#SBATCH --output=logs/log-%j-%x.out
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

#initialize module command
source /etc/profile

#load anaconda
module load anaconda/Python-ML-2025a
module load julia/1.11.3

export JULIA_CONDAPKG_BACKEND="Null"
PROJECT_DIR="$HOME/End2EndThermalImg.jl"

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} -e 'import Pkg; Pkg.instantiate(); Pkg.precompile(io=devnull)'

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} scripts/add_psfs_initial.jl
