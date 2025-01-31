#!/bin/bash

CORES=40  # Define the number of cores
SCRIPT_NAME="time_compute_surrogate_model"

#SBATCH --job-name=$SCRIPT_NAME
#SBATCH --output=logs/log-%j-%x.out
#SBATCH -c $CORES
#initialize module command
source /etc/profile

#load anaconda
module load anaconda/2023b
module load julia/1.10.1  

julia --threads $CORES ${SCRIPT_NAME}.jl
