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

# precompile once, serially, before the driver and workers below start concurrently --
# they share the same networked depot (~/.julia), so without this they each race to
# precompile the same packages at once and end up lock-waiting on each other
julia --project=${PROJECT_DIR} -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} scripts/test_incident_fields_speed.jl
