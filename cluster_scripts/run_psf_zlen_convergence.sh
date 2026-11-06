#!/bin/bash

#SBATCH --job-name=psf_zlen_convergence
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

# precompile once, serially, before the driver and workers below start concurrently -- see
# run_endtoend_opt.sh for why (shared networked ~/.julia depot, precompile-lock contention).
julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} -e 'import Pkg; Pkg.instantiate(); Pkg.precompile(io=devnull)'

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} scripts/psf_zlen_convergence.jl

# Copy the complete job log into its own convergence_studies/psf_zlen_<datetime>/ run folder,
# alongside results.jld2 -- same convenience/rationale as run_endtoend_opt.sh's copy into
# optdata/opt_<datetime>/. No-ops if the job never reached the final save (e.g. a crash mid-sweep).
LOG_FILE="logs/log-${SLURM_JOB_ID}-${SLURM_JOB_NAME}.out"
RUN_FOLDER=$(grep -m 1 "^run_folder" "$LOG_FILE" | sed 's/^run_folder[[:space:]]*=[[:space:]]*//')
if [ -n "$RUN_FOLDER" ] && [ -d "$RUN_FOLDER" ]; then
    cp "$LOG_FILE" "$RUN_FOLDER/log.out"
fi
