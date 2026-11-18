#!/bin/bash

#SBATCH --job-name=checkpoint_report
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

# Same node/cpu configuration as run_endtoend_opt.sh -- scripts/checkpoint_report.jl pays the same
# per-reconstruction forward-model cost (PSF formation across the freq x PSF_zlen grid) as any one
# outer-loop iteration of the original job, just for 3 reconstructions instead of ~1000+, so
# downsizing the allocation would only make each reconstruction slower for no real savings in
# node-hours (see the planning discussion that led to this script).
E2E_RUN_FOLDER="${E2E_RUN_FOLDER:?E2E_RUN_FOLDER must be set to the optdata/opt_<datetime> run folder to report on, e.g. E2E_RUN_FOLDER=optdata/opt_2026-11-13T13:07:16.666 sbatch cluster_scripts/run_checkpoint_report.sh}"

# precompile once, serially -- see run_endtoend_opt.sh for why (shared networked depot lock
# contention across the driver and workers below).
julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} -e 'import Pkg; Pkg.instantiate(); Pkg.precompile(io=devnull)'

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} scripts/checkpoint_report.jl

# Copy the complete job log into the run folder, same as run_endtoend_opt.sh's own log.out --
# unlike that script, E2E_RUN_FOLDER is already known here (not discovered via grep), since it's
# the input to this job rather than something scripts/checkpoint_report.jl generates fresh.
LOG_FILE="logs/log-${SLURM_JOB_ID}-${SLURM_JOB_NAME}.out"
if [ -d "$E2E_RUN_FOLDER" ]; then
    cp "$LOG_FILE" "$E2E_RUN_FOLDER/log.out"
fi
