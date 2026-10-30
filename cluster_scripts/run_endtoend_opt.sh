#!/bin/bash

#SBATCH --job-name=endtoend_opt
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

# precompile once, serially, before the driver and workers below start concurrently -- they share
# the same networked depot (~/.julia), so without this they each race to precompile the same
# packages at once and end up lock-waiting on each other. --threads must match the run below,
# since some packages' caches are invalidated when reloaded under a different thread count, which
# otherwise forces a second precompile pass.
# io=devnull: this step still fills log-%j-%x.out with many lines of Pkg's own
# compilecache-contention warnings (other jobs on the cluster share this same networked depot) --
# harmless (a real failure here still throws and halts the script), but drowns out the actual run.
# Pkg.precompile's `io` kwarg (default stderr) controls exactly this display.
julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} -e 'import Pkg; Pkg.instantiate(); Pkg.precompile(io=devnull)'

julia --project=${PROJECT_DIR} --threads=${SLURM_CPUS_PER_TASK} scripts/endtoend_opt.jl

# Copy the complete job log into its own optdata/opt_<datetime>/ run folder, alongside
# results.jld2/plots/etc., so the run's own log doesn't have to be hunted down separately by job ID
# in logs/. Done here, after the julia run above has fully finished and the log is complete -- not
# from inside the Julia script itself, since anything written to the log after a self-copy (the
# copy's own confirmation, trailing SLURM output) would be missing from that copy. Left in logs/
# either way (not moved): a lookup-by-job-ID convenience that costs negligible disk space, and a
# bug in this grep then just means "no convenience copy got made," not "the log got lost." No-ops
# if the job never reached run_and_save_outer_optimization_distributed (e.g. a precompile failure)
# -- grep finds nothing, nothing gets copied, log stays in logs/ as it always did.
LOG_FILE="logs/log-${SLURM_JOB_ID}-${SLURM_JOB_NAME}.out"
RUN_FOLDER=$(grep -m 1 "^run_folder" "$LOG_FILE" | sed 's/^run_folder[[:space:]]*=[[:space:]]*//')
if [ -n "$RUN_FOLDER" ] && [ -d "$RUN_FOLDER" ]; then
    cp "$LOG_FILE" "$RUN_FOLDER/log.out"
fi
