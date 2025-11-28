#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/whippet/RERUN/logs/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/whippet/RERUN/logs/%x_%j.err

# time limit in minutes
#SBATCH --time=3:00:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=36GB

# job name
#SBATCH --job-name whippet_quantify

# make bash behave more robustly
set -e
set -u
set -o pipefail

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

###################
# set environment #
###################
# Set paths
JULIA_BIN="/users/jferrer/jmidgley/.julia/juliaup/julia-1.6.7+0.x64.linux.gnu/bin/julia"
WHIPPET_PROJECT="/users/jferrer/jmidgley/Whippet.jl"
CUSTOM_DEPOT="/no_backup/jferrer/jmidgley/julia_depot"

# Set environment
export JULIA_DEPOT_PATH="$CUSTOM_DEPOT"
export JULIA_PROJECT="$WHIPPET_PROJECT"

for dir in AFLO BKEK DBHQ GPRL HUEN; do
    for fastq in /no_backup/jferrer/jmidgley/${dir}/*.fastq.gz; do
        $JULIA_BIN $WHIPPET_PROJECT/bin/whippet-quant.jl "$fastq" \
        -o /no_backup/jferrer/jmidgley/whippet/RERUN/results/$(basename "$fastq" .fastq.gz) \
        -x /no_backup/jferrer/jmidgley/whippet/RERUN/index.jls \
        --sam > /no_backup/jferrer/jmidgley/whippet/RERUN/results/$(basename "$fastq" .fastq.gz).sam
    done
done

###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
