#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepC_star/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepC_star/%x_%j.err

# time limit in minutes
#SBATCH --time=5:59:00

# queue
#SBATCH --qos=short

# memory (MB)
#SBATCH --mem=50GB

# job name
#SBATCH --job-name genomegenerate

# make bash behave more robustly
set -e
set -u
set -o pipefail

# email when job is finished
#SBATCH --mail-type=END
#SBATCH --mail-user=jessie.midgley@crg.eu

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

###################
# set environment #
###################
module load STAR

#################################################
# run command 
#################################################
STAR \
--runThreadN 6 --runMode genomeGenerate \
--genomeDir /no_backup/jferrer/jmidgley/VASAseq/RERUN/stepC_star/hg38_index \
--genomeFastaFiles /no_backup/jferrer/mplanas/VASAseq/files/hg38.fa \
--sjdbGTFfile /no_backup/jferrer/jmidgley/VASAseq/jessies_HIT_v3.0.gtf #--sjdbOverhang  99

###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
