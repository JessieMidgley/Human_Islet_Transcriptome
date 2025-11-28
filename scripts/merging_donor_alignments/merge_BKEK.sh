#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/logs/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/logs/%x_%j.err

# time limit in minutes
#SBATCH --time=1:00:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=15GB

# job name
#SBATCH --job-name merge_BKEK

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
module load SAMtools

#################################################
# run command 
#################################################
mkdir -p BKEK

samtools merge BKEK/BKEK_merge_beta.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v033/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_ND/HUB-NG-v033HUB-NG-v033.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v034/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_ND/HUB-NG-v034HUB-NG-v034.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v036/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_ND/HUB-NG-v036HUB-NG-v036.bam \

samtools merge BKEK/BKEK_merge_alpha.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v033/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_ND/HUB-NG-v033HUB-NG-v033.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v034/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_ND/HUB-NG-v034HUB-NG-v034.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v036/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_ND/HUB-NG-v036HUB-NG-v036.bam \



###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
