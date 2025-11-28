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
#SBATCH --job-name merge_DXAK

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
mkdir -p ESYP

samtools merge ESYP/ESYP_merge_beta.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v045/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_T2D/HUB-NG-v045HUB-NG-v045.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v046/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_T2D/HUB-NG-v046HUB-NG-v046.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v047/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_T2D/HUB-NG-v047HUB-NG-v047.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v048/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_T2D/HUB-NG-v048HUB-NG-v048.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v049/bam_hiv2_hgr38/extract/newAnnotations_hg38_Beta_T2D/HUB-NG-v049HUB-NG-v049.bam \


samtools merge ESYP/ESYP_merge_alpha.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v045/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_T2D/HUB-NG-v045HUB-NG-v045.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v046/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_T2D/HUB-NG-v046HUB-NG-v046.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v047/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_T2D/HUB-NG-v047HUB-NG-v047.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v048/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_T2D/HUB-NG-v048HUB-NG-v048.bam \
/no_backup/jferrer/mplanas/VASAseq/processed_2/HUB-NG-v049/bam_hiv2_hgr38/extract/newAnnotations_hg38_Alpha_T2D/HUB-NG-v049HUB-NG-v049.bam \


###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
