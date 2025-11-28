#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepD_align/__SAMPLE__/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepD_align/__SAMPLE__/%x_%j.err

# time limit in minutes
#SBATCH --time=2:59:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=36GB

# job name
#SBATCH --job-name align___SAMPLE__

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
--runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 \
--outMultimapperOrder Random  --outSAMmultNmax 1 --readFilesCommand zcat \
--readFilesIn /no_backup/jferrer/mplanas/VASAseq/processed_2/__SAMPLE__/trimmed.R2.fastq.gz \
--genomeDir /no_backup/jferrer/jmidgley/VASAseq/RERUN/stepC_star/hg38_index/ \
--outFileNamePrefix /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/


###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
