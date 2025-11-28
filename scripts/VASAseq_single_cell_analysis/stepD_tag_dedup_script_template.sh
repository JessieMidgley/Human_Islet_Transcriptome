#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepF_tag_dedup/__SAMPLE__/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepF_tag_dedup/__SAMPLE__/%x_%j.err

# time limit in minutes
#SBATCH --time=5:59:00

# queue
#SBATCH --qos=short

# memory (MB)
#SBATCH --mem=24GB

# job name
#SBATCH --job-name tagD___SAMPLE__

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

#################################################
# run command 
#################################################
module load SAMtools

samtools sort /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.out.bam.featureCounts.bam -o /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.out.bam.featureCounts.sort.bam

samtools index /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.out.bam.featureCounts.sort.bam


/users/jferrer/jmidgley/miniforge3/envs/scm/bin/bamtagmultiome.py \
-method cs_feature_counts -umi_hamming_distance 0 \
/no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.out.bam.featureCounts.sort.bam \
-o /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.featureCounts.sort.tag_dep.out.bam

###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
