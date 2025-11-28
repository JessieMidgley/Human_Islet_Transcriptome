#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepG_count_tables/__SAMPLE__/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepG_count_tables/__SAMPLE__/%x_%j.err

# time limit in minutes
#SBATCH --time=00:59:00

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=12GB

# job name
#SBATCH --job-name countT___SAMPLE__

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
/users/jferrer/jmidgley/miniforge3/envs/scm/bin/bamToCountTable.py \
--noNames /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.featureCounts.sort.tag_dep.out.bam \
-o /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/__SAMPLE__.output.countTable \
-joinedFeatureTags XT,chrom -sampleTags SM 

###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
