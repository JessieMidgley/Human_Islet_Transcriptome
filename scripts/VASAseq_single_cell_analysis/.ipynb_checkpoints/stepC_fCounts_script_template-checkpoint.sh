#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepE_feature_counting/__SAMPLE__/%x_%j.out
#SBATCH --error=/no_backup/jferrer/jmidgley/VASAseq/RERUN/stepE_feature_counting/__SAMPLE__/%x_%j.err

# time limit in minutes
#SBATCH --time=2:59:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=20GB

# job name
#SBATCH --job-name fCounts___SAMPLE__

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
/no_backup/jferrer/jmidgley/VASAseq/stepE_feature_counting/subread-2.1.1-Linux-x86_64/bin/featureCounts \
-s 1 -T 4 -R BAM  \
-a /no_backup/jferrer/jmidgley/VASAseq/jessies_HIT_v3.0.gtf \
-o /no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/fCounts.txt \
/no_backup/jferrer/jmidgley/VASAseq/RERUN/processed/__SAMPLE__/bam/Aligned.sortedByCoord.out.bam


###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $peak_mem bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
