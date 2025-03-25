#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load picard/2.26.2

FILE_NAME=$1

DIR_IN="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/03_mergeBAM"
DIR_OUT="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/04_MarkDuplicates"


# ------------------------------------------------
# Step: Mark duplicates
# ------------------------------------------------
echo "Step: Mark duplicate"

# remove failed output
if [ -f $DIR_OUT/$FILE_NAME.step4.bam ]; then
    rm $DIR_OUT/$FILE_NAME.step4.bam
fi

picard MarkDuplicates \
    I=$DIR_IN/$FILE_NAME.step3.bam \
    O=$DIR_OUT/$FILE_NAME.step4.bam \
    M=$DIR_OUT/$FILE_NAME.step4.marked_dup_metrics.txt \
    VALIDATION_STRINGENCY=STRICT

mv $DIR_IN/$FILE_NAME.step3.bam $DIR_IN/complete/
echo "Mark duplicates successfully!"