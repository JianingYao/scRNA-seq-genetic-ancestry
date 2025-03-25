#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM in /scratch1/yaojiani/processed/HnsccImmuneLandscape/03_mergeBAM/*.step3.bam
do 
    FILE_NAME="$(basename $BAM | awk -F'.step3.bam' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/04_MarkDuplicates/04_MarkDuplicates_HnsccImmuneLandscape.sh $FILE_NAME
done