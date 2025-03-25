#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM in /scratch2/yaojiani/processed/humanPreimplantationEmbryos/03_mergeBAM/*.step3.bam
do 
    FILE_NAME="$(basename $BAM | awk -F'.step3.bam' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/04_MarkDuplicates/04_MarkDuplicates_humanPreimplantationEmbryos.sh $FILE_NAME
done