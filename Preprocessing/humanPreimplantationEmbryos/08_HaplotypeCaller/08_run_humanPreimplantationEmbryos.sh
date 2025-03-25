#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM in /scratch2/yaojiani/processed/humanPreimplantationEmbryos/07_ApplyBQSR/*.step7.bam
do 
    FILE_NAME="$(basename $BAM | awk -F'.step7.bam' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/08_HaplotypeCaller/08_HaplotypeCaller_humanPreimplantationEmbryos.sh $FILE_NAME
done