#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for FILE in /scratch1/yaojiani/processed/eQTLAutoimmune/03_mergeBAM/*
do
    FILE_NAME="$(basename $FILE | awk -F'.step3.bam' '{print $1}')"
    DIR_IN="/scratch1/yaojiani/processed/eQTLAutoimmune/03_mergeBAM/"
    DIR_OUT="/scratch1/yaojiani/processed/eQTLAutoimmune/GATK_part1/"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/GATK_pipeline/part1_bams.sh $FILE_NAME $DIR_IN $DIR_OUT
done


