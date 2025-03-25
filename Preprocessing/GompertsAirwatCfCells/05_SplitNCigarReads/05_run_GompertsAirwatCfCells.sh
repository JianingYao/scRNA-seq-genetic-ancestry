#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM in /scratch1/yaojiani/processed/GompertsAirwatCfCells/04_MarkDuplicates/*.step4.bam
do 
    FILE_NAME="$(basename $BAM | awk -F'.step4.bam' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/05_SplitNCigarReads/05_SplitNCigarReads_GompertsAirwatCfCells.sh $FILE_NAME
done