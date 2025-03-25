#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM in /scratch2/yaojiani/processed/StemCellsInChronicMyeloidLeukemia/04_MarkDuplicates/*.step4.bam
do 
    FILE_NAME="$(basename $BAM | awk -F'.step4.bam' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/StemCellsInChronicMyeloidLeukemia/05_SplitNCigarReads/05_SplitNCigarReads_StemCellsInChronicMyeloidLeukemia.sh $FILE_NAME
done