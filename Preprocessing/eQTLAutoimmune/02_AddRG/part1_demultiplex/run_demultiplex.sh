#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for BAM_FILE in /scratch1/yaojiani/processed/eQTLAutoimmune/01_star_alignment/*.Aligned.sortedByCoord.out.bam
do 
    SRA="$(basename $BAM_FILE | awk -F'.' '{print $1}')"
    for BARCODES in /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/02_AddRG/part1_demultiplex/*_barcodes.txt
    do
        DONOR="$(basename $BARCODES | awk -F'_barcodes.txt' '{print $1}')"
        sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/02_AddRG/part1_demultiplex/eQTLAutoimmune_demultiplex.sh $BAM_FILE $BARCODES $DONOR $SRA
    done
done

