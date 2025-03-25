#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for dir in /scratch1/yaojiani/StemCellsInChronicMyeloidLeukemia/data/*
do 
    file1=$dir/*_R1*
    name="$(basename $file1 | awk -F'[_]' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/StemCellsInChronicMyeloidLeukemia/01_STARSolo/01_STARsolo_StemCellsInChronicMyeloidLeukemia.sh $file1 $name $dir
done