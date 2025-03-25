#!/bin/bash

# set tmp directory
export TMPDIR=/scratch2/yaojiani/job_tmp

for dir in /scratch1/yaojiani/HnsccImmuneLandscape/data/*
do 
    file1=$dir/*_1*
    file2=$dir/*_2*
    name="$(basename $file1 | awk -F'[_]' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/01_STARSolo/01_STARsolo_HnsccImmuneLandscape.sh $file1 $file2 $name $dir
done

