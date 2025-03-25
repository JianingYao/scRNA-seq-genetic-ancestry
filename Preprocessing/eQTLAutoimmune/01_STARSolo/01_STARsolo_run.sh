#!/bin/bash

# set tmp directory
export TMPDIR=/scratch2/yaojiani/job_tmp

for dir in /scratch1/yaojiani/eQTLAutoimmune/pool1/*
do 
    file1=$dir/*_2*
    file2=$dir/*_3*
    name="$(basename $file1 | awk -F'[_]' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/01_STARSolo/01_STARsolo_eQTLAutoimmune.sh $file1 $file2 $name
done

