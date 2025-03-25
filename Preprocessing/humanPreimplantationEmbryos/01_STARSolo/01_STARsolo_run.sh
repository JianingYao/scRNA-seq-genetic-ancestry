#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for dir in /scratch1/yaojiani/humanPreimplantationEmbryos/data/*
do 
    file1=$dir/*
    name="$(basename $file1 | awk -F'.fastq.gz' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/01_STARSolo/01_STARsolo_humanPreimplantationEmbryos.sh $file1 $name $dir
done