#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

for FILE_LIST in /project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/03_mergeBAM/donors/lists/*
do 
    DONOR="$(basename $FILE_LIST | awk -F'.txt' '{print $1}')"
    sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/03_mergeBAM/03_merge_HnsccImmuneLandscape.sh $DONOR $FILE_LIST
done