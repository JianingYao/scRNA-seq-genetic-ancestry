#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

STUDY="eQTLAutoimmune"
DIR_IN="/scratch1/yaojiani/processed/eQTLAutoimmune/GATK_part1"
DIR_OUT="/scratch1/yaojiani/processed/eQTLAutoimmune/GATK_part2"

sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/GATK_pipeline/part2_vcfs.sh $STUDY $DIR_IN $DIR_OUT