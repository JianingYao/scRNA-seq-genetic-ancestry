#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00

 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load samtools

DONOR=$1
FILE_LIST=$2

DIR="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/03_mergeBAM"


# ------------------------------------------------
# Step: Merge and sort bams for each donor
# -----------------------------------------------
echo "Step: Merge and sort bams for each donor"

## Merge
samtools merge -b $FILE_LIST -o $DIR/$DONOR.merged.step3.bam
## Sort
samtools sort $DIR/$DONOR.merged.step3.bam -o $DIR/$DONOR.step3.bam

rm $DIR/$DONOR.merged.step3.bam

mv $FILE_LIST /project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/03_mergeBAM/donors/
echo "Merge and sort bams successfully!"