#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load samtools

DONOR=$1
FILE_LIST=$2

DIR="/scratch1/yaojiani/processed/eQTLAutoimmune/03_mergeBAM"


# ------------------------------------------------
# Step 3: Merge and sort bams for each donor
# -----------------------------------------------
echo "Step 3: Merge and sort bams for each donor"

## Merge
samtools merge -b $FILE_LIST -o $DIR/$DONOR.merged.step3.bam
## Sort
samtools sort $DIR/$DONOR.merged.step3.bam -o $DIR/$DONOR.step3.bam

rm $DIR/$DONOR.merged.step3.bam

echo "Merge and sort bams successfully!"