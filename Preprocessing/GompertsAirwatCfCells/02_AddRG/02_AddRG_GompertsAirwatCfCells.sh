#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00


export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load picard/2.26.2
module load samtools

FILE=$1
RGLB=$2
RGPL=$3
RGPU=$4
RGSM=$5
RGDS=$6

INPUT=/scratch1/yaojiani/processed/GompertsAirwatCfCells/01_star_alignment/$FILE.step1.Aligned.sortedByCoord.out.bam
FIX=/scratch1/yaojiani/processed/GompertsAirwatCfCells/02_AddRG/$FILE.fixed.bam
OUTPUT=/scratch1/yaojiani/processed/GompertsAirwatCfCells/02_AddRG/$FILE.step2.bam

samtools view -h $INPUT | grep -v "^@RG" | sed "s/\tRG:Z:[^\t]*//" | samtools view -bo $FIX -

# ------------------------------------------------
# Step 2: AddOrReplaceReadGroups
# -----------------------------------------------
echo "Step: AddOrReplaceReadGroups"
picard AddOrReplaceReadGroups \
    I=$FIX \
    O=$OUTPUT \
    RGID=$FILE \
    RGLB=$RGLB \
    RGPL=$RGPL \
    RGPU=$RGPU \
    RGSM=$RGSM \
    RGDS=$RGDS

rm $FIX
mv $INPUT /scratch1/yaojiani/processed/GompertsAirwatCfCells/01_star_alignment/complete

echo "Finish Step 2: AddOrReplaceReadGroups successfully!"
