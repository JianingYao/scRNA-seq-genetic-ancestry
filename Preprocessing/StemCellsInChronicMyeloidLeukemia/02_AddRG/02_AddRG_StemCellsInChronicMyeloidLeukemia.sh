#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=30G
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
RGPM=$7

INPUT=/scratch2/yaojiani/processed/StemCellsInChronicMyeloidLeukemia/01_star_alignment/$FILE.step1.Aligned.sortedByCoord.out.bam
FIX=/scratch2/yaojiani/processed/StemCellsInChronicMyeloidLeukemia/02_AddRG/$FILE.fixed.bam
OUTPUT=/scratch2/yaojiani/processed/StemCellsInChronicMyeloidLeukemia/02_AddRG/$FILE.step2.bam

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
    RGDS=$RGDS \
    RGPM=$RGPM

rm $FIX
mv $INPUT /scratch2/yaojiani/processed/StemCellsInChronicMyeloidLeukemia/01_star_alignment/complete

echo "Finish Step 2: AddOrReplaceReadGroups successfully!"
