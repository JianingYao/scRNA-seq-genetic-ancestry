#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00

export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load picard/2.26.2

FILE=$1
DONOR=$2
SRA=$3

RGID=${DONOR}_${SRA}
RGPU=$SRA
RGSM=$DONOR

OUTPUT=/scratch1/yaojiani/processed/eQTLAutoimmune/02_AddRG/$DONOR/${DONOR}_${SRA}.step2.bam

# ------------------------------------------------
# Step 2: AddOrReplaceReadGroups
# -----------------------------------------------
echo "Step: AddOrReplaceReadGroups"
picard AddOrReplaceReadGroups \
    I=$FILE \
    O=$OUTPUT \
    RGID=$RGID \
    RGLB=pool1 \
    RGPL=ILLUMINA \
    RGPU=$RGPU \
    RGSM=$RGSM

rm $FILE

echo "AddOrReplaceReadGroups successfully"
