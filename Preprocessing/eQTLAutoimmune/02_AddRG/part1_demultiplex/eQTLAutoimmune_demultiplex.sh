#!/bin/bash

#SBATCH --partition=main
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00

# set tmp directory
export TMPDIR=/scratch2/yaojiani/job_tmp
# if error occurs, exit
set -e

BAM_FILE=$1
BARCODES=$2
DONOR=$3
SRA=$4
DIR_OUT="/scratch1/yaojiani/processed/eQTLAutoimmune/02_AddRG/$DONOR"

module load samtools

if [ ! -d "$DIR_OUT" ]; then
  mkdir "$DIR_OUT"
fi

samtools view -H $BAM_FILE > $DIR_OUT/${DONOR}_${SRA}.SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $BARCODES > $DIR_OUT/${DONOR}_${SRA}.filtered_SAM_body
cat $DIR_OUT/${DONOR}_${SRA}.SAM_header $DIR_OUT/${DONOR}_${SRA}.filtered_SAM_body > $DIR_OUT/${DONOR}_${SRA}.sam
samtools view -b $DIR_OUT/${DONOR}_${SRA}.sam > $DIR_OUT/${DONOR}_${SRA}.bam

rm $DIR_OUT/${DONOR}_${SRA}.SAM_header $DIR_OUT/${DONOR}_${SRA}.filtered_SAM_body $DIR_OUT/${DONOR}_${SRA}.sam

echo "Demultiplex the pool successfully"