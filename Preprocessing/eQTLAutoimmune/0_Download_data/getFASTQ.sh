#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00

set -e

SRAID=$1

STUDY_DIR="/scratch1/yaojiani/eQTLAutoimmune/pool1/"
SRAFILE=$STUDY_DIR/$SRAID/$SRAID.sra

## download SRA
module load gcc/8.3.0
module load sratoolkit/2.11.0
prefetch $SRAID --max-size u --output-directory $STUDY_DIR

## convert to fastq files
echo "Begin to convert $SRAID"
fasterq-dump -O $STUDY_DIR/$SRAID --split-files --include-technical $SRAFILE
gzip $STUDY_DIR/$SRAID/*fastq 

echo "Finished downloading SRA and converting to fastq files of eQTLAutoimmune pool 1"