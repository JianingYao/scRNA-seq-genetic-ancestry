#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

FILE_NAME=$1

DIR_IN="/scratch1/yaojiani/processed/HnsccImmuneLandscape/04_MarkDuplicates"
DIR_OUT="/scratch1/yaojiani/processed/HnsccImmuneLandscape/05_SplitNCigarReads"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"


# ------------------------------------------------
# Step: SplitNCigarReads 
# ------------------------------------------------
echo "Step: SplitNCigarReads"

gatk SplitNCigarReads \
    --tmp-dir /scratch2/yaojiani/job_tmp \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_IN/$FILE_NAME.step4.bam \
    -O $DIR_OUT/$FILE_NAME.step5.bam


mv $DIR_IN/$FILE_NAME.step4.bam $DIR_IN/complete/
echo "SplitNCigarReads successfully!"