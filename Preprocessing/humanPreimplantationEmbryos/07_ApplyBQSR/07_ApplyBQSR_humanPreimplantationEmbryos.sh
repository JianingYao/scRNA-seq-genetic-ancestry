#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00

 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

FILE_NAME=$1

DIR_IN_1="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/05_SplitNCigarReads/complete"
DIR_IN_2="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/06_BaseRecalibrator"
DIR_OUT="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/07_ApplyBQSR"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"


# ------------------------------------------------
# Step: ApplyBQSR
# ------------------------------------------------
echo "Step: ApplyBQSR"

gatk ApplyBQSR \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_IN_1/$FILE_NAME.step5.bam \
    --bqsr-recal-file $DIR_IN_2/$FILE_NAME.step6.table \
    -O $DIR_OUT/$FILE_NAME.step7.bam


mv $DIR_IN_2/$FILE_NAME.step6.table $DIR_IN_2/complete/
echo "ApplyBQSR successfully"





