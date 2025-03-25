#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

FILE_NAME=$1

DIR_IN="/scratch1/yaojiani/processed/HnsccImmuneLandscape/05_SplitNCigarReads"
DIR_OUT="/scratch1/yaojiani/processed/HnsccImmuneLandscape/06_BaseRecalibrator"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"
SITES_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/known-sites"

# ------------------------------------------------
# Step: BaseRecalibrator 
# ------------------------------------------------
echo "Step: BaseRecalibrator"

gatk BaseRecalibrator \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_IN/$FILE_NAME.step5.bam \
    -O $DIR_OUT/$FILE_NAME.step6.table \
    --known-sites $SITES_DIR/dbsnp_138.hg19.vcf \
    --known-sites $SITES_DIR/1000G_phase1.indels.hg19.sites.vcf \
    --known-sites $SITES_DIR/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf

mv $DIR_IN/$FILE_NAME.step5.bam $DIR_IN/complete/
echo "BaseRecalibrator successfully"


