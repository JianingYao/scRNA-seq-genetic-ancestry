#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00


export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

STUDY="GompertsAirwatCfCells"
DIR_IN="/scratch1/yaojiani/processed/GompertsAirwatCfCells/09_CombineGVCFs"
DIR_OUT="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/10_GenotypeGVCFs"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"

# ------------------------------------------------
# Step: GenotypeGVCFs
# ------------------------------------------------
echo "Step: GenotypeGVCFs"

gatk GenotypeGVCFs \
   -R $HG19_DIR/hg19.fa \
   -V $DIR_IN/$STUDY.g.vcf.gz \
   -O $DIR_OUT/$STUDY.vcf.gz \
   --tmp-dir $TMPDIR

echo "GenotypeGVCFs successfully"

