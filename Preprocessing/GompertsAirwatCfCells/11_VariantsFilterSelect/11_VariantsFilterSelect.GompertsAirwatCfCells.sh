#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00


export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

STUDY="GompertsAirwatCfCells"
DIR_IN="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/10_GenotypeGVCFs"
DIR_OUT="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/11_VariantsFilterSelect"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"

# ------------------------------------------------
# Step: Variants Filtration and Selection
# ------------------------------------------------
echo "Step: Variants Filtration and Selection"

gatk SelectVariants \
   -R $HG19_DIR/hg19.fa \
   -V $DIR_IN/$STUDY.vcf.gz \
   --select-type-to-include SNP \
   -O $DIR_OUT/$STUDY.SNP.vcf.gz

gatk VariantFiltration \
   -R $HG19_DIR/hg19.fa \
   -V $DIR_OUT/$STUDY.SNP.vcf.gz \
   -filter "QD < 2.0" --filter-name "QD2" \
   -filter "QUAL < 30.0" --filter-name "QUAL30" \
   -filter "SOR > 3.0" --filter-name "SOR3" \
   -filter "FS > 60.0" --filter-name "FS60" \
   -filter "MQ < 40.0" --filter-name "MQ40" \
   -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
   -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
   -O $DIR_OUT/$STUDY.SNP.Filtered.vcf.gz

gatk SelectVariants \
   -R $HG19_DIR/hg19.fa \
   -V $DIR_OUT/$STUDY.SNP.Filtered.vcf.gz \
   --exclude-filtered true \
   -O $DIR_OUT/$STUDY.SNP.Filtered.SV.vcf.gz

echo "Variants Filtration and Selection successfully"




