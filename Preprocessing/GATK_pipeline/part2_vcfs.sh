#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

STUDY=$1
DIR_IN=$2
DIR_OUT=$3

HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"

find $DIR_IN -type f -name "*.interval.hgdp_1kg.g.vcf.gz" | sort > $DIR_OUT/$STUDY.gvcf.list

# ------------------------------------------------
# Step 9: CombineGVCFs
# ------------------------------------------------
echo "Step 9: CombineGVCFs"

gatk CombineGVCFs \
    -R $HG19_DIR/hg19.fa \
    $(awk '{print "--variant " $0}' $DIR_OUT/$STUDY.gvcf.list) \
    -O $DIR_OUT/$STUDY.g.vcf.gz

echo "Step 9: CombineGVCFs successfully"


# ------------------------------------------------
# Step 10: GenotypeGVCFs
# ------------------------------------------------
echo "Step 10: GenotypeGVCFs"

gatk GenotypeGVCFs \
   -R $HG19_DIR/hg19.fa \
   -V $DIR_OUT/$STUDY.g.vcf.gz \
   -O $DIR_OUT/$STUDY.vcf.gz \
   --tmp-dir $TMPDIR

# rm $DIR_OUT/$STUDY.g.vcf.gz 
echo "Step 10: GenotypeGVCFs successfully"


# ------------------------------------------------
# Step 11: Variants Filtration and Selection
# ------------------------------------------------
echo "Step 11: Variants Filtration and Selection"

gatk SelectVariants \
   -R $HG19_DIR/hg19.fa \
   -V $DIR_OUT/$STUDY.vcf.gz \
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

echo "Step 11: Variants Filtration and Selection successfully"
