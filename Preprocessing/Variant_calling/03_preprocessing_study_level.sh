#!/bin/bash

export TMPDIR=./job_tmp

module load gatk/4.2.6.1

STUDY=$1

HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"

find . -type f -name "*.interval.hgdp_1kg.g.vcf.gz" | sort > $STUDY.gvcf.list

# ------------------------------------------------
# Step 9: CombineGVCFs
# ------------------------------------------------
echo "Step 9: CombineGVCFs"

gatk CombineGVCFs \
    -R $HG19_DIR/hg19.fa \
    $(awk '{print "--variant " $0}' $STUDY.gvcf.list) \
    -O $STUDY.g.vcf.gz

echo "Step 9: CombineGVCFs successfully"


# ------------------------------------------------
# Step 10: GenotypeGVCFs
# ------------------------------------------------
echo "Step 10: GenotypeGVCFs"

gatk GenotypeGVCFs \
   -R $HG19_DIR/hg19.fa \
   -V $STUDY.g.vcf.gz \
   -O $STUDY.vcf.gz \
   --tmp-dir $TMPDIR

echo "Step 10: GenotypeGVCFs successfully"


# ------------------------------------------------
# Step 11: Variants Filtration and Selection
# ------------------------------------------------
echo "Step 11: Variants Filtration and Selection"

gatk SelectVariants \
   -R $HG19_DIR/hg19.fa \
   -V $STUDY.vcf.gz \
   --select-type-to-include SNP \
   -O $STUDY.SNP.vcf.gz

gatk VariantFiltration \
   -R $HG19_DIR/hg19.fa \
   -V $STUDY.SNP.vcf.gz \
   -filter "QD < 2.0" --filter-name "QD2" \
   -filter "QUAL < 30.0" --filter-name "QUAL30" \
   -filter "SOR > 3.0" --filter-name "SOR3" \
   -filter "FS > 60.0" --filter-name "FS60" \
   -filter "MQ < 40.0" --filter-name "MQ40" \
   -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
   -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
   -O $STUDY.SNP.Filtered.vcf.gz

gatk SelectVariants \
   -R $HG19_DIR/hg19.fa \
   -V $STUDY.SNP.Filtered.vcf.gz \
   --exclude-filtered true \
   -O $STUDY.SNP.Filtered.SV.vcf.gz

echo "Step 11: Variants Filtration and Selection successfully"
