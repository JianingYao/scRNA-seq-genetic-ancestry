#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load picard/2.26.2
module load gatk/4.2.6.1

FILE_NAME=$1
DIR_IN=$2
DIR_OUT=$3

HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"
SITES_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/known-sites"


# ------------------------------------------------
# Step 4: Mark duplicates
# ------------------------------------------------
echo "Step 4: Mark duplicate"

picard MarkDuplicates \
    I=$DIR_IN/$FILE_NAME.step3.bam \
    O=$DIR_OUT/$FILE_NAME.step4.bam \
    M=$DIR_OUT/$FILE_NAME.step4.marked_dup_metrics.txt \
    VALIDATION_STRINGENCY=STRICT

# rm $DIR_IN/$FILE_NAME.step3* 
echo "Step 4: Mark duplicates successfully!"


# ------------------------------------------------
# Step 5: SplitNCigarReads 
# ------------------------------------------------
echo "Step 5: SplitNCigarReads"

gatk SplitNCigarReads \
    --tmp-dir /scratch2/yaojiani/job_tmp \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_OUT/$FILE_NAME.step4.bam \
    -O $DIR_OUT/$FILE_NAME.step5.bam

# rm $DIR_OUT/$FILE_NAME.step4.*
echo "Step 5: SplitNCigarReads successfully!"


# ------------------------------------------------
# Step 6: BaseRecalibrator 
# ------------------------------------------------
echo "Step 6: BaseRecalibrator"

gatk BaseRecalibrator \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_OUT/$FILE_NAME.step5.bam \
    -O $DIR_OUT/$FILE_NAME.step6.table \
    --known-sites $SITES_DIR/dbsnp_138.hg19.vcf \
    --known-sites $SITES_DIR/1000G_phase1.indels.hg19.sites.vcf \
    --known-sites $SITES_DIR/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf

echo "Step 6: BaseRecalibrator successfully"


# ------------------------------------------------
# Step 7: ApplyBQSR
# ------------------------------------------------
echo "Step 7: ApplyBQSR"

gatk ApplyBQSR \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_OUT/$FILE_NAME.step5.bam \
    --bqsr-recal-file $DIR_OUT/$FILE_NAME.step6.table \
    -O $DIR_OUT/$FILE_NAME.step7.bam

# rm $DIR_OUT/$FILE_NAME.step5.* $DIR_OUT/$FILE_NAME.step6.*
echo "Step 7: ApplyBQSR successfully"


# ------------------------------------------------
# Step 8: HaplotypeCaller
# ------------------------------------------------
echo "Step 8: HaplotypeCaller"

gatk HaplotypeCaller \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_OUT/$FILE_NAME.step7.bam \
    -O $DIR_OUT/$FILE_NAME.interval.hgdp_1kg.g.vcf.gz \
    -L /project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.interval.bed \
    -ERC GVCF

# rm $DIR_OUT/$FILE_NAME.step7.*
echo "Step 8: HaplotypeCaller successfully"

