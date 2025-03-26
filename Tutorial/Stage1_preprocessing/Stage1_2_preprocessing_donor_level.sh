#!/bin/bash

export TMPDIR=./job_tmp

module load picard/2.26.2
module load gatk/4.2.6.1
module load samtools

DONOR=$1 # one donor name
FILE_LIST=$2 # bam files from step 2 to be merged for this specific donor

HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"
SITES_DIR="./known-sites"


# ------------------------------------------------
# Step 3: Merge and sort bams for each donor
# ------------------------------------------------
echo "Step: Merge and sort bams for each donor"

samtools merge -b $FILE_LIST -o $DONOR.merged.step3.bam
samtools sort $DONOR.merged.step3.bam -o $DONOR.step3.bam
rm $DONOR.merged.step3.bam

echo "Finish Step 3: Merge and sort bams successfully!"


# ------------------------------------------------
# Step 4: Mark duplicates
# ------------------------------------------------
echo "Step 4: Mark duplicate"

picard MarkDuplicates \
    I=$DONOR.step3.bam \
    O=$DONOR.step4.bam \
    M=$DONOR.step4.marked_dup_metrics.txt \
    VALIDATION_STRINGENCY=STRICT

echo "Step 4: Mark duplicates successfully!"


# ------------------------------------------------
# Step 5: SplitNCigarReads 
# ------------------------------------------------
echo "Step 5: SplitNCigarReads"

gatk SplitNCigarReads \
    --tmp-dir $TMPDIR \
    -R $HG19_DIR/hg19.fa \
    -I $DONOR.step4.bam \
    -O $DONOR.step5.bam

echo "Step 5: SplitNCigarReads successfully!"


# ------------------------------------------------
# Step 6: BaseRecalibrator 
# ------------------------------------------------
echo "Step 6: BaseRecalibrator"

gatk BaseRecalibrator \
    -R $HG19_DIR/hg19.fa \
    -I $DONOR.step5.bam \
    -O $DONOR.step6.table \
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
    -I $DONOR.step5.bam \
    --bqsr-recal-file $DONOR.step6.table \
    -O $DONOR.step7.bam

echo "Step 7: ApplyBQSR successfully"


# ------------------------------------------------
# Step 8: HaplotypeCaller
# ------------------------------------------------
echo "Step 8: HaplotypeCaller"

gatk HaplotypeCaller \
    -R $HG19_DIR/hg19.fa \
    -I $DONOR.step7.bam \
    -O $DONOR.interval.hgdp_1kg.g.vcf.gz \
    -L /project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.interval.bed \
    -ERC GVCF

echo "Step 8: HaplotypeCaller successfully"