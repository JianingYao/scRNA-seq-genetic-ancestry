#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

 
export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

FILE_NAME=$1

DIR_IN="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/07_ApplyBQSR"
DIR_OUT="/scratch2/yaojiani/processed/humanPreimplantationEmbryos/08_HaplotypeCaller"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"


# ------------------------------------------------
# Step: HaplotypeCaller
# ------------------------------------------------
echo "Step: HaplotypeCaller"

gatk HaplotypeCaller \
    -R $HG19_DIR/hg19.fa \
    -I $DIR_IN/$FILE_NAME.step7.bam \
    -O $DIR_OUT/$FILE_NAME.interval.hgdp_1kg.g.vcf.gz \
    -L /project/gazal_569/DATA/HGDP_1KG_gnomad/plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.interval.bed \
    -ERC GVCF


mv $DIR_IN/$FILE_NAME.step7.bam $DIR_IN/complete/
echo "HaplotypeCaller successfully"


