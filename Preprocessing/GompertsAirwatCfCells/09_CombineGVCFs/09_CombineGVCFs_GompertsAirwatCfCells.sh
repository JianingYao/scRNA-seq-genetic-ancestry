#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00


export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

module load gatk/4.2.6.1

STUDY="GompertsAirwatCfCells"
DIR_IN="/scratch1/yaojiani/processed/GompertsAirwatCfCells/08_HaplotypeCaller"
DIR_OUT="/scratch1/yaojiani/processed/GompertsAirwatCfCells/09_CombineGVCFs"
HG19_DIR="/project/gazal_569/DATA/gene_info/hg19"

find $DIR_IN -type f -name "*.interval.hgdp_1kg.g.vcf.gz" | sort > $DIR_OUT/$STUDY.gvcf.list

# ------------------------------------------------
# Step: CombineGVCFs
# ------------------------------------------------
echo "Step: CombineGVCFs"

gatk CombineGVCFs \
    -R $HG19_DIR/hg19.fa \
    $(awk '{print "--variant " $0}' $DIR_OUT/$STUDY.gvcf.list) \
    -O $DIR_OUT/$STUDY.g.vcf.gz

while read file; do
    mv "$file" $DIR_IN/complete/
done < $DIR_OUT/$STUDY.gvcf.list

echo "CombineGVCFs successfully"

