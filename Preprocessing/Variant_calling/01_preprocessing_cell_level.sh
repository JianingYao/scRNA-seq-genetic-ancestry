#!/bin/bash

# Required inputs
FASTQ_INR1=$1 # read 1
FASTQ_INR2=$2 # read 2
FILE=$3 # file output name
RGLB=$4 # read group information
RGPL=$5 # read group information
RGPU=$6 # read group information
RGSM=$7 # read group information
RGDS=$8 # read group information

module load picard/2.26.2
module load samtools
STAR="./bin/STAR/source/STAR"

# set tmp directory
export TMPDIR=./job_tmp

# ------------------------------------------------
# Step 0: Generating genome indexes
# ------------------------------------------------
# Note: This step only needs to be performed ONCE before using STARsolo for mapping. Comment out this section after the first round. 
DIR="/project/gazal_569/DATA/gene_info/hg19"

$STAR \
--runMode genomeGenerate \
--genomeDir star_index \
--genomeFastaFiles $DIR/hg19.fa \
--sjdbGTFfile $DIR/hg19.ensGene.gtf \
--sjdbOverhang 100

echo "Finish Step 0: Generating genome indexes successfully!"


# ------------------------------------------------
# Step 1: STARsolo mapping
# ------------------------------------------------
# Note: subject to change for different sequencing technologies.
STAR_DIR="./star_index"
SOLO="CB_UMI_Simple"
WHITELIST="./bin/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt"
UMIlen=10
UMIstart=17
CBlen=16
CBstart=1
STRAND="Forward"
NODE=1
echo $NODE

$STAR \
--runMode alignReads \
--genomeDir $STAR_DIR \
--readFilesCommand zcat \
--readFilesIn $FASTQ_INR2 $FASTQ_INR1 \
--runDirPerm All_RWX \
--runThreadN $NODE \
--outFileNamePrefix "$FILE.step1." \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattrRGline ID:$FILE \
--twopassMode Basic \
--twopass1readsN -1 \
--outSAMattributes NH HI nM AS MD XS CR CY UR UY CB UB sS sM sQ \
--soloType $SOLO \
--soloCBwhitelist $WHITELIST \
--soloUMIlen $UMIlen \
--soloUMIstart $UMIstart \
--soloCBlen $CBlen \
--soloCBstart $CBstart \
--soloBarcodeReadLength 0 \
--soloStrand $STRAND \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--soloCellFilter EmptyDrops_CR \
--clipAdapterType CellRanger4

echo "Finish Step 1: STARsolo mapping successfully!"


# ------------------------------------------------
# Step 2: AddOrReplaceReadGroups
# ------------------------------------------------
INPUT=$FILE.step1.Aligned.sortedByCoord.out.bam
FIX=$FILE.fixed.bam
OUTPUT=$FILE.step2.bam

samtools view -h $INPUT | grep -v "^@RG" | sed "s/\tRG:Z:[^\t]*//" | samtools view -bo $FIX -
echo "Step 2: AddOrReplaceReadGroups"
picard AddOrReplaceReadGroups \
    I=$FIX \
    O=$OUTPUT \
    RGID=$FILE \
    RGLB=$RGLB \
    RGPL=$RGPL \
    RGPU=$RGPU \
    RGSM=$RGSM \
    RGDS=$RGDS

rm $FIX
echo "Finish Step 2: AddOrReplaceReadGroups successfully!"
