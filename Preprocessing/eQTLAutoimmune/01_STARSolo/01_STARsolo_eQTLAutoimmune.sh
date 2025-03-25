#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00

# set tmp directory
export TMPDIR=/scratch2/yaojiani/job_tmp
echo "New TMPDIR: $TMPDIR"
# if error occurs, exit
set -e

FASTQ_INR1=$1
FASTQ_INR2=$2
FILE_OUT=$3

STAR="/project/gazal_569/jianing/sc-RNA_ancestry/bin/STAR/source/STAR"
STAR_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/star_index_eQTLAutoimmune"
SOLO="CB_UMI_Simple"
WHITELIST="/home1/yaojiani/bin/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt"
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
--outFileNamePrefix "/scratch1/yaojiani/processed/eQTLAutoimmune/01_star_alignment/$FILE_OUT.step1." \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattrRGline ID:$FILE_OUT LB:$FILE_OUT PL:ILLUMINA SM:pool1 PU:$FILE_OUT DS:10x3v2 \
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


