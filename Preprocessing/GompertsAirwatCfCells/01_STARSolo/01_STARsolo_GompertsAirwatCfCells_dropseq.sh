#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00


# set tmp directory
export TMPDIR=/scratch2/yaojiani/job_tmp
# if error occurs, exit
set -e

FASTQ_INR1=$1
FASTQ_INR2=$2
FILE_OUT=$3
DIR_R1=$4
DIR_R2=$5

STAR="/project/gazal_569/jianing/sc-RNA_ancestry/bin/STAR/source/STAR"
STAR_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/star_index_GompertsAirwatCfCells"
SOLO="CB_UMI_Simple"
WHITELIST=None
UMIlen=8
UMIstart=13
CBlen=12
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
--outFileNamePrefix "/scratch1/yaojiani/processed/GompertsAirwatCfCells/01_star_alignment/$FILE_OUT.step1." \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattrRGline ID:$FILE_OUT \
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
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--soloCellFilter EmptyDrops_CR \
--clipAdapterType CellRanger4


mv $DIR_R1 /scratch1/yaojiani/GompertsAirwatCfCells/
mv $DIR_R2 /scratch1/yaojiani/GompertsAirwatCfCells/








