#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00


export TMPDIR=/scratch2/yaojiani/job_tmp
set -e

FASTQ_INR1=$1
FILE_OUT=$2
DIR=$3

STAR="/project/gazal_569/jianing/sc-RNA_ancestry/bin/STAR/source/STAR"
STAR_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/star_index_humanPreimplantationEmbryos"
NODE=1
echo $NODE
 
$STAR \
--genomeDir $STAR_DIR \
--readFilesCommand zcat \
--readFilesIn $FASTQ_INR1 \
--outSAMattrRGline ID:$FILE_OUT \
--outFileNamePrefix "/scratch2/yaojiani/processed/humanPreimplantationEmbryos/01_star_alignment/$FILE_OUT.step1." \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS MD XS \
--twopassMode Basic \
--twopass1readsN -1 \
--soloType SmartSeq \
--soloUMIdedup Exact \
--soloStrand Unstranded \
--soloCellFilter None \
--runDirPerm All_RWX \
--runMode alignReads \
--runThreadN $NODE \


mv $DIR /scratch1/yaojiani/humanPreimplantationEmbryos/




