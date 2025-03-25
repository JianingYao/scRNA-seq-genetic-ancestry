#!/bin/bash

#SBATCH --partition=epyc-64
#SBATCH --mem=100G
#SBATCH --time=10:00:00


# Generating genome indexes

DIR="/project/gazal_569/DATA/gene_info/hg19"
STAR="/project/gazal_569/jianing/sc-RNA_ancestry/bin/STAR/source/STAR"
$STAR \
--runMode genomeGenerate \
--genomeDir star_index_HnsccImmuneLandscape \
--genomeFastaFiles $DIR/hg19.fa \
--sjdbGTFfile $DIR/hg19.ensGene.gtf \
--sjdbOverhang 97