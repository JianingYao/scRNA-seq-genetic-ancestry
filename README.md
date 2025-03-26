# Evaluating genetic-ancestry inference from single-cell RNA-seq data

Characterizing the ancestry of donors in single-cell RNA sequencing (scRNA-seq) studies is critical to ensure the genetic homogeneity of the dataset and reduce biases in analyses, to identify ancestry-specific regulatory mechanisms and understand their downstream role in diseases, and to ensure that existing datasets are representative of human genetic diversity. While scRNA-seq is now widely available, the information on the ancestry of the donors is often missing, hindering further analysis. Here we propose a framework to evaluate methods for inferring genetic-ancestry from genetic polymorphisms detected from scRNA-seq reads. We demonstrate that widely used tools (e.g., ADMIXTURE) provide accurate inference of genetic-ancestry and admixture proportions despite the limited number of genetic polymorphisms identified and imperfect variant calling from scRNA-seq reads. We inferred genetic-ancestry for 196 donors from four scRNA-seq datasets from the Human Cell Atlas and highlighted an extremely large proportion of donors of European ancestry. For researchers generating single-cell datasets, we recommend reporting genetic-ancestry inference for all donors and generating datasets that represent diverse ancestries.

**This repository documents all necessary scripts for the analysis in *Evaluating genetic-ancestry inference from single-cell RNA-seq data* paper.**


## Tutorial: from scRNA-seq reads to genetic-ancestry inference

[This tutorial](https://github.com/JianingYao/scRNA-seq-genetic-ancestry/tree/main/Tutorial) documents the steps to infer genetic-ancestry for donors from a scRNA-seq dataset. Please refer to `README_Tutorial.md` for a detailed description of how to use these scripts.

In [stage 1](https://github.com/JianingYao/scRNA-seq-genetic-ancestry/tree/main/Tutorial/Stage1_preprocessing), we preprocess FASTQ files and perform variant calling using the GATK RNAseq short variant discovery workflow. First, the sequencing reads are aligned to the human hg19 reference genome using STARSolo in `Stage1_1_preprocessing_cell_level.sh`. We demonstrate the mapping procedure using a 10X 3' v2 dataset -- please adjust accordingly if your sequencing reads are generated with a different technology. Second, all BAM files corresponding to the same donor are merged and processed to generate variants in `Stage1_2_preprocessing_donor_level.sh`. Finally, joint genotyping of all donors in the dataset is performed, and SNPs are filtered prior to downstream analyses in `Stage1_3_preprocessing_study_level.sh`. 

In [stage 2](https://github.com/JianingYao/scRNA-seq-genetic-ancestry/tree/main/Tutorial/Stage2_ancestry_inference), we infer the genetic ancestry of donors using SNPs identified from scRNA-seq reads. In the tutorial, we consider the harmonized HGDP+1kGP dataset as a reference population dataset for genetic-ancestry inference and select six genetic-ancestry groups (Africa, America, Europe, Middle East, East Asia, and South Asia). *Note, please follow the script at `../../Analysis/00_TGP/run_TGP.sh` to first create and QC the reference HGDP+1kGP dataset, or select an alternative dataset of your choice, before running `Stage2_ancestry_inference.sh`.* First, we convert the VCF file from Stage 1 into PLINK format and perform data quality control. After restricting to common SNPs present in both the scRNA-seq and reference datasets and performing genetic pruning, we use the three approaches PCA-Distance, PC-RandomForest, and ADMIXTURE to infer genetic-ancestry of donors. 

**Requirements in the tutorial:**
- STARsolo (v2.7.11a)
- picard (2.26.2)
- samtools
- gatk (4.2.6.1)
- PLINK (1.9)
- ADMIXTURE (1.3.0)
- R (RColorBrewer, randomForest, dplyr)

## Citation

To complete

## Contact 
Jianing Yao: jyao37@jhmi.edu

Steven Gazal: gazal@usc.edu
