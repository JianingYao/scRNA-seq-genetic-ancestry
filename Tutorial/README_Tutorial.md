## Tutorial: from scRNA-seq reads to genetic-ancestry inference

[This tutorial](https://github.com/JianingYao/scRNA-seq-genetic-ancestry/tree/main/Tutorial) documents the steps to infer genetic-ancestry for donors from a scRNA-seq dataset. 

In [stage 1](https://github.com/JianingYao/scRNA-seq-genetic-ancestry/tree/main/Tutorial/Stage1_preprocessing), we preprocess FASTQ files and perform variant calling using the GATK RNAseq short variant discovery workflow. First, the sequencing reads are aligned to the human hg19 reference genome using STARSolo in `Stage1_1_preprocessing_cell_level.sh`. We demonstrate the mapping procedure using a 10X 3' v2 dataset -- please adjust accordingly if your sequencing reads are generated with a different technology. Second, all BAM files corresponding to the same donor are merged and processed to generate variants in `Stage1_2_preprocessing_donor_level.sh`. Finally, joint genotyping of all donors in the dataset is performed, and SNPs are filtered prior to downstream analyses in `Stage1_3_preprocessing_study_level.sh`. 

In [stage 2](https://github.com/JianingYao/scRNA-seq-genetic-ancestry/tree/main/Tutorial/Stage2_ancestry_inference), we infer the genetic ancestry of donors using SNPs identified from scRNA-seq reads. In the tutorial, we consider the harmonized HGDP+1kGP dataset as a reference population dataset for genetic-ancestry inference and select six genetic-ancestry groups (Africa, America, Europe, Middle East, East Asia, and South Asia). *Note, please follow the script at `../../Analysis/00_TGP/run_TGP.sh` to first create and QC the reference HGDP+1kGP dataset, or select an alternative dataset of your choice, before running `Stage2_ancestry_inference.sh`.* First, we convert the VCF file from Stage 1 into PLINK format and perform data quality control. After restricting to common SNPs present in both the scRNA-seq and reference datasets and performing genetic pruning, we use the three approaches PCA-Distance, PC-RandomForest, and ADMIXTURE to infer genetic-ancestry of donors. 

### Requirements in the tutorial:
- STARsolo (v2.7.11a)
- picard (2.26.2)
- samtools
- gatk (4.2.6.1)
- PLINK (1.9)
- ADMIXTURE (1.3.0)
- R (RColorBrewer, randomForest, dplyr)

### Procedure

**Stage 1: Preprocessing scRNA-seq reads and generating VCF files:** run `Stage1_1_preprocessing_cell_level.sh`, `Stage1_2_preprocessing_donor_level.sh`, and `Stage1_3_preprocessing_study_level.sh` in order. Note the required inputs for each script. Please adjust `Stage1_1_preprocessing_cell_level.sh` accordingly if your sequencing reads are generated with a different technology.

**Stage 2: Inferring genetic-ancestry using SNPs identified from scRNA-seq reads:** Please follow the script at ../../Analysis/00_TGP/run_TGP.sh to first create and QC the reference HGDP+1kGP dataset, or select an alternative dataset of your choice. Then, run `Stage2_ancestry_inference.sh`. The required inputs are study name and path to the corresponding VCF file (this is the output from Stage 1). 