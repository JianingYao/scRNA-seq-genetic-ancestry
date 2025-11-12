module load plink2
PLINK19="/project/gazal_569/soft/plink1.9/plink"
VCF="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/VCF/GTEx_donor_tissue.SNP.Filtered.SV.vcf.gz"
ORI="/project/gazal_569/DATA/GTEx/GTEx_v9/GTEx_Analysis_2021-02-11_v9_WGS_VCF_files/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz"

zcat $VCF | grep chr | grep -v contig | cut -f1-2 | awk '{print $1, $2-1, $2, $2}' OFS="\t" > gtex.hg19.bed
#create gtex.hg38.bed from https://genome.ucsc.edu/cgi-bin/hgLiftOver

#get exon_UTRS
cp /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/Results/exon_UTRS.bed exon_UTRS.hg19.bed
#create exon_UTRS.hg38.bed from https://genome.ucsc.edu/cgi-bin/hgLiftOver

#get sc data
$PLINK19 --vcf $VCF --double-id --make-bed --out sc_data 
perl convert_to_hg38.pl
cp sc_data.bim sc_data.hg19.bim
cp sc_data.hg38.bim sc_data.bim

#get wgs data
plink2 --vcf $ORI --make-bed --extract bed0 gtex.hg38.bed --keep donors.txt --out wgs_data
awk '{print $1, $1":"$4":"$5":"$6, $3, $4, $5, $6}' wgs_data.bim > wgs_data.1.bim
awk '{print $1, $1":"$4":"$6":"$5, $3, $4, $5, $6}' wgs_data.bim > wgs_data.2.bim

#SNPs in both files
R
data1 = as.character(read.table("sc_data.bim")[,2])
data2 = as.character(read.table("wgs_data.1.bim")[,2])
data3 = as.character(read.table("wgs_data.2.bim")[,2])
inter = intersect(data1,c(data2,data3))
write.table(inter,file="snps.txt",quote=F,col.names=F,row.names=F)
q()
n

#
plink2 --bfile sc_data  --extract snps.txt --make-bed --out tmp
grep C tmp.bim | grep G > AT_CG.txt
grep A tmp.bim | grep T >> AT_CG.txt
plink2 --bfile tmp --exclude AT_CG.txt --make-bed --out sc_data_clean
rm tmp.* 
#
cp wgs_data.bim wgs_data.old.bim
cp wgs_data.1.bim wgs_data.bim
plink2 --bfile wgs_data --extract sc_data_clean.bim --make-bed --out tmp1
cp wgs_data.2.bim wgs_data.bim
plink2 --bfile wgs_data --extract sc_data_clean.bim --make-bed --out tmp2
$PLINK19 --bfile tmp1 --bmerge tmp2 --make-bed --out wgs_data_clean

rm tmp* *.nosex
mkdir tmp
mv sc_data.* tmp
mv wgs_data.* tmp
mv snps.txt tmp
mv AT_CG.txt tmp
mv *log tmp
mv *hg*bed tmp/

module load plink2
PLINK19="/project/gazal_569/soft/plink1.9/plink"

#loop per cell type
mkdir merged/
for CELLTYPE in Breast_Mammary_Tissue Esophagus_Mucosa Esophagus_Muscularis Heart_Left_Ventricle Lung Muscle_Skeletal Prostate Skin_Sun_Exposed_Lower_leg; do
    grep $CELLTYPE sc_data_clean.fam > tmp1.list
    cut -f1 tmp1.list | sed s/\\./\\t/g | cut -f1 | awk '{print 0, "GTEX-"$1}'  > tmp2.list
    $PLINK19 --bfile sc_data_clean --keep tmp1.list --geno 0.1 --make-bed --out tmp1
    $PLINK19 --bfile wgs_data_clean --keep tmp2.list --extract tmp1.bim --make-bed --out tmp2
    $PLINK19 --bfile tmp1 --bmerge tmp2 --make-bed --out $CELLTYPE
    rm tmp?.*
    #keep UTRs SNPs
    awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' $CELLTYPE.bim > $CELLTYPE.bim.bed
    module load bedtools2
    bedtools intersect -a $CELLTYPE.bim.bed -b exon_UTRS.hg38.bed | cut -f4 > $CELLTYPE.exon_UTRS.list
    #
    $PLINK19 --bfile $CELLTYPE --recode A --out merged/$CELLTYPE
    $PLINK19 --bfile $CELLTYPE --recode A --extract $CELLTYPE.exon_UTRS.list --out merged/$CELLTYPE.exon_UTRS
    #
    $PLINK19 --bfile $CELLTYPE --make-bed --extract $CELLTYPE.exon_UTRS.list --out $CELLTYPE.exon_UTRS
    plink2 --bfile $CELLTYPE.exon_UTRS --recode A --extract-intersect bed0 ../00_TGP/Results/TGP_PLINK_HG38_maf05.bim.bed --out merged/$CELLTYPE.exon_UTRS.maf05
    plink2 --bfile $CELLTYPE.exon_UTRS --recode A --extract-intersect bed0 ../00_TGP/Results/TGP_PLINK_HG38_maf02.bim.bed --out merged/$CELLTYPE.exon_UTRS.maf02
    plink2 --bfile $CELLTYPE.exon_UTRS --recode A --extract-intersect bed0 ../00_TGP/Results/TGP_PLINK_HG38_maf01.bim.bed --out merged/$CELLTYPE.exon_UTRS.maf01
    rm $CELLTYPE.*
done
rm merged/*log merged/*nosex

R

get_error_rate <- function(celltype,mytype=".exon_UTRS"){
#get_error_rate <- function(celltype,mytype=""){

    data = read.table(paste0("merged/",celltype,mytype,".raw"),h=T)[,-(1:6)]
    
    if (celltype=="Esophagus_Mucosa"){
        data = data[c(1,2,3,6,5,4),]
    }

    x=matrix(0,3,3)
    nbind = nrow(data)/2
    for (i in 1:nbind){
        x = x + table(as.integer(data[i,]),as.integer(data[i+nbind,]))[1:3,1:3]
    }
    error_rate = sum(diag(x))/sum(x)
    #x
    #error_rate
    #i=1; x[i,i]/sum(x[i,])
    #i=2; x[i,i]/sum(x[i,])
    #i=3; x[i,i]/sum(x[i,])
    out = NULL
    out$error_rate = 1-error_rate
    out$x = x
    out$snps = ncol(data)
    out
}

get_all_error_rate  <- function(celltype){
    res1 = get_error_rate(celltype,"")
    res2 = get_error_rate(celltype,".exon_UTRS")
    c(res1$snps,res2$snps,res1$error_rate,res2$error_rate)
}

all_tissues = c(get_error_rate("Breast_Mammary_Tissue")$error_rate, get_error_rate("Esophagus_Mucosa")$error_rate, get_error_rate("Esophagus_Muscularis")$error_rate, get_error_rate("Heart_Left_Ventricle")$error_rate, get_error_rate("Lung")$error_rate, get_error_rate("Muscle_Skeletal")$error_rate, get_error_rate("Prostate")$error_rate, get_error_rate("Skin_Sun_Exposed_Lower_leg")$error_rate)


Breast_Mammary_Tissue      = get_all_error_rate("Breast_Mammary_Tissue")
Esophagus_Mucosa           = get_all_error_rate("Esophagus_Mucosa")
Esophagus_Muscularis       = get_all_error_rate("Esophagus_Muscularis")
Heart_Left_Ventricle       = get_all_error_rate("Heart_Left_Ventricle")
Lung                       = get_all_error_rate("Lung")
Muscle_Skeletal            = get_all_error_rate("Muscle_Skeletal")
Prostate                   = get_all_error_rate("Prostate")
Skin_Sun_Exposed_Lower_leg = get_all_error_rate("Skin_Sun_Exposed_Lower_leg")

Breast_Mammary_Tissue
Esophagus_Mucosa
Esophagus_Muscularis
Heart_Left_Ventricle
Lung
Muscle_Skeletal
Prostate
Skin_Sun_Exposed_Lower_leg


#error rate for AC AG CT GT
pairs     = c("A_C", "A_G", "C_T", "G_T")
pairs_inv = c("C_A", "G_A", "T_C", "T_G")

get_error_rate_perSNPtype <- function(celltype,mytype=".exon_UTRS"){

    data = read.table(paste0("merged/",celltype,mytype,".raw"),h=T)[,-(1:6)]
    info = matrix(unlist(strsplit(colnames(data), "\\.")), ncol = 4, byrow = TRUE)
    
    if (celltype=="Esophagus_Mucosa"){
        data = data[c(1,2,3,6,5,4),]
    }

    error_rate = NULL
    for (i in 1:4){
        tokeep = c(which(info[,4] %in% pairs[i]),which(info[,4] %in% pairs_inv[i]))
        mydata = data[,tokeep]
        x=matrix(0,3,3)
        nbind = nrow(mydata)/2
        for (i in 1:nbind){
            x = x + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
        }
        error_rate = c(error_rate, 1 - sum(diag(x))/sum(x))
    }
    names(error_rate) = pairs
    error_rate
}
Breast_Mammary_Tissue      = get_error_rate_perSNPtype("Breast_Mammary_Tissue")
Esophagus_Mucosa           = get_error_rate_perSNPtype("Esophagus_Mucosa")
Esophagus_Muscularis       = get_error_rate_perSNPtype("Esophagus_Muscularis")
Heart_Left_Ventricle       = get_error_rate_perSNPtype("Heart_Left_Ventricle")
Lung                       = get_error_rate_perSNPtype("Lung")
Muscle_Skeletal            = get_error_rate_perSNPtype("Muscle_Skeletal")
Prostate                   = get_error_rate_perSNPtype("Prostate")
Skin_Sun_Exposed_Lower_leg = get_error_rate_perSNPtype("Skin_Sun_Exposed_Lower_leg")

res = rbind(Breast_Mammary_Tissue,
Esophagus_Mucosa,
Esophagus_Muscularis,
Heart_Left_Ventricle,
Lung,
Muscle_Skeletal,
Prostate,
Skin_Sun_Exposed_Lower_leg)

write.csv(res, file = "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Plot_Table/Genotype_error/genotype_error_pair.csv", row.names = TRUE)

#
get_error_rate_perallele <- function(celltype,mytype=".exon_UTRS"){

    data = read.table(paste0("merged/",celltype,mytype,".raw"),h=T)[,-(1:6)]
    info = matrix(unlist(strsplit(colnames(data), "\\.")), ncol = 4, byrow = TRUE)

    if (celltype=="Esophagus_Mucosa"){
        data = data[c(1,2,3,6,5,4),]
    }
    
    #A_to_C and C_to_A
    tokeep = which(info[,4] %in% "A_C") # counting the number of C
    mydata = data[,tokeep]
    x1=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x1 = x1 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    tokeep = which(info[,4] %in% "C_A") #counting the number of A
    mydata = data[,tokeep]
    x2=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x2 = x2 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    error_rate_A_to_C = (x1[2,1]+2*x1[3,1]+x1[3,2]+x2[1,2]+2*x2[1,3]+x2[2,3])/(sum(x1[,1])*2+sum(x1[,2])+sum(x2[,3])*2+sum(x2[,2]))
    error_rate_C_to_A = (x1[1,2]+2*x1[1,3]+x1[2,3]+x2[2,1]+2*x2[3,1]+x2[3,2])/(sum(x1[,3])*2+sum(x1[,2])+sum(x2[,1])*2+sum(x2[,2]))

    #A_to_G and G_to_A
    tokeep = which(info[,4] %in% "A_G") # counting the number of G
    mydata = data[,tokeep]
    x1=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x1 = x1 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    tokeep = which(info[,4] %in% "G_A") #counting the number of A
    mydata = data[,tokeep]
    x2=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x2 = x2 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    error_rate_A_to_G = (x1[2,1]+2*x1[3,1]+x1[3,2]+x2[1,2]+2*x2[1,3]+x2[2,3])/(sum(x1[,1])*2+sum(x1[,2])+sum(x2[,3])*2+sum(x2[,2]))
    error_rate_G_to_A = (x1[1,2]+2*x1[1,3]+x1[2,3]+x2[2,1]+2*x2[3,1]+x2[3,2])/(sum(x1[,3])*2+sum(x1[,2])+sum(x2[,1])*2+sum(x2[,2]))
    
    #C_to_T and T_to_C
    tokeep = which(info[,4] %in% "C_T") # counting the number of T
    mydata = data[,tokeep]
    x1=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x1 = x1 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    tokeep = which(info[,4] %in% "T_C") #counting the number of C
    mydata = data[,tokeep]
    x2=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x2 = x2 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    error_rate_C_to_T = (x1[2,1]+2*x1[3,1]+x1[3,2]+x2[1,2]+2*x2[1,3]+x2[2,3])/(sum(x1[,1])*2+sum(x1[,2])+sum(x2[,3])*2+sum(x2[,2]))
    error_rate_T_to_C = (x1[1,2]+2*x1[1,3]+x1[2,3]+x2[2,1]+2*x2[3,1]+x2[3,2])/(sum(x1[,3])*2+sum(x1[,2])+sum(x2[,1])*2+sum(x2[,2]))

    #G_to_T and T_to_G
    tokeep = which(info[,4] %in% "G_T") # counting the number of T
    mydata = data[,tokeep]
    x1=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x1 = x1 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    tokeep = which(info[,4] %in% "T_G") #counting the number of G
    mydata = data[,tokeep]
    x2=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in 1:nbind){
        x2 = x2 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    #
    error_rate_G_to_T = (x1[2,1]+2*x1[3,1]+x1[3,2]+x2[1,2]+2*x2[1,3]+x2[2,3])/(sum(x1[,1])*2+sum(x1[,2])+sum(x2[,3])*2+sum(x2[,2]))
    error_rate_T_to_G = (x1[1,2]+2*x1[1,3]+x1[2,3]+x2[2,1]+2*x2[3,1]+x2[3,2])/(sum(x1[,3])*2+sum(x1[,2])+sum(x2[,1])*2+sum(x2[,2]))

    error_rate = c(error_rate_A_to_C,error_rate_A_to_G,error_rate_C_to_A,error_rate_C_to_T,error_rate_G_to_A,error_rate_G_to_T,error_rate_T_to_C,error_rate_T_to_G)
    names(error_rate) = c("A_to_C","A_to_G","C_to_A","C_to_T","G_to_A","G_to_T","T_to_C","T_to_G")
    error_rate
}

Breast_Mammary_Tissue      = get_error_rate_perallele("Breast_Mammary_Tissue")
Esophagus_Mucosa           = get_error_rate_perallele("Esophagus_Mucosa")
Esophagus_Muscularis       = get_error_rate_perallele("Esophagus_Muscularis")
Heart_Left_Ventricle       = get_error_rate_perallele("Heart_Left_Ventricle")
Lung                       = get_error_rate_perallele("Lung")
Muscle_Skeletal            = get_error_rate_perallele("Muscle_Skeletal")
Prostate                   = get_error_rate_perallele("Prostate")
Skin_Sun_Exposed_Lower_leg = get_error_rate_perallele("Skin_Sun_Exposed_Lower_leg")

res = rbind(Breast_Mammary_Tissue,
Esophagus_Mucosa,
Esophagus_Muscularis,
Heart_Left_Ventricle,
Lung,
Muscle_Skeletal,
Prostate,
Skin_Sun_Exposed_Lower_leg)

write.csv(res, file = "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Plot_Table/Genotype_error/genotype_error_allele.csv", row.names = TRUE)

