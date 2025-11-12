# Tutorial on genetic-ancestry inference using scRNA data
Jianing Yao

- [<span class="toc-section-number">0.1</span> QC](#qc)
- [<span class="toc-section-number">0.2</span> PCA
  analysis](#pca-analysis)
- [<span class="toc-section-number">0.3</span> ADMIXTURE
  analysis](#admixture-analysis)

#### This tutorial will showcase genetic-ancestry inference on a specific HCA study after obtaining SNPs from the scRNA data following steps described in `Stage1_preprocessing`.

### QC

``` bash
# Set up
PLINK2="../../bin/plink1.9/plink"
PCA_FILE="HGDP_1kGP.exon_UTRS"
tar -xJvf ../Stage2_ancestry_inference/$PCA_FILE.tar.xz
STUDY='StemCellsInChronicMyeloidLeukemia'
VCF='StemCellsInChronicMyeloidLeukemia.SNP.Filtered.SV.vcf.gz'
DIR=$STUDY
if [ ! -d "$DIR" ]; then
    mkdir -p "$DIR"
fi

#### Step 1: Convert your vcf into plink format and QC the data
$PLINK2 --const-fid 0 --vcf $VCF --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step1
$PLINK2 --bfile $DIR/$STUDY.plink.step1 --mind 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink.step2
$PLINK2 --const-fid 0 --vcf $VCF --keep $DIR/$STUDY.plink.step2.fam --geno 0.10 --make-bed --allow-extra-chr --out $DIR/$STUDY.plink
rm $DIR/$STUDY.plink.step*

#### Step 2: Find rs ID and create a list with rs id in both scRNA dataset and the reference dataset
perl ../Stage2_ancestry_inference/annotate_rs_vcf.pl $DIR/$STUDY.plink $PCA_FILE
grep -v toremove $DIR/$STUDY.plink.bim | cut -f2 > $DIR/$STUDY.plink.list


#### Step 3: Merge the 2 files and perform pruning
$PLINK2 --bfile $DIR/$STUDY.plink --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp1 --allow-extra-chr
$PLINK2 --bfile $PCA_FILE         --extract $DIR/$STUDY.plink.list --make-bed --out $DIR/$STUDY.tmp2
$PLINK2 --bfile $DIR/$STUDY.tmp1  --bmerge $DIR/$STUDY.tmp2 --make-bed --allow-no-sex --out $DIR/$STUDY.HGDP_1kGP
$PLINK2 --bfile $DIR/$STUDY.HGDP_1kGP   --indep-pairwise 50 10 0.1 --out $DIR/$STUDY.HGDP_1kGP.pca
$PLINK2 --bfile $DIR/$STUDY.HGDP_1kGP   --extract $DIR/$STUDY.HGDP_1kGP.pca.prune.in --make-bed --out $DIR/$STUDY.HGDP_1kGP.pca
rm $DIR/$STUDY.tmp*


#### Step 4: Run PCA
awk '{print $1, $2, "MYID"}' $DIR/$STUDY.plink.fam > $DIR/$STUDY.HGDP_1kGP.cluster
awk '{print $1, $2, "HGDP_1kGP"}' $PCA_FILE.fam >> $DIR/$STUDY.HGDP_1kGP.cluster
$PLINK2 --bfile $DIR/$STUDY.HGDP_1kGP.pca --pca --within $DIR/$STUDY.HGDP_1kGP.cluster --pca-cluster-names HGDP_1kGP --out $DIR/$STUDY.HGDP_1kGP.pca
```

### PCA analysis

``` r
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(readr)

STUDY = "StemCellsInChronicMyeloidLeukemia"
data=read.table(paste0(STUDY,"/",STUDY,".HGDP_1kGP.pca.eigenvec"),h=F)
rownames(data) = data$V2
info=read.table("info.txt",h=T,sep="\t")
rownames(info) = info$IID
mydata  = data[as.character(read.table(paste0(STUDY,"/",STUDY,".plink.fam"))[,2]),]
data_HGDP_1kGP = data[as.character(info$IID),]

topPCtokeep = 5  #let's keep the top 5 PCs
tokeep = topPCtokeep+2

#compute center of each reference population
center_afr = apply(subset(data_HGDP_1kGP[,3:tokeep],info$myREG=="Africa"),2,mean)
center_amr = apply(subset(data_HGDP_1kGP[,3:tokeep],info$myREG=="American"),2,mean)
center_eas = apply(subset(data_HGDP_1kGP[,3:tokeep],info$myREG=="East Asia"),2,mean)
center_eur = apply(subset(data_HGDP_1kGP[,3:tokeep],info$myREG=="European" | info$myREG=="Finnish"),2,mean)
center_mea = apply(subset(data_HGDP_1kGP[,3:tokeep],info$myREG=="Middle-East"),2,mean)
center_sas = apply(subset(data_HGDP_1kGP[,3:tokeep],info$myREG=="South Asia"),2,mean)
out=NULL
for (i in 1:nrow(mydata)){
    thisiid = c(
    sum((mydata[i,3:tokeep]-center_afr)**2),
    sum((mydata[i,3:tokeep]-center_amr)**2),
    sum((mydata[i,3:tokeep]-center_eas)**2),
    sum((mydata[i,3:tokeep]-center_eur)**2),
    sum((mydata[i,3:tokeep]-center_mea)**2),
    sum((mydata[i,3:tokeep]-center_sas)**2))
    thisiid = as.integer(thisiid == min(thisiid))
    out=rbind(out,thisiid)
}
colnames(out) = c("afr","amr","eas","eur","mea","sas")
rownames(out)= mydata[,2]
write.table(out,file=paste0(STUDY,"/",STUDY,".pca_center.txt"),col.names=T,row.names=F,quote=F,sep="\t")

mylegend = c("Africa", "America", "East Asia", "Europe", "Middle East", "South Asia")
mycol = c("#1B9E77","#D95F02","#7570B3","#E6AB02","#66A61E","#E7298A")
populations <- c("Africa", "American", "East Asia", "European", "Middle-East", "South Asia")
info$Color <- mycol[match(info$myREG, populations)]
#
x1=c(data_HGDP_1kGP$V3,mydata$V3); y1=c(data_HGDP_1kGP$V4,mydata$V4)
x2=c(data_HGDP_1kGP$V5,mydata$V5); y2=c(data_HGDP_1kGP$V6,mydata$V6)

df_ref <- data.frame(
  PC1 = data_HGDP_1kGP$V3, PC2 = data_HGDP_1kGP$V4,
  PC3 = data_HGDP_1kGP$V5, PC4 = data_HGDP_1kGP$V6,
  REG = factor(info$myREG, levels = populations)
)

df_my <- data.frame(
  PC1 = mydata$V3, PC2 = mydata$V4,
  PC3 = mydata$V5, PC4 = mydata$V6
)

pop_cols <- setNames(mycol, populations)

p12 <- ggplot() +
  geom_point(data = df_ref, aes(PC1, PC2, color = REG), size = 1.6, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = df_my,  aes(PC1, PC2), color = "black", size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = pop_cols, breaks = populations, labels = mylegend, drop = FALSE) +
  coord_cartesian(xlim = range(x1, na.rm = TRUE), ylim = range(y1, na.rm = TRUE), expand = TRUE) +
  labs(x = "PC1", y = "PC2", color = "Population") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

p34 <- ggplot() +
  geom_point(data = df_ref, aes(PC3, PC4, color = REG), size = 1.6, alpha = 0.7, na.rm = TRUE) +
  geom_point(data = df_my,  aes(PC3, PC4), color = "black", size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = pop_cols, breaks = populations, labels = mylegend, drop = FALSE) +
  coord_cartesian(xlim = range(x2, na.rm = TRUE), ylim = range(y2, na.rm = TRUE), expand = TRUE) +
  labs(x = "PC3", y = "PC4", color = "Population") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

legend_row <- guides(color = guide_legend(nrow = 1, byrow = TRUE))

p12 <- p12 + legend_row + theme(legend.position = "top")
p34 <- p34 + legend_row + theme(legend.position = "top")

(p12 + p34 + plot_layout(guides = "collect", widths = c(1, 1))) &
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box.just = "left")
```

![](tutorial_ancestry_files/figure-commonmark/pca-1.png)

### ADMIXTURE analysis

``` r
STUDY = "StemCellsInChronicMyeloidLeukemia"

bed=read.table(paste0(STUDY,"/",STUDY,".HGDP_1kGP.fam"),h=F)[,1:2]
info=read.table("info.txt",h=T,sep="\t")
pop=rep(".",nrow(bed)); names(pop)=bed$V2

for (mypop in c("eas","amr","sas","afr","mid","eur")){
    pop[as.character(info$IID[which(info$REG==mypop)])] = mypop
}

write.table(pop,file=paste0(STUDY,"/",STUDY,".HGDP_1kGP.pca.pop"),sep="\t",quote=F,col.names=F,row.names=F)
```

``` bash
STUDY='StemCellsInChronicMyeloidLeukemia'
DIR=$STUDY
cd $DIR
../../../bin/admixture_linux-1.3.0/admixture --supervised $STUDY.HGDP_1kGP.pca.bed 6
```

    ****                   ADMIXTURE Version 1.3.0                  ****
    ****                    Copyright 2008-2015                     ****
    ****           David Alexander, Suyash Shringarpure,            ****
    ****                John  Novembre, Ken Lange                   ****
    ****                                                            ****
    ****                 Please cite our paper!                     ****
    ****   Information at www.genetics.ucla.edu/software/admixture  ****

    Random seed: 43
    Point estimation method: Block relaxation algorithm
    Convergence acceleration algorithm: QuasiNewton, 3 secant conditions
    Point estimation will terminate when objective function delta < 0.0001
    Estimation of standard errors disabled; will compute point estimates only.
    Supervised analysis mode.  Examining .pop file...
    Size of G: 3518x2968
    Performing five EM steps to prime main algorithm
    1 (EM)  Elapsed: 0.679  Loglikelihood: -9.42888e+06 (delta): 1.15548e+07
    2 (EM)  Elapsed: 0.679  Loglikelihood: -9.42747e+06 (delta): 1413.29
    3 (EM)  Elapsed: 0.679  Loglikelihood: -9.42721e+06 (delta): 253.83
    4 (EM)  Elapsed: 0.679  Loglikelihood: -9.42697e+06 (delta): 237.982
    5 (EM)  Elapsed: 0.68   Loglikelihood: -9.42675e+06 (delta): 224.917
    Initial loglikelihood: -9.42675e+06
    Starting main algorithm
    1 (QN/Block)    Elapsed: 0.996  Loglikelihood: -9.42174e+06 (delta): 5009.35
    2 (QN/Block)    Elapsed: 1.027  Loglikelihood: -9.4217e+06  (delta): 40.3882
    3 (QN/Block)    Elapsed: 2.37   Loglikelihood: -9.42169e+06 (delta): 6.80939
    4 (QN/Block)    Elapsed: 1.604  Loglikelihood: -9.42169e+06 (delta): 4.30699
    5 (QN/Block)    Elapsed: 1.602  Loglikelihood: -9.42169e+06 (delta): 1.1896
    6 (QN/Block)    Elapsed: 1.603  Loglikelihood: -9.42169e+06 (delta): 0.0206683
    7 (QN/Block)    Elapsed: 1.603  Loglikelihood: -9.42169e+06 (delta): 9.36352e-06
    Summary: 
    Converged in 7 iterations (15.668 sec)
    Loglikelihood: -9421687.257330
    Fst divergences between estimated populations: 
        Pop0    Pop1    Pop2    Pop3    Pop4    
    Pop0    
    Pop1    0.102   
    Pop2    0.097   0.085   
    Pop3    0.032   0.066   0.081   
    Pop4    0.114   0.136   0.152   0.101   
    Pop5    0.015   0.099   0.101   0.032   0.095   
    Writing output files.

``` r
library(RColorBrewer)

STUDY = 'StemCellsInChronicMyeloidLeukemia'
info=read.table("info.txt",h=T,sep="\t")

admixture = read.table(paste0(STUDY,"/",STUDY,".HGDP_1kGP.pca.6.Q"),h=F)
admixture.pop = read.table(paste0(STUDY,"/",STUDY,".HGDP_1kGP.pca.pop"))
admixture.ind = read.table(paste0(STUDY,"/",STUDY,".HGDP_1kGP.pca.fam"))[,2]

colnames(admixture)=c("eur","eas","amr","sas","afr","mid")
rownames(admixture)=admixture.ind

myadmixture  = admixture[as.character(read.table(paste0(STUDY,"/",STUDY,".plink.fam"))[,2]),] # subset only mydata
toplot=myadmixture[,order(apply(myadmixture,2,mean),decreasing=T)]
mycol=c("#E6AB02","#7570B3","#D95F02","#E7298A","#1B9E77","#66A61E")
mycol=mycol[order(apply(myadmixture,2,mean),decreasing=T)]

for (i in 6:1){
    toplot=toplot[order(toplot[,i],decreasing=T),]
}

barplot(
  t(as.matrix(toplot)),
  col = mycol,
  xlab = "",
  ylab = "Admixture proportion",
  border = NA,
  xaxt = "n",
  ylim = c(0, 1)
)
```

![](tutorial_ancestry_files/figure-commonmark/ADM-1.png)
