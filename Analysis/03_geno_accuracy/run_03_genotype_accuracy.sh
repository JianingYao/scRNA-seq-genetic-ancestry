PLINK2="/project/gazal_569/soft/plink1.9/plink"
VCF="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/VCF/eQTLAutoimmune.SNP.Filtered.SV.vcf.gz"
ORI="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/raw_genotype/OneK1K_AllChr"

#get sc data
$PLINK2 --vcf $VCF --geno 0.10 --make-bed --out sc_data
#edit bim file
awk '{print $1, $1":"$4, $3, $4, $5, $6}' sc_data.bim > tmp
mv tmp sc_data.bim

#get geno data
cp $ORI.bed geno_data.bed
cp $ORI.bim geno_data.bim
cp $ORI.fam geno_data.fam
awk '{print $1, $1":"$4, $3, $4, $5, $6}' geno_data.bim > tmp
mv tmp geno_data.bim

#SNPs in both files
R
data1 = read.table("sc_data.bim")[,2]
data2 = read.table("geno_data.bim")[,2]
inter = intersect(data1,data2)
write.table(inter,file="snps.txt",quote=F,col.names=F,row.names=F)
q()
n

#IID to keep
awk '{print $1, $2}' sc_data.fam > iids.txt

#extract snps and iids from the two files
$PLINK2 --bfile sc_data   --extract snps.txt                 --make-bed --out sc_data.filtered
$PLINK2 --bfile geno_data --extract snps.txt --keep iids.txt --make-bed --out geno_data.filtered
#rename fam info
awk '{print "sc", $2, $3, $4, $5, $6}' sc_data.filtered.fam > tmp
mv tmp sc_data.filtered.fam
awk '{print "geno", $2, $3, $4, 0, $6}' geno_data.filtered.fam > tmp
mv tmp geno_data.filtered.fam

#merge the files 
cp /project/gazal_569/DATA/yazar2022/genotype/geno_n981_plink/allchr.* .
sed -e s/\_/\\t/g allchr.fam | awk '{print "geno", $3,$4,$5,$6,$7}' > tmp
mv tmp allchr.fam
$PLINK2 --bfile allchr --extract geno_data.filtered.bim --keep geno_data.filtered.fam --make-bed --out geno_data.filtered
$PLINK2 --bfile sc_data.filtered --bmerge geno_data.filtered --make-bed --out merged
#$PLINK2 --bfile sc_data.filtered --make-bed --flip merged-merge.missnp --out sc_data.filtered.flip
#$PLINK2 --bfile sc_data.filtered.flip --bmerge geno_data.filtered --make-bed --out merged
#$PLINK2 --bfile sc_data.filtered.flip --exclude merged-merge.missnp --make-bed --out sc_data.qc
$PLINK2 --bfile sc_data.filtered      --exclude merged-merge.missnp --make-bed --out sc_data.qc
$PLINK2 --bfile geno_data.filtered    --exclude merged-merge.missnp --make-bed --out geno_data.qc
#create list of AT/CG variants
grep C sc_data.qc.bim | grep G | cut -f2 > AT_CG.txt
grep A sc_data.qc.bim | grep T | cut -f2 >> AT_CG.txt
$PLINK2 --bfile sc_data.qc --bmerge geno_data.qc --exclude AT_CG.txt --make-bed --out merged
rm *filtered* *log *.nosex AT_CG.txt

#keep UTRs SNPs
cp /project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/Results/exon_UTRS.bed .
awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' merged.bim > merged.bim.bed
module load bedtools2
bedtools intersect -a merged.bim.bed -b exon_UTRS.bed | cut -f4 > merged.exon_UTRS.list

$PLINK2 --bfile merged --recode A --extract merged.exon_UTRS.list --out merged
#$PLINK2 --bfile merged --recode transpose --extract merged.exon_UTRS.list --out merged      #within exon_UTRS
$PLINK2 --bfile merged --recode A --out merged_all      #all


$PLINK2 --bfile merged --make-bed --extract merged.exon_UTRS.list --out merged.exon_UTRS
awk '{ print $1, $1":"$4":"$5":"$6, $3, $4, $5, $6}' merged.exon_UTRS.bim > tmp.bim
mv tmp.bim merged.exon_UTRS.bim
$PLINK2 --bfile merged.exon_UTRS --recode A --out merged.exon_UTRS

R


get_error_rate <- function(mytype){

    data = read.table(paste0("merged",mytype,".raw"),h=T)[,-(1:6)]
    
    x=matrix(0,3,3)
    nbind = nrow(data)/2
    for (i in 1:nbind){
        x = x + table(as.integer(data[i,]),as.integer(data[i+nbind,]))[1:3,1:3]
    }
    error_rate = sum(diag(x))/sum(x)
    out = NULL
    out$error_rate = 1-error_rate
    out$x = x
    out$snps = ncol(data)
    out
}

get_all_error_rate  <- function(){
    res1 = get_error_rate("_all")
    res2 = get_error_rate("")
    c(res1$snps,res2$snps,res1$error_rate,res2$error_rate)
}

get_all_error_rate()


#error rate for AC AG CT GT
pairs     = c("A_C", "A_G", "C_T", "G_T")
pairs_inv = c("C_A", "G_A", "T_C", "T_G")

get_error_rate_perSNPtype <- function(){

    data = read.table("merged.exon_UTRS.raw",h=T)[,-(1:6)]
    info = matrix(unlist(strsplit(colnames(data), "\\.")), ncol = 4, byrow = TRUE)
    
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
get_error_rate_perSNPtype()

#
get_error_rate_perallele <- function(){

    data = read.table("merged.exon_UTRS.raw",h=T)[,-(1:6)]
    info = matrix(unlist(strsplit(colnames(data), "\\.")), ncol = 4, byrow = TRUE)
    
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
    for (i in (1:nbind)[-c(3,7,8,12)]){
        x1 = x1 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    for (i in c(3,7,8,12)){
        x = table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))
        x1[1:nrow(x),1:ncol(x)] = x1[1:nrow(x),1:ncol(x)] + x
    }
    #
    tokeep = which(info[,4] %in% "T_G") #counting the number of G
    mydata = data[,tokeep]
    x2=matrix(0,3,3)
    nbind = nrow(mydata)/2
    for (i in (1:nbind)[-c(4,9,10)]){
        x2 = x2 + table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))[1:3,1:3]
    }
    for (i in c(4,9,10)){
        x = table(as.integer(mydata[i,]),as.integer(mydata[i+nbind,]))
        x2[1:nrow(x),1:ncol(x)] = x2[1:nrow(x),1:ncol(x)] + x
    }
    #
    error_rate_G_to_T = (x1[2,1]+2*x1[3,1]+x1[3,2]+x2[1,2]+2*x2[1,3]+x2[2,3])/(sum(x1[,1])*2+sum(x1[,2])+sum(x2[,3])*2+sum(x2[,2]))
    error_rate_T_to_G = (x1[1,2]+2*x1[1,3]+x1[2,3]+x2[2,1]+2*x2[3,1]+x2[3,2])/(sum(x1[,3])*2+sum(x1[,2])+sum(x2[,1])*2+sum(x2[,2]))

    error_rate = c(error_rate_A_to_C,error_rate_A_to_G,error_rate_C_to_A,error_rate_C_to_T,error_rate_G_to_A,error_rate_G_to_T,error_rate_T_to_C,error_rate_T_to_G)
    names(error_rate) = c("A_to_C","A_to_G","C_to_A","C_to_T","G_to_A","G_to_T","T_to_C","T_to_G")
    error_rate
}

get_error_rate_perallele()

data = read.table("merged.raw",h=T)[,-(1:6)]

x=matrix(0,3,3)
for (i in 1:12){
    x = x + table(as.integer(data[i,]),as.integer(data[i+12,]))[1:3,1:3]
}
x
error_rate = sum(diag(x))/sum(x)
i=1; x[i,i]/sum(x[i,])
i=2; x[i,i]/sum(x[i,])
i=3; x[i,i]/sum(x[i,])

# table(apply(data[1:12,],2,sum)-apply(data[13:24,],2,sum))
# which(abs(apply(data[1:12,],2,sum)-apply(data[13:24,],2,sum))>=10)

error_rate = 0.9222302
# or
error_rate = rbind(
c(8733 , 201 , 10),
c( 481 ,3475 ,217),
c(   6 , 166 ,611))
for (i in 1:3){error_rate[i,]=error_rate[i,]/sum(error_rate[i,])}

# individual results
out = NULL
for (i in 1:12){
  tmp = c(ncol(data)-sum(is.na(data[i,]+data[12+i,])),sum(data[i,]==data[12+i,],na.rm=T))
  out = rbind(out, c(tmp,(tmp[1]-tmp[2])/tmp[1]))
}
tmp = apply(out,2,sum)
out = rbind(out, c(tmp[1:2],(tmp[1]-tmp[2])/tmp[1]))
out
