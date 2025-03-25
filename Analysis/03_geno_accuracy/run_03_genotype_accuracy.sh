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
cp /project/gazal_569/jianing/Archive/HCA/ancestry/exon_UTRS.bed .
awk -v OFS="\t" '{ print "chr"$1, $4-1, $4, $2, $1, $5, $6 }' merged.bim > merged.bim.bed
module load bedtools2
bedtools intersect -a merged.bim.bed -b exon_UTRS.bed | cut -f4 > merged.exon_UTRS.list

$PLINK2 --bfile merged --recode A --extract merged.exon_UTRS.list --out merged
$PLINK2 --bfile merged --recode transpose --extract merged.exon_UTRS.list --out merged

