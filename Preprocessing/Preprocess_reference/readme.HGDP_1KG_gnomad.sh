for CHR in {1..22}; do
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.vcf.bgz
done
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.hgdp_1kg_subset.sample_meta.tsv.gz
gunzip gnomad.genomes.v3.1.hgdp_1kg_subset.sample_meta.tsv.gz


mkdir plink_maf01
module load gcc/8.3.0  openblas/0.3.8 plink2/2.00a2.3-openblas
for CHR in {1..22}; do
    sbatch --mem=50000 -t 0-10:00 --wrap="plink2 --vcf gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.vcf.bgz --maf 0.01 --make-bed --out plink_maf01/gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR"
done
cd plink_maf01
for chr in {1..22}; do
    perl qc.pl $chr
done
for chr in {2..22}; do
    echo "gnomad.genomes.v3.1.hgdp_1kg_subset.qc.chr$chr.bed gnomad.genomes.v3.1.hgdp_1kg_subset.qc.chr$chr.bim gnomad.genomes.v3.1.hgdp_1kg_subset.qc.chr$chr.fam" >> list.bed
done
PLINK2="/project/gazal_569/soft/plink1.9/plink"
$PLINK2 --bfile gnomad.genomes.v3.1.hgdp_1kg_subset.qc.chr1  --merge-list list.bed --make-bed --allow-no-sex --out gnomad.genomes.v3.1.hgdp_1kg_subset.qc
rm gnomad.genomes.v3.1.hgdp_1kg_subset*chr* *log *nosex list.bed

##### lift over to hg19

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
chmod 777 liftOver
mkdir unlifted

awk '{print "chr"$1,$4-1,$4,$2}' gnomad.genomes.v3.1.hgdp_1kg_subset.qc.bim > hg38.bed
./liftOver \
    hg38.bed \
    hg38ToHg19.over.chain.gz \
    hg19.bed \
    hg19.unlifted.bed
grep -v new hg19.unlifted.bed | awk '{print $4}' > hg19.rs_to_remove.txt

$PLINK2 \
--bfile gnomad.genomes.v3.1.hgdp_1kg_subset.qc \
--exclude hg19.rs_to_remove.txt \
--make-bed \
--out gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19

cut -f1-3 gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.bim > temp1
cut -f3   hg19.bed                                   > temp2
cut -f5-  gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.bim > temp3
paste temp1 temp2 temp3 > gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.bim

rm hg*.* temp* *log *nosex *.over.chain.gz  liftOver 

awk '{print $1 "\t0\t" $2}' ../../gene_info/hg19/hg19.fa.fai > human_hg19.bed
awk '{print "chr"$1"\t"$4-1"\t"$4}' gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.bim | sort -k1,1 -k2,2n | uniq > tmp
bedtools merge -i tmp > tmp2
bedtools intersect -a tmp2 -b human_hg19.bed > gnomad.genomes.v3.1.hgdp_1kg_subset.qc_hg19.interval.bed
rm tmp* human_hg19.bed

