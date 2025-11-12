$CHR=$ARGV[0];

open(IN ,"gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.bim") || die;
open(OUT1,">gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.new.bim") || die;
open(OUT2,">gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.torm.bim") || die;
while (<IN>) {
	chomp $_; @line=split;
    #
    if ($line[1] eq ".") {$line[1] = "$line[0]_$line[3]_$line[4]_$line[5]"}
    printf OUT1 "@line\n";
    #
	if(defined($rs{$line[1]})) {
		printf OUT2 "$line[1]\n";
	} else {
		$rs{$line[1]}=1;
	}
}
close IN;
close OUT;

system("/project/gazal_569/soft/plink1.9/plink --bed gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.bed --bim gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.new.bim --fam gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.fam --geno 0.01 --exclude gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.torm.bim --make-bed --out gnomad.genomes.v3.1.hgdp_1kg_subset.qc.chr$CHR");

system("rm -rf gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.new.bim gnomad.genomes.v3.1.hgdp_1kg_subset.chr$CHR.torm.bim gnomad.genomes.v3.1.hgdp_1kg_subset.qc.chr$CHR.nosex");

