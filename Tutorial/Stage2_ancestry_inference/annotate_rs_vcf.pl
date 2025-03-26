$plink_file = $ARGV[0]; 
$pca_file   = $ARGV[1]; 

#
print "1/ Stock SNPs from $plink_file.bim \n";
open(IN,"$plink_file.bim") || die;
$cpt_in  = 0;
while (<IN>) {
	chomp $_; @line=split;
    $rs{$line[0]}{$line[3]}{$line[4]}{$line[5]} = "toremove";
    $cpt_in++;
}
close IN;
print "   $cpt_in SNPs in $plink_file.bim.\n\n";

print "2/ Stock rs ID from $pca_file.bed\n";
open(IN,"$pca_file.bed") || die;
$cpt=0;
while (<IN>) {
	chomp $_; @line=split;
    if (defined($rs{$line[4]}{$line[2]}{$line[5]}{$line[6]})){ $rs{$line[4]}{$line[2]}{$line[5]}{$line[6]} = $line[3]; }
    if (defined($rs{$line[4]}{$line[2]}{$line[6]}{$line[5]})){ $rs{$line[4]}{$line[2]}{$line[6]}{$line[5]} = $line[3]; }
    $cpt++;
}
close IN;
print "   $cpt SNPs in $pca_file.bed.\n\n";

print "3/ Update $plink_file.bim with rs id\n";
system "cp $plink_file.bim $plink_file.tmp.bim";
open(IN,"$plink_file.tmp.bim") || die;
open(OUT,">$plink_file.bim") || die;
while (<IN>) {
	chomp $_; @line=split;
    print OUT "$line[0]\t$rs{$line[0]}{$line[3]}{$line[4]}{$line[5]}\t0\t$line[3]\t$line[4]\t$line[5]\n";
}
close IN;
close OUT;
system "rm $plink_file.tmp.bim";

