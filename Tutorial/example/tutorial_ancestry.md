# Tutorial on genetic-ancestry inference using scRNA data
Jianing Yao

- [<span class="toc-section-number">0.1</span> ](#section)

### 

#### This tutorial will showcase genetic-ancestry inference on a specific HCA study after obtaining SNPs from the scRNA data following steps described in `Stage1_preprocessing`.

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

    HGDP_1kGP.exon_UTRS.bed
    HGDP_1kGP.exon_UTRS.bim
    HGDP_1kGP.exon_UTRS.fam
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1.log.
    Options in effect:
      --allow-extra-chr
      --const-fid 0
      --geno 0.10
      --make-bed
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1
      --vcf StemCellsInChronicMyeloidLeukemia.SNP.Filtered.SV.vcf.gz

    257405 MB RAM detected; reserving 128702 MB for main workspace.

    --vcf: 1k variants complete.
    --vcf: 2k variants complete.
    --vcf: 3k variants complete.
    --vcf: 4k variants complete.
    --vcf: 5k variants complete.
    --vcf: 6k variants complete.
    --vcf: 7k variants complete.
    --vcf: 8k variants complete.
    --vcf: 9k variants complete.
    --vcf: 10k variants complete.
    --vcf: 11k variants complete.
    --vcf: 12k variants complete.
    --vcf: 13k variants complete.
    --vcf: 14k variants complete.
    --vcf: 15k variants complete.
    --vcf: 16k variants complete.
    --vcf: 17k variants complete.
    --vcf: 18k variants complete.
    --vcf: 19k variants complete.
    --vcf: 20k variants complete.
    --vcf: 21k variants complete.
    --vcf: 22k variants complete.
    --vcf: 23k variants complete.
    --vcf: 24k variants complete.
    --vcf: 25k variants complete.
    --vcf: 26k variants complete.
    --vcf: 27k variants complete.
    --vcf: 28k variants complete.
    --vcf: 29k variants complete.
    --vcf: 30k variants complete.
    --vcf: 31k variants complete.
    --vcf: 32k variants complete.
    --vcf: 33k variants complete.
    --vcf: 34k variants complete.
    --vcf: 35k variants complete.
    --vcf: 36k variants complete.
    --vcf: 37k variants complete.
    --vcf: 38k variants complete.
    --vcf: 39k variants complete.
    --vcf: 40k variants complete.
    --vcf: 41k variants complete.
    --vcf: 42k variants complete.
    --vcf: 43k variants complete.
    --vcf: 44k variants complete.
    --vcf: 45k variants complete.
    --vcf: 46k variants complete.
    --vcf: 47k variants complete.
    --vcf: 48k variants complete.
    --vcf: 49k variants complete.
    --vcf: 50k variants complete.
    --vcf: 51k variants complete.
    --vcf: 52k variants complete.
    --vcf: 53k variants complete.
    --vcf: 54k variants complete.
    --vcf: 55k variants complete.
    --vcf: 56k variants complete.
    --vcf: 57k variants complete.
    --vcf: 58k variants complete.
    --vcf: 59k variants complete.
    --vcf: 60k variants complete.
    --vcf: 61k variants complete.
    --vcf: 62k variants complete.
    --vcf: 63k variants complete.
    --vcf: 64k variants complete.
    --vcf: 65k variants complete.
    --vcf: 66k variants complete.
    --vcf: 67k variants complete.
    --vcf: 68k variants complete.
    --vcf: 69k variants complete.
    --vcf: 70k variants complete.
    --vcf: 71k variants complete.
    --vcf: 72k variants complete.
    --vcf: 73k variants complete.
    --vcf: 74k variants complete.
    --vcf: 75k variants complete.
    --vcf: 76k variants complete.
    --vcf: 77k variants complete.
    --vcf: 78k variants complete.
    --vcf: 79k variants complete.
    --vcf: 80k variants complete.
    --vcf: 81k variants complete.
    --vcf: 82k variants complete.
    --vcf: 83k variants complete.
    --vcf: 84k variants complete.
    --vcf: 85k variants complete.
    --vcf: 86k variants complete.
    --vcf: 87k variants complete.
    --vcf: 88k variants complete.
    --vcf: 89k variants complete.
    --vcf: 90k variants complete.
    --vcf: 91k variants complete.
    --vcf: 92k variants complete.
    --vcf: 93k variants complete.
    --vcf: 94k variants complete.
    --vcf: 95k variants complete.
    --vcf: 96k variants complete.
    --vcf: 97k variants complete.
    --vcf: 98k variants complete.
    --vcf: 99k variants complete.
    --vcf: 100k variants complete.
    --vcf: 101k variants complete.
    --vcf: 102k variants complete.
    --vcf: 103k variants complete.
    --vcf: 104k variants complete.
    --vcf: 105k variants complete.
    --vcf: 106k variants complete.
    --vcf: 107k variants complete.
    --vcf: 108k variants complete.
    --vcf: 109k variants complete.
    --vcf: 110k variants complete.
    --vcf: 111k variants complete.
    --vcf: 112k variants complete.
    --vcf: 113k variants complete.
    --vcf: 114k variants complete.
    --vcf: 115k variants complete.
    --vcf: 116k variants complete.
    --vcf: 117k variants complete.
    --vcf: 118k variants complete.
    --vcf: 119k variants complete.
    --vcf: 120k variants complete.
    --vcf: 121k variants complete.
    --vcf: 122k variants complete.
    --vcf: 123k variants complete.
    --vcf: 124k variants complete.
    --vcf: 125k variants complete.
    --vcf: 126k variants complete.
    --vcf: 127k variants complete.
    --vcf: 128k variants complete.
    --vcf: 129k variants complete.
    --vcf: 130k variants complete.
    --vcf: 131k variants complete.
    --vcf: 132k variants complete.
    --vcf: 133k variants complete.
    --vcf: 134k variants complete.
    --vcf: 135k variants complete.
    --vcf: 136k variants complete.
    --vcf: 137k variants complete.
    --vcf: 138k variants complete.
    --vcf: 139k variants complete.
    --vcf: 140k variants complete.
    --vcf: 141k variants complete.
    --vcf: 142k variants complete.
    --vcf: 143k variants complete.
    --vcf: 144k variants complete.
    --vcf: 145k variants complete.
    --vcf: 146k variants complete.
    --vcf: 147k variants complete.
    --vcf: 148k variants complete.
    --vcf: 149k variants complete.
    --vcf: 150k variants complete.
    --vcf: 151k variants complete.
    --vcf: 152k variants complete.
    --vcf: 153k variants complete.
    --vcf: 154k variants complete.
    --vcf: 155k variants complete.
    --vcf: 156k variants complete.
    --vcf: 157k variants complete.
    --vcf: 158k variants complete.
    --vcf: 159k variants complete.
    --vcf: 160k variants complete.
    --vcf: 161k variants complete.
    --vcf: 162k variants complete.
    --vcf: 163k variants complete.
    --vcf: 164k variants complete.
    --vcf: 165k variants complete.
    --vcf: 166k variants complete.
    --vcf: 167k variants complete.
    --vcf: 168k variants complete.
    --vcf: 169k variants complete.
    --vcf: 170k variants complete.
    --vcf: 171k variants complete.
    --vcf: 172k variants complete.
    --vcf: 173k variants complete.
    --vcf: 174k variants complete.
    --vcf: 175k variants complete.
    --vcf: 176k variants complete.
    --vcf: 177k variants complete.
    --vcf: 178k variants complete.
    --vcf: 179k variants complete.
    --vcf: 180k variants complete.
    --vcf: 181k variants complete.
    --vcf: 182k variants complete.
    --vcf: 183k variants complete.
    --vcf: 184k variants complete.
    --vcf: 185k variants complete.
    --vcf: 186k variants complete.
    --vcf: 187k variants complete.
    --vcf: 188k variants complete.
    --vcf: 189k variants complete.
    --vcf: 190k variants complete.
    --vcf: 191k variants complete.
    --vcf: 192k variants complete.
    --vcf: 193k variants complete.
    --vcf: 194k variants complete.
    --vcf: 195k variants complete.
    --vcf: 196k variants complete.
    --vcf: 197k variants complete.
    --vcf: 198k variants complete.
    --vcf: 199k variants complete.
    --vcf: 200k variants complete.
    --vcf: 201k variants complete.
    --vcf: 202k variants complete.
    --vcf: 203k variants complete.
    --vcf: 204k variants complete.
    --vcf: 205k variants complete.
    --vcf: 206k variants complete.
    --vcf: 207k variants complete.
    --vcf: 208k variants complete.
    --vcf: 209k variants complete.
    --vcf: 210k variants complete.
    --vcf: 211k variants complete.
    --vcf: 212k variants complete.
    --vcf: 213k variants complete.
    --vcf: 214k variants complete.
    --vcf: 215k variants complete.
    --vcf: 216k variants complete.
    --vcf: 217k variants complete.
    --vcf: 218k variants complete.
    --vcf: 219k variants complete.
    --vcf: 220k variants complete.
    --vcf: 221k variants complete.
    --vcf: 222k variants complete.
    --vcf: 223k variants complete.
    --vcf: 224k variants complete.
    --vcf: 225k variants complete.
    --vcf: 226k variants complete.
    --vcf: 227k variants complete.
    --vcf: 228k variants complete.
    --vcf: 229k variants complete.
    --vcf: 230k variants complete.
    --vcf: 231k variants complete.
    --vcf: 232k variants complete.
    --vcf: 233k variants complete.
    --vcf: 234k variants complete.
    --vcf: 235k variants complete.
    --vcf: 236k variants complete.
    --vcf: 237k variants complete.
    --vcf: 238k variants complete.
    --vcf: 239k variants complete.
    --vcf: 240k variants complete.
    --vcf: 241k variants complete.
    --vcf: 242k variants complete.
    --vcf: 243k variants complete.
    --vcf: 244k variants complete.
    --vcf: 245k variants complete.
    --vcf: 246k variants complete.
    --vcf: 247k variants complete.
    --vcf: 248k variants complete.
    --vcf: 249k variants complete.
    --vcf: 250k variants complete.
    --vcf: 251k variants complete.
    --vcf: 252k variants complete.
    --vcf: 253k variants complete.
    --vcf: 254k variants complete.
    --vcf: 255k variants complete.
    --vcf: 256k variants complete.
    --vcf: 257k variants complete.
    --vcf: 258k variants complete.
    --vcf: 259k variants complete.
    --vcf: 260k variants complete.
    --vcf: 261k variants complete.
    --vcf: 262k variants complete.
    --vcf: 263k variants complete.
    --vcf: 264k variants complete.
    --vcf: 265k variants complete.
    --vcf: 266k variants complete.
    --vcf: 267k variants complete.
    --vcf: 268k variants complete.
    --vcf: 269k variants complete.
    --vcf: 270k variants complete.
    --vcf: 271k variants complete.
    --vcf: 272k variants complete.
    --vcf: 273k variants complete.
    --vcf: 274k variants complete.
    --vcf: 275k variants complete.
    --vcf: 276k variants complete.
    --vcf: 277k variants complete.
    --vcf: 278k variants complete.
    --vcf: 279k variants complete.
    --vcf: 280k variants complete.
    --vcf: 281k variants complete.
    --vcf: 282k variants complete.
    --vcf: 283k variants complete.
    --vcf: 284k variants complete.
    --vcf: 285k variants complete.
    --vcf: 286k variants complete.
    --vcf: 287k variants complete.
    --vcf: 288k variants complete.
    --vcf: 289k variants complete.
    --vcf: 290k variants complete.
    --vcf: 291k variants complete.
    --vcf: 292k variants complete.
    --vcf: 293k variants complete.
    --vcf: 294k variants complete.
    --vcf: 295k variants complete.
    --vcf: 296k variants complete.
    --vcf: 297k variants complete.
    --vcf: 298k variants complete.
    --vcf: 299k variants complete.
    --vcf: 300k variants complete.
    --vcf: 301k variants complete.
    --vcf: 302k variants complete.
    --vcf: 303k variants complete.
    --vcf: 304k variants complete.
    --vcf: 305k variants complete.
    --vcf: 306k variants complete.
    --vcf: 307k variants complete.
    --vcf: 308k variants complete.
    --vcf: 309k variants complete.
    --vcf: 310k variants complete.
    --vcf: 311k variants complete.
    --vcf: 312k variants complete.
    --vcf: 313k variants complete.
    --vcf: 314k variants complete.
    --vcf: 315k variants complete.
    --vcf: 316k variants complete.
    --vcf: 317k variants complete.
    --vcf: 318k variants complete.
    --vcf: 319k variants complete.
    --vcf: 320k variants complete.
    --vcf: 321k variants complete.
    --vcf: 322k variants complete.
    --vcf: 323k variants complete.
    --vcf: 324k variants complete.
    --vcf: 325k variants complete.
    --vcf: 326k variants complete.
    --vcf: 327k variants complete.
    --vcf: 328k variants complete.
    --vcf: 329k variants complete.
    --vcf: 330k variants complete.
    --vcf: 331k variants complete.
    --vcf: 332k variants complete.
    --vcf: 333k variants complete.
    --vcf: 334k variants complete.
    --vcf: 335k variants complete.
    --vcf: 336k variants complete.
    --vcf: 337k variants complete.
    --vcf: 338k variants complete.
    --vcf: 339k variants complete.
    --vcf: 340k variants complete.
    --vcf: 341k variants complete.
    --vcf: 342k variants complete.
    --vcf: 343k variants complete.
    --vcf: 344k variants complete.
    --vcf: 345k variants complete.
    --vcf: 346k variants complete.
    --vcf: 347k variants complete.
    --vcf: 348k variants complete.
    --vcf: 349k variants complete.
    --vcf: 350k variants complete.
    --vcf: 351k variants complete.
    --vcf: 352k variants complete.
    --vcf: 353k variants complete.
    --vcf: 354k variants complete.
    --vcf: 355k variants complete.
    --vcf: 356k variants complete.
    --vcf: 357k variants complete.
    --vcf: 358k variants complete.
    --vcf: 359k variants complete.
    --vcf: 360k variants complete.
    --vcf: 361k variants complete.
    --vcf: 362k variants complete.
    --vcf: 363k variants complete.
    --vcf: 364k variants complete.
    --vcf: 365k variants complete.
    --vcf: 366k variants complete.
    --vcf: 367k variants complete.
    --vcf: 368k variants complete.
    --vcf: 369k variants complete.
    --vcf: 370k variants complete.
    --vcf: 371k variants complete.
    --vcf: 372k variants complete.
    --vcf: 373k variants complete.
    --vcf: 374k variants complete.
    --vcf: 375k variants complete.
    --vcf: 376k variants complete.
    --vcf: 377k variants complete.
    --vcf: 378k variants complete.
    --vcf: 379k variants complete.
    --vcf: 380k variants complete.
    --vcf: 381k variants complete.
    --vcf: 382k variants complete.
    --vcf: 383k variants complete.
    --vcf: 384k variants complete.
    --vcf: 385k variants complete.
    --vcf: 386k variants complete.
    --vcf: 387k variants complete.
    --vcf: 388k variants complete.
    --vcf: 389k variants complete.
    --vcf: 390k variants complete.
    --vcf: 391k variants complete.
    --vcf: 392k variants complete.
    --vcf: 393k variants complete.
    --vcf: 394k variants complete.
    --vcf: 395k variants complete.
    --vcf: 396k variants complete.
    --vcf: 397k variants complete.
    --vcf: 398k variants complete.
    --vcf: 399k variants complete.
    --vcf: 400k variants complete.
    --vcf: 401k variants complete.
    --vcf: 402k variants complete.
    --vcf: 403k variants complete.
    --vcf: 404k variants complete.
    --vcf: 405k variants complete.
    --vcf: 406k variants complete.
    --vcf: 407k variants complete.
    --vcf: 408k variants complete.
    --vcf: 409k variants complete.
    --vcf: 410k variants complete.
    --vcf: 411k variants complete.
    --vcf: 412k variants complete.
    --vcf: 413k variants complete.
    --vcf: 414k variants complete.
    --vcf: 415k variants complete.
    --vcf: 416k variants complete.
    --vcf: 417k variants complete.
    --vcf: 418k variants complete.
    --vcf: 419k variants complete.
    --vcf: 420k variants complete.
    --vcf: 421k variants complete.
    --vcf: 422k variants complete.
    --vcf: 423k variants complete.
    --vcf: 424k variants complete.
    --vcf: 425k variants complete.
    --vcf: 426k variants complete.
    --vcf: 427k variants complete.
    --vcf: 428k variants complete.
    --vcf: 429k variants complete.
    --vcf: 430k variants complete.
    --vcf: 431k variants complete.
    --vcf: 432k variants complete.
    --vcf: 433k variants complete.
    --vcf: 434k variants complete.
    --vcf: 435k variants complete.
    --vcf: 436k variants complete.
    --vcf: 437k variants complete.
    --vcf: 438k variants complete.
    --vcf: 439k variants complete.
    --vcf: 440k variants complete.
    --vcf: 441k variants complete.
    --vcf: 442k variants complete.
    --vcf: 443k variants complete.
    --vcf: 444k variants complete.
    --vcf: 445k variants complete.
    --vcf: 446k variants complete.
    --vcf: 447k variants complete.
    --vcf: 448k variants complete.
    --vcf: 449k variants complete.
    --vcf: 450k variants complete.
    --vcf: 451k variants complete.
    --vcf: 452k variants complete.
    --vcf: 453k variants complete.
    --vcf: 454k variants complete.
    --vcf: 455k variants complete.
    --vcf: 456k variants complete.
    --vcf: 457k variants complete.
    --vcf: 458k variants complete.
    --vcf: 459k variants complete.
    --vcf: 460k variants complete.
    --vcf: 461k variants complete.
    --vcf: 462k variants complete.
    --vcf: 463k variants complete.
    --vcf: 464k variants complete.
    --vcf: 465k variants complete.
    --vcf: 466k variants complete.
    --vcf: 467k variants complete.
    --vcf: 468k variants complete.
    --vcf: 469k variants complete.
    --vcf: 470k variants complete.
    --vcf: 471k variants complete.
    --vcf: 472k variants complete.
    --vcf: 473k variants complete.
    --vcf: 474k variants complete.
    --vcf: 475k variants complete.
    --vcf: 476k variants complete.
    --vcf: 477k variants complete.
    --vcf: 478k variants complete.
    --vcf: 479k variants complete.
    --vcf: 480k variants complete.
    --vcf: 481k variants complete.
    --vcf: 482k variants complete.
    --vcf: 483k variants complete.
    --vcf: 484k variants complete.
    --vcf: 485k variants complete.
    --vcf: 486k variants complete.
    --vcf: 487k variants complete.
    --vcf: 488k variants complete.
    --vcf: 489k variants complete.
    --vcf: 490k variants complete.
    --vcf: 491k variants complete.
    --vcf: 492k variants complete.
    --vcf: 493k variants complete.
    --vcf: 494k variants complete.
    --vcf:
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1-temporary.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1-temporary.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1-temporary.fam
    written.
    494384 variants loaded from .bim file.
    42 people (0 males, 0 females, 42 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1.nosex
    .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 42 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.163577.
    486495 variants removed due to missing genotype data (--geno).
    7889 variants and 42 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1.fam
    ... 0%0%1%1%2%2%3%3%4%4%5%5%6%6%7%7%8%8%9%9%10%10%11%11%12%12%13%13%14%14%15%15%16%16%17%17%18%18%19%19%20%20%21%21%22%22%23%23%24%24%25%25%26%26%27%27%28%28%29%29%30%30%31%31%32%32%33%33%34%34%35%35%36%36%37%37%38%38%39%39%40%40%41%41%42%42%43%43%44%44%45%45%46%46%47%47%48%48%49%49%50%50%51%51%52%52%53%53%54%54%55%55%56%56%57%57%58%58%59%59%60%60%61%61%62%62%63%63%64%64%65%65%66%66%67%67%68%68%69%69%70%70%71%71%72%72%73%73%74%74%75%75%76%76%77%77%78%78%79%79%80%80%81%81%82%82%83%83%84%84%85%85%86%86%87%87%88%88%89%89%90%90%91%91%92%92%93%93%94%94%95%95%96%96%97%97%98%98%99%done.
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.log.
    Options in effect:
      --allow-extra-chr
      --bfile StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step1
      --make-bed
      --mind 0.10
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    7889 variants loaded from .bim file.
    42 people (0 males, 0 females, 42 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.nosex
    .
    5 people removed due to missing genotype data (--mind).
    IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.irem
    .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 37 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate in remaining samples is 0.979969.
    7889 variants and 37 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.fam
    ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.log.
    Options in effect:
      --allow-extra-chr
      --const-fid 0
      --geno 0.10
      --keep StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.step2.fam
      --make-bed
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink
      --vcf StemCellsInChronicMyeloidLeukemia.SNP.Filtered.SV.vcf.gz

    257405 MB RAM detected; reserving 128702 MB for main workspace.

    --vcf: 1k variants complete.
    --vcf: 2k variants complete.
    --vcf: 3k variants complete.
    --vcf: 4k variants complete.
    --vcf: 5k variants complete.
    --vcf: 6k variants complete.
    --vcf: 7k variants complete.
    --vcf: 8k variants complete.
    --vcf: 9k variants complete.
    --vcf: 10k variants complete.
    --vcf: 11k variants complete.
    --vcf: 12k variants complete.
    --vcf: 13k variants complete.
    --vcf: 14k variants complete.
    --vcf: 15k variants complete.
    --vcf: 16k variants complete.
    --vcf: 17k variants complete.
    --vcf: 18k variants complete.
    --vcf: 19k variants complete.
    --vcf: 20k variants complete.
    --vcf: 21k variants complete.
    --vcf: 22k variants complete.
    --vcf: 23k variants complete.
    --vcf: 24k variants complete.
    --vcf: 25k variants complete.
    --vcf: 26k variants complete.
    --vcf: 27k variants complete.
    --vcf: 28k variants complete.
    --vcf: 29k variants complete.
    --vcf: 30k variants complete.
    --vcf: 31k variants complete.
    --vcf: 32k variants complete.
    --vcf: 33k variants complete.
    --vcf: 34k variants complete.
    --vcf: 35k variants complete.
    --vcf: 36k variants complete.
    --vcf: 37k variants complete.
    --vcf: 38k variants complete.
    --vcf: 39k variants complete.
    --vcf: 40k variants complete.
    --vcf: 41k variants complete.
    --vcf: 42k variants complete.
    --vcf: 43k variants complete.
    --vcf: 44k variants complete.
    --vcf: 45k variants complete.
    --vcf: 46k variants complete.
    --vcf: 47k variants complete.
    --vcf: 48k variants complete.
    --vcf: 49k variants complete.
    --vcf: 50k variants complete.
    --vcf: 51k variants complete.
    --vcf: 52k variants complete.
    --vcf: 53k variants complete.
    --vcf: 54k variants complete.
    --vcf: 55k variants complete.
    --vcf: 56k variants complete.
    --vcf: 57k variants complete.
    --vcf: 58k variants complete.
    --vcf: 59k variants complete.
    --vcf: 60k variants complete.
    --vcf: 61k variants complete.
    --vcf: 62k variants complete.
    --vcf: 63k variants complete.
    --vcf: 64k variants complete.
    --vcf: 65k variants complete.
    --vcf: 66k variants complete.
    --vcf: 67k variants complete.
    --vcf: 68k variants complete.
    --vcf: 69k variants complete.
    --vcf: 70k variants complete.
    --vcf: 71k variants complete.
    --vcf: 72k variants complete.
    --vcf: 73k variants complete.
    --vcf: 74k variants complete.
    --vcf: 75k variants complete.
    --vcf: 76k variants complete.
    --vcf: 77k variants complete.
    --vcf: 78k variants complete.
    --vcf: 79k variants complete.
    --vcf: 80k variants complete.
    --vcf: 81k variants complete.
    --vcf: 82k variants complete.
    --vcf: 83k variants complete.
    --vcf: 84k variants complete.
    --vcf: 85k variants complete.
    --vcf: 86k variants complete.
    --vcf: 87k variants complete.
    --vcf: 88k variants complete.
    --vcf: 89k variants complete.
    --vcf: 90k variants complete.
    --vcf: 91k variants complete.
    --vcf: 92k variants complete.
    --vcf: 93k variants complete.
    --vcf: 94k variants complete.
    --vcf: 95k variants complete.
    --vcf: 96k variants complete.
    --vcf: 97k variants complete.
    --vcf: 98k variants complete.
    --vcf: 99k variants complete.
    --vcf: 100k variants complete.
    --vcf: 101k variants complete.
    --vcf: 102k variants complete.
    --vcf: 103k variants complete.
    --vcf: 104k variants complete.
    --vcf: 105k variants complete.
    --vcf: 106k variants complete.
    --vcf: 107k variants complete.
    --vcf: 108k variants complete.
    --vcf: 109k variants complete.
    --vcf: 110k variants complete.
    --vcf: 111k variants complete.
    --vcf: 112k variants complete.
    --vcf: 113k variants complete.
    --vcf: 114k variants complete.
    --vcf: 115k variants complete.
    --vcf: 116k variants complete.
    --vcf: 117k variants complete.
    --vcf: 118k variants complete.
    --vcf: 119k variants complete.
    --vcf: 120k variants complete.
    --vcf: 121k variants complete.
    --vcf: 122k variants complete.
    --vcf: 123k variants complete.
    --vcf: 124k variants complete.
    --vcf: 125k variants complete.
    --vcf: 126k variants complete.
    --vcf: 127k variants complete.
    --vcf: 128k variants complete.
    --vcf: 129k variants complete.
    --vcf: 130k variants complete.
    --vcf: 131k variants complete.
    --vcf: 132k variants complete.
    --vcf: 133k variants complete.
    --vcf: 134k variants complete.
    --vcf: 135k variants complete.
    --vcf: 136k variants complete.
    --vcf: 137k variants complete.
    --vcf: 138k variants complete.
    --vcf: 139k variants complete.
    --vcf: 140k variants complete.
    --vcf: 141k variants complete.
    --vcf: 142k variants complete.
    --vcf: 143k variants complete.
    --vcf: 144k variants complete.
    --vcf: 145k variants complete.
    --vcf: 146k variants complete.
    --vcf: 147k variants complete.
    --vcf: 148k variants complete.
    --vcf: 149k variants complete.
    --vcf: 150k variants complete.
    --vcf: 151k variants complete.
    --vcf: 152k variants complete.
    --vcf: 153k variants complete.
    --vcf: 154k variants complete.
    --vcf: 155k variants complete.
    --vcf: 156k variants complete.
    --vcf: 157k variants complete.
    --vcf: 158k variants complete.
    --vcf: 159k variants complete.
    --vcf: 160k variants complete.
    --vcf: 161k variants complete.
    --vcf: 162k variants complete.
    --vcf: 163k variants complete.
    --vcf: 164k variants complete.
    --vcf: 165k variants complete.
    --vcf: 166k variants complete.
    --vcf: 167k variants complete.
    --vcf: 168k variants complete.
    --vcf: 169k variants complete.
    --vcf: 170k variants complete.
    --vcf: 171k variants complete.
    --vcf: 172k variants complete.
    --vcf: 173k variants complete.
    --vcf: 174k variants complete.
    --vcf: 175k variants complete.
    --vcf: 176k variants complete.
    --vcf: 177k variants complete.
    --vcf: 178k variants complete.
    --vcf: 179k variants complete.
    --vcf: 180k variants complete.
    --vcf: 181k variants complete.
    --vcf: 182k variants complete.
    --vcf: 183k variants complete.
    --vcf: 184k variants complete.
    --vcf: 185k variants complete.
    --vcf: 186k variants complete.
    --vcf: 187k variants complete.
    --vcf: 188k variants complete.
    --vcf: 189k variants complete.
    --vcf: 190k variants complete.
    --vcf: 191k variants complete.
    --vcf: 192k variants complete.
    --vcf: 193k variants complete.
    --vcf: 194k variants complete.
    --vcf: 195k variants complete.
    --vcf: 196k variants complete.
    --vcf: 197k variants complete.
    --vcf: 198k variants complete.
    --vcf: 199k variants complete.
    --vcf: 200k variants complete.
    --vcf: 201k variants complete.
    --vcf: 202k variants complete.
    --vcf: 203k variants complete.
    --vcf: 204k variants complete.
    --vcf: 205k variants complete.
    --vcf: 206k variants complete.
    --vcf: 207k variants complete.
    --vcf: 208k variants complete.
    --vcf: 209k variants complete.
    --vcf: 210k variants complete.
    --vcf: 211k variants complete.
    --vcf: 212k variants complete.
    --vcf: 213k variants complete.
    --vcf: 214k variants complete.
    --vcf: 215k variants complete.
    --vcf: 216k variants complete.
    --vcf: 217k variants complete.
    --vcf: 218k variants complete.
    --vcf: 219k variants complete.
    --vcf: 220k variants complete.
    --vcf: 221k variants complete.
    --vcf: 222k variants complete.
    --vcf: 223k variants complete.
    --vcf: 224k variants complete.
    --vcf: 225k variants complete.
    --vcf: 226k variants complete.
    --vcf: 227k variants complete.
    --vcf: 228k variants complete.
    --vcf: 229k variants complete.
    --vcf: 230k variants complete.
    --vcf: 231k variants complete.
    --vcf: 232k variants complete.
    --vcf: 233k variants complete.
    --vcf: 234k variants complete.
    --vcf: 235k variants complete.
    --vcf: 236k variants complete.
    --vcf: 237k variants complete.
    --vcf: 238k variants complete.
    --vcf: 239k variants complete.
    --vcf: 240k variants complete.
    --vcf: 241k variants complete.
    --vcf: 242k variants complete.
    --vcf: 243k variants complete.
    --vcf: 244k variants complete.
    --vcf: 245k variants complete.
    --vcf: 246k variants complete.
    --vcf: 247k variants complete.
    --vcf: 248k variants complete.
    --vcf: 249k variants complete.
    --vcf: 250k variants complete.
    --vcf: 251k variants complete.
    --vcf: 252k variants complete.
    --vcf: 253k variants complete.
    --vcf: 254k variants complete.
    --vcf: 255k variants complete.
    --vcf: 256k variants complete.
    --vcf: 257k variants complete.
    --vcf: 258k variants complete.
    --vcf: 259k variants complete.
    --vcf: 260k variants complete.
    --vcf: 261k variants complete.
    --vcf: 262k variants complete.
    --vcf: 263k variants complete.
    --vcf: 264k variants complete.
    --vcf: 265k variants complete.
    --vcf: 266k variants complete.
    --vcf: 267k variants complete.
    --vcf: 268k variants complete.
    --vcf: 269k variants complete.
    --vcf: 270k variants complete.
    --vcf: 271k variants complete.
    --vcf: 272k variants complete.
    --vcf: 273k variants complete.
    --vcf: 274k variants complete.
    --vcf: 275k variants complete.
    --vcf: 276k variants complete.
    --vcf: 277k variants complete.
    --vcf: 278k variants complete.
    --vcf: 279k variants complete.
    --vcf: 280k variants complete.
    --vcf: 281k variants complete.
    --vcf: 282k variants complete.
    --vcf: 283k variants complete.
    --vcf: 284k variants complete.
    --vcf: 285k variants complete.
    --vcf: 286k variants complete.
    --vcf: 287k variants complete.
    --vcf: 288k variants complete.
    --vcf: 289k variants complete.
    --vcf: 290k variants complete.
    --vcf: 291k variants complete.
    --vcf: 292k variants complete.
    --vcf: 293k variants complete.
    --vcf: 294k variants complete.
    --vcf: 295k variants complete.
    --vcf: 296k variants complete.
    --vcf: 297k variants complete.
    --vcf: 298k variants complete.
    --vcf: 299k variants complete.
    --vcf: 300k variants complete.
    --vcf: 301k variants complete.
    --vcf: 302k variants complete.
    --vcf: 303k variants complete.
    --vcf: 304k variants complete.
    --vcf: 305k variants complete.
    --vcf: 306k variants complete.
    --vcf: 307k variants complete.
    --vcf: 308k variants complete.
    --vcf: 309k variants complete.
    --vcf: 310k variants complete.
    --vcf: 311k variants complete.
    --vcf: 312k variants complete.
    --vcf: 313k variants complete.
    --vcf: 314k variants complete.
    --vcf: 315k variants complete.
    --vcf: 316k variants complete.
    --vcf: 317k variants complete.
    --vcf: 318k variants complete.
    --vcf: 319k variants complete.
    --vcf: 320k variants complete.
    --vcf: 321k variants complete.
    --vcf: 322k variants complete.
    --vcf: 323k variants complete.
    --vcf: 324k variants complete.
    --vcf: 325k variants complete.
    --vcf: 326k variants complete.
    --vcf: 327k variants complete.
    --vcf: 328k variants complete.
    --vcf: 329k variants complete.
    --vcf: 330k variants complete.
    --vcf: 331k variants complete.
    --vcf: 332k variants complete.
    --vcf: 333k variants complete.
    --vcf: 334k variants complete.
    --vcf: 335k variants complete.
    --vcf: 336k variants complete.
    --vcf: 337k variants complete.
    --vcf: 338k variants complete.
    --vcf: 339k variants complete.
    --vcf: 340k variants complete.
    --vcf: 341k variants complete.
    --vcf: 342k variants complete.
    --vcf: 343k variants complete.
    --vcf: 344k variants complete.
    --vcf: 345k variants complete.
    --vcf: 346k variants complete.
    --vcf: 347k variants complete.
    --vcf: 348k variants complete.
    --vcf: 349k variants complete.
    --vcf: 350k variants complete.
    --vcf: 351k variants complete.
    --vcf: 352k variants complete.
    --vcf: 353k variants complete.
    --vcf: 354k variants complete.
    --vcf: 355k variants complete.
    --vcf: 356k variants complete.
    --vcf: 357k variants complete.
    --vcf: 358k variants complete.
    --vcf: 359k variants complete.
    --vcf: 360k variants complete.
    --vcf: 361k variants complete.
    --vcf: 362k variants complete.
    --vcf: 363k variants complete.
    --vcf: 364k variants complete.
    --vcf: 365k variants complete.
    --vcf: 366k variants complete.
    --vcf: 367k variants complete.
    --vcf: 368k variants complete.
    --vcf: 369k variants complete.
    --vcf: 370k variants complete.
    --vcf: 371k variants complete.
    --vcf: 372k variants complete.
    --vcf: 373k variants complete.
    --vcf: 374k variants complete.
    --vcf: 375k variants complete.
    --vcf: 376k variants complete.
    --vcf: 377k variants complete.
    --vcf: 378k variants complete.
    --vcf: 379k variants complete.
    --vcf: 380k variants complete.
    --vcf: 381k variants complete.
    --vcf: 382k variants complete.
    --vcf: 383k variants complete.
    --vcf: 384k variants complete.
    --vcf: 385k variants complete.
    --vcf: 386k variants complete.
    --vcf: 387k variants complete.
    --vcf: 388k variants complete.
    --vcf: 389k variants complete.
    --vcf: 390k variants complete.
    --vcf: 391k variants complete.
    --vcf: 392k variants complete.
    --vcf: 393k variants complete.
    --vcf: 394k variants complete.
    --vcf: 395k variants complete.
    --vcf: 396k variants complete.
    --vcf: 397k variants complete.
    --vcf: 398k variants complete.
    --vcf: 399k variants complete.
    --vcf: 400k variants complete.
    --vcf: 401k variants complete.
    --vcf: 402k variants complete.
    --vcf: 403k variants complete.
    --vcf: 404k variants complete.
    --vcf: 405k variants complete.
    --vcf: 406k variants complete.
    --vcf: 407k variants complete.
    --vcf: 408k variants complete.
    --vcf: 409k variants complete.
    --vcf: 410k variants complete.
    --vcf: 411k variants complete.
    --vcf: 412k variants complete.
    --vcf: 413k variants complete.
    --vcf: 414k variants complete.
    --vcf: 415k variants complete.
    --vcf: 416k variants complete.
    --vcf: 417k variants complete.
    --vcf: 418k variants complete.
    --vcf: 419k variants complete.
    --vcf: 420k variants complete.
    --vcf: 421k variants complete.
    --vcf: 422k variants complete.
    --vcf: 423k variants complete.
    --vcf: 424k variants complete.
    --vcf: 425k variants complete.
    --vcf: 426k variants complete.
    --vcf: 427k variants complete.
    --vcf: 428k variants complete.
    --vcf: 429k variants complete.
    --vcf: 430k variants complete.
    --vcf: 431k variants complete.
    --vcf: 432k variants complete.
    --vcf: 433k variants complete.
    --vcf: 434k variants complete.
    --vcf: 435k variants complete.
    --vcf: 436k variants complete.
    --vcf: 437k variants complete.
    --vcf: 438k variants complete.
    --vcf: 439k variants complete.
    --vcf: 440k variants complete.
    --vcf: 441k variants complete.
    --vcf: 442k variants complete.
    --vcf: 443k variants complete.
    --vcf: 444k variants complete.
    --vcf: 445k variants complete.
    --vcf: 446k variants complete.
    --vcf: 447k variants complete.
    --vcf: 448k variants complete.
    --vcf: 449k variants complete.
    --vcf: 450k variants complete.
    --vcf: 451k variants complete.
    --vcf: 452k variants complete.
    --vcf: 453k variants complete.
    --vcf: 454k variants complete.
    --vcf: 455k variants complete.
    --vcf: 456k variants complete.
    --vcf: 457k variants complete.
    --vcf: 458k variants complete.
    --vcf: 459k variants complete.
    --vcf: 460k variants complete.
    --vcf: 461k variants complete.
    --vcf: 462k variants complete.
    --vcf: 463k variants complete.
    --vcf: 464k variants complete.
    --vcf: 465k variants complete.
    --vcf: 466k variants complete.
    --vcf: 467k variants complete.
    --vcf: 468k variants complete.
    --vcf: 469k variants complete.
    --vcf: 470k variants complete.
    --vcf: 471k variants complete.
    --vcf: 472k variants complete.
    --vcf: 473k variants complete.
    --vcf: 474k variants complete.
    --vcf: 475k variants complete.
    --vcf: 476k variants complete.
    --vcf: 477k variants complete.
    --vcf: 478k variants complete.
    --vcf: 479k variants complete.
    --vcf: 480k variants complete.
    --vcf: 481k variants complete.
    --vcf: 482k variants complete.
    --vcf: 483k variants complete.
    --vcf: 484k variants complete.
    --vcf: 485k variants complete.
    --vcf: 486k variants complete.
    --vcf: 487k variants complete.
    --vcf: 488k variants complete.
    --vcf: 489k variants complete.
    --vcf: 490k variants complete.
    --vcf: 491k variants complete.
    --vcf: 492k variants complete.
    --vcf: 493k variants complete.
    --vcf: 494k variants complete.
    --vcf:
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink-temporary.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink-temporary.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink-temporary.fam
    written.
    494384 variants loaded from .bim file.
    42 people (0 males, 0 females, 42 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.nosex
    .
    --keep: 37 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 37 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate in remaining samples is 0.178942.
    483819 variants removed due to missing genotype data (--geno).
    10565 variants and 37 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.bed +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.bim +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.fam
    ... 0%0%1%1%2%2%3%3%4%4%5%5%6%6%7%7%8%8%9%9%10%10%11%11%12%12%13%13%14%14%15%15%16%16%17%17%18%18%19%20%20%21%21%22%22%23%23%24%24%25%25%26%26%27%27%28%28%29%29%30%30%31%31%32%32%33%33%34%34%35%35%36%36%37%37%38%38%39%40%40%41%41%42%42%43%43%44%44%45%45%46%46%47%47%48%48%49%49%50%50%51%51%52%52%53%53%54%54%55%55%56%56%57%57%58%58%59%60%60%61%61%62%62%63%63%64%64%65%65%66%66%67%67%68%68%69%69%70%70%71%71%72%72%73%73%74%74%75%75%76%76%77%77%78%78%79%80%80%81%81%82%82%83%83%84%84%85%85%86%86%87%87%88%88%89%89%90%90%91%91%92%92%93%93%94%94%95%95%96%96%97%97%98%98%99%done.
    1/ Stock SNPs from StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.bim 
       10565 SNPs in StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.bim.

    2/ Stock rs ID from HGDP_1kGP.exon_UTRS.bim
       302006 SNPs in HGDP_1kGP.exon_UTRS.bim.

    3/ Update StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.bim with rs id
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.log.
    Options in effect:
      --allow-extra-chr
      --bfile StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink
      --extract StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.list
      --make-bed
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    10565 variants loaded from .bim file.
    37 people (0 males, 0 females, 37 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.nosex
    .
    --extract: 6589 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 37 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.969737.
    6589 variants and 37 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.bed +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.bim +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.fam
    ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.log.
    Options in effect:
      --bfile HGDP_1kGP.exon_UTRS
      --extract StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.plink.list
      --make-bed
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    302006 variants loaded from .bim file.
    3481 people (0 males, 0 females, 3481 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.nosex
    .
    --extract: 6589 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 3481 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999677.
    6589 variants and 3481 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.bed +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.bim +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.fam
    ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.log.
    Options in effect:
      --allow-no-sex
      --bfile StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1
      --bmerge StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2
      --make-bed
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    37 people loaded from
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.fam.
    3481 people to be merged from
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.fam.
    Of these, 3481 are new, while 0 are present in the base dataset.
    6589 markers loaded from
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp1.bim.
    6589 markers to be merged from
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.tmp2.bim.
    Of these, 0 are new, while 6589 are present in the base dataset.
    Performing single-pass merge (3518 people, 6589 variants).

    Pass 1: fileset #1 complete.
                                                  
    Merged fileset written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP-merge.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP-merge.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP-merge.fam
    .
    6589 variants loaded from .bim file.
    3518 people (0 males, 0 females, 3518 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.nosex
    .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 3518 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999362.
    6589 variants and 3518 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.fam
    ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.log.
    Options in effect:
      --bfile StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP
      --indep-pairwise 50 10 0.1
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    6589 variants loaded from .bim file.
    3518 people (0 males, 0 females, 3518 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.nosex
    .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 3518 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999362.
    6589 variants and 3518 people pass filters and QC.
    Note: No phenotypes present.

    1%
    2%
    4%
    5%
    7%
    8%
    10%
    11%
    13%
    14%
    16%
    17%
    19%
    20%
    21%
    23%
    24%
    26%
    27%
    29%
    30%
    32%
    33%
    35%
    36%
    38%
    39%
    41%
    42%
    43%
    45%
    46%
    48%
    49%
    51%
    52%
    54%
    55%
    57%
    58%
    60%
    61%
    63%
    64%
    65%
    67%
    68%
    70%
    71%
    73%
    74%
    76%
    77%
    79%
    80%
    82%
    83%
    85%
    86%
    87%
    89%
    90%
    92%
    93%
    95%
    96%
    98%
    99%
    Pruned 389 variants from chromosome 1, leaving 293.

    1%
    3%
    5%
    7%
    9%
    11%
    13%
    15%
    17%
    19%
    21%
    22%
    24%
    26%
    28%
    30%
    32%
    34%
    36%
    38%
    40%
    42%
    44%
    45%
    47%
    49%
    51%
    53%
    55%
    57%
    59%
    61%
    63%
    65%
    67%
    68%
    70%
    72%
    74%
    76%
    78%
    80%
    82%
    84%
    86%
    88%
    90%
    91%
    93%
    95%
    97%
    99%
    Pruned 293 variants from chromosome 2, leaving 229.

    2%
    4%
    6%
    9%
    11%
    13%
    16%
    18%
    20%
    23%
    25%
    27%
    30%
    32%
    34%
    37%
    39%
    41%
    44%
    46%
    48%
    51%
    53%
    55%
    58%
    60%
    62%
    65%
    67%
    69%
    72%
    74%
    76%
    79%
    81%
    83%
    86%
    88%
    90%
    93%
    95%
    97%
    Pruned 234 variants from chromosome 3, leaving 195.

    3%
    6%
    9%
    13%
    16%
    19%
    23%
    26%
    29%
    32%
    36%
    39%
    42%
    46%
    49%
    52%
    55%
    59%
    62%
    65%
    69%
    72%
    75%
    78%
    82%
    85%
    88%
    92%
    95%
    98%
    Pruned 164 variants from chromosome 4, leaving 140.

    2%
    5%
    8%
    11%
    14%
    17%
    20%
    23%
    26%
    29%
    32%
    34%
    37%
    40%
    43%
    46%
    49%
    52%
    55%
    58%
    61%
    64%
    67%
    69%
    72%
    75%
    78%
    81%
    84%
    87%
    90%
    93%
    96%
    99%
    Pruned 172 variants from chromosome 5, leaving 171.

    2%
    5%
    8%
    11%
    14%
    16%
    19%
    22%
    25%
    28%
    30%
    33%
    36%
    39%
    42%
    44%
    47%
    50%
    53%
    56%
    58%
    61%
    64%
    67%
    70%
    73%
    75%
    78%
    81%
    84%
    87%
    89%
    92%
    95%
    98%
    Pruned 185 variants from chromosome 6, leaving 171.

    3%
    6%
    9%
    13%
    16%
    19%
    23%
    26%
    29%
    33%
    36%
    39%
    42%
    46%
    49%
    52%
    56%
    59%
    62%
    66%
    69%
    72%
    75%
    79%
    82%
    85%
    89%
    92%
    95%
    99%
    Pruned 162 variants from chromosome 7, leaving 141.

    3%
    7%
    11%
    15%
    19%
    23%
    27%
    31%
    35%
    39%
    43%
    47%
    51%
    55%
    59%
    63%
    67%
    71%
    75%
    79%
    83%
    87%
    91%
    95%
    99%
    Pruned 149 variants from chromosome 8, leaving 102.

    5%
    10%
    15%
    20%
    25%
    30%
    35%
    40%
    45%
    50%
    55%
    60%
    65%
    70%
    75%
    80%
    85%
    90%
    95%
    Pruned 94 variants from chromosome 9, leaving 106.

    3%
    6%
    9%
    13%
    16%
    19%
    23%
    26%
    29%
    33%
    36%
    39%
    43%
    46%
    49%
    53%
    56%
    59%
    63%
    66%
    69%
    73%
    76%
    79%
    83%
    86%
    89%
    93%
    96%
    99%
    Pruned 158 variants from chromosome 10, leaving 143.

    2%
    5%
    8%
    10%
    13%
    16%
    18%
    21%
    24%
    26%
    29%
    32%
    34%
    37%
    40%
    42%
    45%
    48%
    50%
    53%
    56%
    58%
    61%
    64%
    66%
    69%
    72%
    74%
    77%
    80%
    82%
    85%
    88%
    90%
    93%
    96%
    98%
    Pruned 198 variants from chromosome 11, leaving 176.

    2%
    4%
    7%
    9%
    11%
    14%
    16%
    18%
    21%
    23%
    26%
    28%
    30%
    33%
    35%
    37%
    40%
    42%
    45%
    47%
    49%
    52%
    54%
    56%
    59%
    61%
    63%
    66%
    68%
    71%
    73%
    75%
    78%
    80%
    82%
    85%
    87%
    90%
    92%
    94%
    97%
    99%
    Pruned 235 variants from chromosome 12, leaving 187.

    7%
    14%
    21%
    28%
    35%
    42%
    49%
    56%
    63%
    70%
    78%
    85%
    92%
    99%
    Pruned 81 variants from chromosome 13, leaving 60.

    3%
    7%
    11%
    15%
    19%
    23%
    27%
    31%
    35%
    39%
    43%
    47%
    51%
    55%
    59%
    63%
    67%
    71%
    75%
    79%
    83%
    87%
    91%
    95%
    99%
    Pruned 159 variants from chromosome 14, leaving 92.

    5%
    11%
    16%
    22%
    27%
    33%
    38%
    44%
    50%
    55%
    61%
    66%
    72%
    77%
    83%
    88%
    94%
    Pruned 94 variants from chromosome 15, leaving 86.

    3%
    7%
    11%
    15%
    19%
    23%
    27%
    31%
    35%
    39%
    43%
    47%
    51%
    55%
    59%
    63%
    67%
    71%
    75%
    79%
    83%
    87%
    91%
    95%
    99%
    Pruned 148 variants from chromosome 16, leaving 103.

    2%
    5%
    7%
    10%
    13%
    15%
    18%
    20%
    23%
    26%
    28%
    31%
    33%
    36%
    39%
    41%
    44%
    46%
    49%
    52%
    54%
    57%
    59%
    62%
    65%
    67%
    70%
    72%
    75%
    78%
    80%
    83%
    85%
    88%
    91%
    93%
    96%
    98%
    Pruned 228 variants from chromosome 17, leaving 156.

    8%
    17%
    26%
    35%
    44%
    53%
    61%
    70%
    79%
    88%
    97%
    Pruned 48 variants from chromosome 18, leaving 65.

    2%
    5%
    7%
    10%
    12%
    15%
    18%
    20%
    23%
    25%
    28%
    31%
    33%
    36%
    38%
    41%
    44%
    46%
    49%
    51%
    54%
    56%
    59%
    62%
    64%
    67%
    69%
    72%
    75%
    77%
    80%
    82%
    85%
    88%
    90%
    93%
    95%
    98%
    Pruned 202 variants from chromosome 19, leaving 184.

    6%
    12%
    19%
    25%
    32%
    38%
    45%
    51%
    58%
    64%
    71%
    77%
    84%
    90%
    97%
    Pruned 82 variants from chromosome 20, leaving 72.

    10%
    21%
    31%
    42%
    52%
    63%
    73%
    84%
    94%
    Pruned 63 variants from chromosome 21, leaving 32.

    6%
    13%
    20%
    27%
    34%
    40%
    47%
    54%
    61%
    68%
    74%
    81%
    88%
    95%
    Pruned 83 variants from chromosome 22, leaving 64.
    Pruning complete.  3621 of 6589 variants removed.
    Writing...
    Marker lists written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.prune.in
    and
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.prune.out
    .
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.log.
    Options in effect:
      --bfile StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP
      --extract StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.prune.in
      --make-bed
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    6589 variants loaded from .bim file.
    3518 people (0 males, 0 females, 3518 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.nosex
    .
    --extract: 2968 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 3518 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999321.
    2968 variants and 3518 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.bed
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.bim
    +
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.fam
    ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.log.
    Options in effect:
      --bfile StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca
      --out StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca
      --pca
      --pca-cluster-names HGDP_1kGP
      --within StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.cluster

    257405 MB RAM detected; reserving 128702 MB for main workspace.
    2968 variants loaded from .bim file.
    3518 people (0 males, 0 females, 3518 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.nosex
    .
    --within: 2 clusters loaded, covering a total of 3518 people.
    Using up to 63 threads (change this with --threads).
    Before main variant filters, 3518 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999321.
    2968 variants and 3518 people pass filters and QC.
    Note: No phenotypes present.
    --pca-cluster-names/--pca-clusters: 3481 samples specified.

    60 markers complete.
    120 markers complete.
    180 markers complete.
    240 markers complete.
    300 markers complete.
    360 markers complete.
    420 markers complete.
    480 markers complete.
    540 markers complete.
    600 markers complete.
    660 markers complete.
    720 markers complete.
    780 markers complete.
    840 markers complete.
    900 markers complete.
    960 markers complete.
    1020 markers complete.
    1080 markers complete.
    1140 markers complete.
    1200 markers complete.
    1260 markers complete.
    1320 markers complete.
    1380 markers complete.
    1440 markers complete.
    1500 markers complete.
    1560 markers complete.
    1620 markers complete.
    1680 markers complete.
    1740 markers complete.
    1800 markers complete.
    1860 markers complete.
    1920 markers complete.
    1980 markers complete.
    2040 markers complete.
    2100 markers complete.
    2160 markers complete.
    2220 markers complete.
    2280 markers complete.
    2340 markers complete.
    2400 markers complete.
    2460 markers complete.
    2520 markers complete.
    2580 markers complete.
    2640 markers complete.
    2700 markers complete.
    2760 markers complete.
    2820 markers complete.
    2880 markers complete.
    2940 markers complete.
    2968 markers complete.
    Relationship matrix calculation complete.
    [extracting eigenvalues and eigenvectors]
    --pca: Results saved to
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.eigenval
    and
    StemCellsInChronicMyeloidLeukemia/StemCellsInChronicMyeloidLeukemia.HGDP_1kGP.pca.eigenvec
    .

``` r
library(RColorBrewer)

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
populations <- c("Africa", "American", "East Asia", "European", "Middle East", "South Asia")
info$Color <- mycol[match(info$myREG, populations)]
#
x1=c(data_HGDP_1kGP$V3,mydata$V3); y1=c(data_HGDP_1kGP$V4,mydata$V4)
x2=c(data_HGDP_1kGP$V5,mydata$V5); y2=c(data_HGDP_1kGP$V6,mydata$V6)

op <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
on.exit(par(op), add = TRUE)

# PC1 vs PC2
plot(
  data_HGDP_1kGP$V3, data_HGDP_1kGP$V4,
  col = as.character(info$Color), pch = info$PCH,
  xlab = "PC1", ylab = "PC2",
  xlim = c(min(x1), max(x1)), ylim = c(min(y1), max(y1))
)
points(mydata$V3, mydata$V4, pch = 16)
```

![](tutorial_ancestry_files/figure-commonmark/unnamed-chunk-2-1.png)

``` r
# PC3 vs PC4
plot(
  data_HGDP_1kGP$V5, data_HGDP_1kGP$V6,
  col = as.character(info$Color), pch = info$PCH,
  xlab = "PC3", ylab = "PC4",
  xlim = c(min(x2), max(x2)), ylim = c(min(y2), max(y2))
)
points(mydata$V5, mydata$V6, pch = 16)

legend("topright", cex = 0.7, legend = mylegend, col = mycol, pch = 16, bty = "n")
```

![](tutorial_ancestry_files/figure-commonmark/unnamed-chunk-2-2.png)
