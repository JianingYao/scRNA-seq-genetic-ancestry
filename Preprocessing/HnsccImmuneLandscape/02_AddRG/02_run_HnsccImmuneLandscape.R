info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/02_AddRG/HnsccImmuneLandscape_RG.tsv"
HnsccImmuneLandscape = read.table(info, sep = "\t", header = TRUE) 

# submit a job for each bam file
for (i in 1:dim(HnsccImmuneLandscape)[1]) {
    FILE =  HnsccImmuneLandscape[i, "RGID"]
    RGLB = HnsccImmuneLandscape[i, "RGLB"]
    RGPL = HnsccImmuneLandscape[i, "RGPL"]
    RGPU = HnsccImmuneLandscape[i, "RGPU"]
    RGSM = HnsccImmuneLandscape[i, "RGSM"]
    RGDS = HnsccImmuneLandscape[i, "RGDS"]
    
    # Set environment variables
    Sys.setenv(FILE = FILE)
    Sys.setenv(RGLB = RGLB)
    Sys.setenv(RGPL = RGPL)
    Sys.setenv(RGPU = RGPU)
    Sys.setenv(RGSM = RGSM)
    Sys.setenv(RGDS = RGDS)

    system("sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/02_AddRG/02_AddRG_HnsccImmuneLandscape.sh $FILE $RGLB $RGPL $RGPU $RGSM $RGDS")
}