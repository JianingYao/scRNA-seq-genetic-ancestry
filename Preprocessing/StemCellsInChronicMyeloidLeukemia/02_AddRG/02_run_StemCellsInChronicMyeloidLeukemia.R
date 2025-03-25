info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/StemCellsInChronicMyeloidLeukemia/02_AddRG/StemCellsInChronicMyeloidLeukemia_RG.tsv"
StemCellsInChronicMyeloidLeukemia = read.table(info, sep = "\t", header = TRUE) 

# submit a job for each bam file
for (i in 1:dim(StemCellsInChronicMyeloidLeukemia)[1]) {
    FILE =  StemCellsInChronicMyeloidLeukemia[i, "RGID"]
    RGLB = StemCellsInChronicMyeloidLeukemia[i, "RGLB"]
    RGPL = StemCellsInChronicMyeloidLeukemia[i, "RGPL"]
    RGPU = StemCellsInChronicMyeloidLeukemia[i, "RGPU"]
    RGSM = StemCellsInChronicMyeloidLeukemia[i, "RGSM"]
    RGDS = paste0(StemCellsInChronicMyeloidLeukemia[i, "RGDS1"], '_', StemCellsInChronicMyeloidLeukemia[i, "RGDS2"])
    RGPM = StemCellsInChronicMyeloidLeukemia[i, "RGPM"] 

    # Set environment variables
    Sys.setenv(FILE = FILE)
    Sys.setenv(RGLB = RGLB)
    Sys.setenv(RGPL = RGPL)
    Sys.setenv(RGPU = RGPU)
    Sys.setenv(RGSM = RGSM)
    Sys.setenv(RGDS = RGDS)
    Sys.setenv(RGPM = RGPM)

    system("sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/StemCellsInChronicMyeloidLeukemia/02_AddRG/02_AddRG_StemCellsInChronicMyeloidLeukemia.sh $FILE $RGLB $RGPL $RGPU $RGSM $RGDS $RGPM")
}

