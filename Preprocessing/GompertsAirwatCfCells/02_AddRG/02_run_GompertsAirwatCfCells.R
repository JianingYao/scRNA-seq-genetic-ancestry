info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/02_AddRG/GompertsAirwatCfCells_RG.tsv"
GompertsAirwatCfCells = read.table(info, sep = "\t", header = TRUE) 

# submit a job for each bam file
for (i in 1:dim(GompertsAirwatCfCells)[1]) {
    FILE =  GompertsAirwatCfCells[i, "RGID"]
    RGLB = GompertsAirwatCfCells[i, "RGLB"]
    RGPL = GompertsAirwatCfCells[i, "RGPL"]
    RGPU = GompertsAirwatCfCells[i, "RGPU"]
    RGSM = GompertsAirwatCfCells[i, "RGSM"]
    RGDS = paste0(GompertsAirwatCfCells[i, "RGDS1"], '_', GompertsAirwatCfCells[i, "RGDS2"])
    
    # Set environment variables
    Sys.setenv(FILE = FILE)
    Sys.setenv(RGLB = RGLB)
    Sys.setenv(RGPL = RGPL)
    Sys.setenv(RGPU = RGPU)
    Sys.setenv(RGSM = RGSM)
    Sys.setenv(RGDS = RGDS)

    system("sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/02_AddRG/02_AddRG_GompertsAirwatCfCells.sh $FILE $RGLB $RGPL $RGPU $RGSM $RGDS")
}