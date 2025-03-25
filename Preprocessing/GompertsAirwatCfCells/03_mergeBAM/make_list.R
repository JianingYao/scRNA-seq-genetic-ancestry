info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/02_AddRG/GompertsAirwatCfCells_RG.tsv"
GompertsAirwatCfCells = read.table(info, sep = "\t", header = TRUE) 

donors = unique(GompertsAirwatCfCells[, "RGSM"])

for (donor in donors){
    files = GompertsAirwatCfCells[GompertsAirwatCfCells$RGSM == donor, "RGID"]
    bams = paste0("/scratch1/yaojiani/processed/GompertsAirwatCfCells/02_AddRG/", files, ".step2.bam")
    list = paste0("./donors/lists/", donor, ".txt")
    write(bams, file = list)
}