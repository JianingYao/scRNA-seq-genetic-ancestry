info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/StemCellsInChronicMyeloidLeukemia/02_AddRG/StemCellsInChronicMyeloidLeukemia_RG.tsv"
StemCellsInChronicMyeloidLeukemia = read.table(info, sep = "\t", header = TRUE) 

donors = unique(StemCellsInChronicMyeloidLeukemia[, "RGSM"])

for (donor in donors){
    files = StemCellsInChronicMyeloidLeukemia[StemCellsInChronicMyeloidLeukemia$RGSM == donor, "RGID"]
    bams = paste0("/scratch2/yaojiani/processed/StemCellsInChronicMyeloidLeukemia/02_AddRG/", files, ".step2.bam")
    list = paste0("./donors/lists/", donor, ".txt")
    write(bams, file = list)
}