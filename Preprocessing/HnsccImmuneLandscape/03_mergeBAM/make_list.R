info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/HnsccImmuneLandscape/02_AddRG/HnsccImmuneLandscape_RG.tsv"
HnsccImmuneLandscape = read.table(info, sep = "\t", header = TRUE) 

donors = unique(HnsccImmuneLandscape[, "RGSM"])

for (donor in donors){
    files = HnsccImmuneLandscape[HnsccImmuneLandscape$RGSM == donor, "RGID"]
    bams = paste0("/scratch1/yaojiani/processed/HnsccImmuneLandscape/02_AddRG/", files, ".step2.bam")
    list = paste0("./donors/lists/", donor, ".txt")
    write(bams, file = list)
}