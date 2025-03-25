info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/02_AddRG/humanPreimplantationEmbryos_RG.tsv"
humanPreimplantationEmbryos = read.table(info, sep = "\t", header = TRUE) 

donors = unique(humanPreimplantationEmbryos[, "RGSM"])

for (donor in donors){
    files = humanPreimplantationEmbryos[humanPreimplantationEmbryos$RGSM == donor, "RGID"]
    bams = paste0("/scratch2/yaojiani/processed/humanPreimplantationEmbryos/02_AddRG/", files, ".step2.bam")
    list = paste0("./donors/lists/", donor, ".txt")
    write(bams, file = list)
}