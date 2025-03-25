individual_barcode = read.csv("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/download_data/GSM5899873_OneK1K_scRNA_Sample1_Individual_Barcodes.csv", header = T)
for (individual in unique(individual_barcode[,1])) {
    barcodes <- individual_barcode[individual_barcode$Individual.ID == individual,2]
    barcodes <- gsub("-1$", "", barcodes)
    barcodes <- paste0("CB:Z:", barcodes)
    write.table(barcodes, paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Preprocessing/eQTLAutoimmune/02_AddRG/", individual, "_barcodes.txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}
