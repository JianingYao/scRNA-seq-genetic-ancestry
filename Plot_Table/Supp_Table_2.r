library(dplyr)

INFO_DIR <- "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
info <- as.data.frame(read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t"))

supp_table1 <- info %>% dplyr::select(IID, myREG, POP)
colnames(supp_table1) <- c("IID", "Continent", "Population")

write.csv(supp_table1, "supp_table_2.csv", row.names = FALSE)