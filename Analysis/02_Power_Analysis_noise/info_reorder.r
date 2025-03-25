INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")
sorted_info <- info[order(info$REG, info$POP, info$IID), ]
write.table(sorted_info,file="sorted_info.txt",col.names=T,row.names=F,quote=F,sep="\t")

unique_combinations <- unique(sorted_info[c("POP", "REG", "myREG")])
sorted_unique_pop <- unique_combinations[order(unique_combinations$REG, unique_combinations$POP), ]
write.table(sorted_unique_pop,file="sorted_listpop.txt",col.names=F,row.names=F, quote=F,sep="\t")

tokeep_IDs <- sorted_info[, c("FID", "IID")]
write.table(tokeep_IDs,file="sorted_idtokeep.txt",col.names=T,row.names=F,quote=F,sep="\t")