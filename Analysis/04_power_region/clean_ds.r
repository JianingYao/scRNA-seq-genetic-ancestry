info = read.table("./Results/myinfo.txt",h=T,sep="\t")
#
myt = table(info$POP)
pop.to.keep = names(myt[myt>10])
pop.to.keep = c(pop.to.keep, c("bantukenya", "tuscan", "dai"))
info = info[which(info$POP %in% pop.to.keep),]
#AFR outliers
info = subset(info,info$POP!="acb" & info$POP!="asw")
info = subset(info,info$IID!="NA18912" & info$IID!="NA18913" & info$IID!="NA18914" & info$IID!="NA19240" & info$IID!="v3.1::NA19238" & info$IID!="v3.1::NA19239")
# AMR 
exclude_ids = c("HG02298", "HG02300", "HG02302", "HG02303", "HG01937", "HG01983", "HG01984", "HGDP01009")
info <- info[!info$IID %in% exclude_ids, ]
#EAS outliers
info = subset(info,info$IID!="HG00578" & info$IID!="HG00579" & info$IID!="HG00581" & info$IID!="HG00582" & info$IID!="HG00635" & info$IID!="HG00636" & info$IID!="HG01798" & info$IID!="HG02122" & info$IID!="HG02015" & info$IID!="HG02016" & info$IID!="HG02131" & info$IID!="NA18943" & info$IID!="NA18954" & info$IID!="NA18976" & info$IID!="NA18964" & info$IID!="HGDP00969" & info$IID!="HGDP00959")
#EUR outliers
info = subset(info,info$POP!="ibs" & info$POP!="french")
exclude_ids = c("NA10835", "NA12248", "NA12249", "NA06997", "NA07045", "NA12801", "NA12813", "NA06986", "NA12812", "NA06991", "NA06993", "NA07019", "NA07056", "NA12752", "NA07014", 
                "NA10857", "NA12329", "NA12336", "HGDP01149", "HGDP01152")
info <- info[!info$IID %in% exclude_ids, ]
#MID
info = subset(info,info$IID!="HGDP01273" & info$IID!="HGDP01279")
#SAS
info = subset(info,info$POP!="gih" & info$POP!="itu" & info$POP!="stu" & info$POP!="pjl")
info = subset(info,info$IID!="HG04193" & info$IID!="HGDP00070" & info$IID!="HGDP00199")
write.table(info,file="./Results/myinfo.QC0.txt",quote=F,sep="\t",col.names=T,row.names=F)

