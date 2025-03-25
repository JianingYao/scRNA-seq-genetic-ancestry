data = read.table("../../Analysis/03_geno_accuracy/merged.raw",h=T)[,-(1:6)]

x=matrix(0,3,3)
for (i in 1:12){
    x = x + table(as.integer(data[i,]),as.integer(data[i+12,]))[1:3,1:3]
}

error_rate = sum(diag(x))/sum(x)
# or
error_rate_table = x
for (i in 1:3){error_rate_table[i,]=error_rate_table[i,]/sum(error_rate_table[i,])}

out = NULL
for (i in 1:12){
  tmp = c(ncol(data)-sum(is.na(data[i,]+data[12+i,])),sum(data[i,]==data[12+i,],na.rm=T))
  out = rbind(out, c(tmp,(tmp[1]-tmp[2])/tmp[1]))
}
tmp = apply(out,2,sum)
out = as.data.frame(rbind(out, c(tmp[1:2],(tmp[1]-tmp[2])/tmp[1])))

iid <- c("682_683","683_684","684_685","685_686","686_687","687_688","688_689","689_690","690_691","691_692","692_693","693_694", "All IDs")
supp_table_2 <- cbind(iid, out)
colnames(supp_table_2) <- c("Sample ID", "# of genotypes in both SNP-array and detected in scRNA-seq reads", "# of consistent genotypes", "Genotype error rate")

write.csv(supp_table_2, "supp_table_2.csv", row.names = FALSE)