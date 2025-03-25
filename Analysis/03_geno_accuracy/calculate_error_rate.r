data = read.table("merged.raw",h=T)[,-(1:6)]

accuracy_table=matrix(0,3,3)
for (i in 1:12){
    accuracy_table = accuracy_table + table(as.integer(data[i,]),as.integer(data[i+12,]))[1:3,1:3]
}

error_rate = 1 - sum(diag(accuracy_table))/sum(accuracy_table)
for (i in 1:3){accuracy_table[i,]=accuracy_table[i,]/sum(accuracy_table[i,])}

out = NULL
for (i in 1:12){
  tmp = c(ncol(data)-sum(is.na(data[i,]+data[12+i,])),sum(data[i,]==data[12+i,],na.rm=T))
  out = rbind(out, c(tmp,(tmp[1]-tmp[2])/tmp[1]))
}
tmp = apply(out,2,sum)
out = rbind(out, c(tmp[1:2],(tmp[1]-tmp[2])/tmp[1]))