mynames = c("Breast Mammary Tissue", "Esophagus Mucosa", "Esophagus Muscularis", "Heart Left Ventricle", "Lung", "Muscle Skeletal", "Prostate", "Skin Sun Exposed Lower Leg", "PBMCs")


error = c(0.095, 0.065, 0.071, 0.060, 0.083, 0.078, 0.067, 0.126, 0.078) 
nb = c(187313, 345258, 365079, 434793, 409038, 324407, 184319, 87656, 2503)
se = sqrt(error*(1-error)/nb)

png("Figure_1.png", width = 8, height = 8, units = "in", res = 300, bg = "white")
par(mar = c(5, 14, 4, 2))
# par(mar=c(5,12,4,2))
barplot(rev(error),horiz = T,names=rev(mynames),xlab="Genotype error rate",border=0,las=2,col="#1B9E77",xlim=c(0,0.13))
abline(v=mean(error),col=2,lty=2)
for(i in 1:9){
    pos = (9-i)*1.2+0.7
    lines(error[i]+1.96*c(-1,1)*se[i],rep(pos,2))
}
dev.off()


