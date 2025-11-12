acb = read.table("../Supp_Table_7/supp_table_7_acb_avg.csv",sep=",",h=T)
asw = read.table("../Supp_Table_7/supp_table_7_asw_avg.csv",sep=",",h=T)
mxl = read.table("../Supp_Table_7/supp_table_7_mxl_avg.csv",sep=",",h=T)
bal = read.table("../Supp_Table_7/supp_table_7_balochi_avg.csv",sep=",",h=T)

myplot <- function(data,myanc,mytitle,mycol,myylab=""){
    mydata = subset(data,data$Ancestry==myanc)
    mylim=max(c(mydata$allSNPs,mydata$sc.SNPs,mydata$sc.SNPs.8.error))
    mylim=1
    plot(mydata$allSNPs,mydata$sc.SNPs,pch=16,col="grey",xlab="all-SNPs",ylab=myylab,main=mytitle,bty="l",xlim=c(0,mylim),ylim=c(0,mylim),cex.axis=1.25,cex.lab=1.25)
    abline(0,1,col="grey",lty=2)
    points(mydata$allSNPs,mydata$sc.SNPs.8.error,pch=16,col=mycol)
    legend("topleft",c(
    "sc-SNPs",
    paste0("(r=",round(cor(mydata$allSNPs,mydata$sc.SNPs),2),", s=",round(glm(mydata$sc.SNPs~mydata$allSNPs+0)$coefficient[1],2),")"),
    "sc-SNPs-8%error",
    paste0("(r=",round(cor(mydata$allSNPs,mydata$sc.SNPs.8.error),2),", s=",round(glm(mydata$sc.SNPs.8.error~mydata$allSNPs+0)$coefficient[1],2),")")
    ),col=c("grey",0,mycol,0),pch=16,bty="n")
}

# png("Figure_4.png", width = 12, height = 9, units = "in", res = 300, bg = "white")
pdf("Figure_4.pdf",width=12,height=9)
par(mar=c(4,4.5,5,0.5))
par(mfrow=c(3,4))
myplot(asw,"afr","\n\nAfrican","#1B9E77","sc-SNPs / sc-SNPs-8%error"); mtext("A - ASW\n\n\n\n", side=2, line=2, at=par('usr')[4], las=2, font=2, adj=0)
myplot(asw,"eur","\n\nEuropean","#E6AB02"); 
myplot(acb,"afr","\n\nAfrican","#1B9E77"); mtext("B - ACB\n\n\n\n", side=2, line=2, at=par('usr')[4], las=2, font=2, adj=0)
myplot(acb,"eur","\n\nEuropean","#E6AB02"); 
myplot(mxl,"afr","\n\nAfrican","#1B9E77","sc-SNPs / sc-SNPs-8%error"); mtext("C - MXL\n\n\n\n", side=2, line=2, at=par('usr')[4], las=2, font=2, adj=0)
myplot(mxl,"amr","\n\nAmerican","#D95F02")
myplot(mxl,"eur","\n\nEuropean","#E6AB02")
myplot(mxl,"mid","\n\nMiddle Eastern","#66A61E")
myplot(bal,"afr","\n\nAfrican","#1B9E77","sc-SNPs / sc-SNPs-8%error"); mtext("D - Balochi\n\n\n\n", side=2, line=2, at=par('usr')[4], las=2, font=2, adj=0)
myplot(bal,"eur","\n\nEuropean","#E6AB02")
myplot(bal,"mid","\n\nMiddle Eastern","#66A61E")
myplot(bal,"sas","\n\nSouth Asian","#E7298A")
dev.off()

mylegend = c("Africa", "America", "East Asia", "Europe", "Middle East", "South Asia")
mycol = c("#1B9E77","#D95F02","#7570B3","#E6AB02","#66A61E","#E7298A")

