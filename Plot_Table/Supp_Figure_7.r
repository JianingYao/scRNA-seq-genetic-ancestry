acb = read.table("../Supp_Table_8/supp_table_8_acb_by_study.csv",sep=",",h=T)
asw = read.table("../Supp_Table_8/supp_table_8_asw_by_study.csv",sep=",",h=T)
mxl = read.table("../Supp_Table_8/supp_table_8_mxl_by_study.csv",sep=",",h=T)
bal = read.table("../Supp_Table_8/supp_table_8_balochi_by_study.csv",sep=",",h=T)

studies = unique(acb$Study_cite)

myplot <- function(study,data,myanc,mytitle,mycol,myylab=""){
    mydata = subset(data,data$Study_cite==study & data$Ancestry==myanc)
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

pdf("Supp_Figure_7.pdf",width=12,height=18)
par(mfrow=c(6,4))
par(mar=c(4,4.5,5,0.5))
for (study in studies){
    myplot(study,asw,"eur","\n\nASW","#E6AB02","sc-SNPs / sc-SNPs-8%error"); mtext(paste0(study,"\n\n\n\n"), side=2, line=2, at=par('usr')[4], las=2, font=2, adj=0)
    myplot(study,acb,"eur","\n\nACB","#E6AB02")
    myplot(study,mxl,"eur","\n\nMXL","#E6AB02")
    myplot(study,bal,"eur","\n\nBalochi","#E6AB02")
}
dev.off() 
