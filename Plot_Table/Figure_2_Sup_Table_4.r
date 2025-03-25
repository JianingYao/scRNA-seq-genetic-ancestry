library(ggplot2)

study_names <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")

asw.allSNPs <- as.data.frame(read.table("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis/Results/allSNPs/admixture_results/allSNPs.noasw.txt",h=F))
colnames(asw.allSNPs) <- c("eur", "eas", "amr", "sas", "afr", "mid")

asw.scRNA <- NULL
asw.scRNA.noise <- NULL
for (STUDY in study_names){
    asw.scRNA.study <- read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis/Results/", STUDY, "/admixture_results/scSNPs.noasw.txt"), h=F)
    asw.scRNA <- cbind(asw.scRNA, asw.scRNA.study[,1])
    asw.scRNA.noise.study <- read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise/Results_0.08/", STUDY, "/admixture_results/scSNPs.noasw.txt"), h=F)
    asw.scRNA.noise <- cbind(asw.scRNA.noise, asw.scRNA.noise.study[,4])
}
scRNA <- rowMeans(asw.scRNA)
scRNA.noise <- rowMeans(asw.scRNA.noise)
result <- as.data.frame(cbind(asw.allSNPs$eur, scRNA, scRNA.noise))
colnames(result) <- c("allSNPs", "sc-SNPs", "sc-SNPs-8%error")

INFO_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
info=read.table(paste0(INFO_DIR, "info.txt"),h=T,sep="\t")
POP="asw"
result$`ASW IID` = info[info$POP==POP,]$IID
result <- result[, c("ASW IID", setdiff(names(result), "ASW IID"))]

threshold <- 1/8
selected_indices <- which(result$allSNPs >= threshold)
result <- result[selected_indices, ]

write.csv(result, "supp_table_4.csv", row.names = FALSE)

## scatter plot
p <- ggplot(result, aes(x = allSNPs)) +
  geom_point(aes(y = `sc-SNPs`, color = "sc-SNPs")) +
  geom_point(aes(y = `sc-SNPs-8%error`, color = "sc-SNPs-8%error")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "European-admixture proportions from all-SNPs genotypes", 
       y = "European-admixture proportions from sc-SNPs and sc-SNPs-error8% genotypes", 
       color = "Genotype") +
  scale_color_manual(values = c("sc-SNPs" = "grey", "sc-SNPs-8%error" = "blue")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.45)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45)) +
  theme_minimal()
ggsave("Figure_2.png", plot = p, width = 8, height = 8)

