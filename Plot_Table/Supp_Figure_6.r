library(ggplot2)
library(dplyr)

study_names <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")
study_cite <- c("Carraro et al.", "Cillo et al.", "Petropoulos et al.", "Giustacchini et al.")

# all-SNPs results
asw.allSNPs <- as.data.frame(read.table("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis/Results/allSNPs/admixture_results/allSNPs.noasw.txt", h=F))
colnames(asw.allSNPs) <- c("eur", "eas", "amr", "sas", "afr", "mid")

study_data <- list()

for (i in seq_along(study_names)) {
    STUDY <- study_names[i]
    STUDY_CITE <- study_cite[i]
    asw.scRNA.study <- read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis/Results/", STUDY, "/admixture_results/scSNPs.noasw.txt"), h=F)
    asw.scRNA.noise.study <- read.table(paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise/Results_0.08/", STUDY, "/admixture_results/scSNPs.noasw.txt"), h=F)

    study_result <- data.frame(
        allSNPs = asw.allSNPs$eur,  
        sc_SNPs = asw.scRNA.study[,1],  
        sc_SNPs_8error = asw.scRNA.noise.study[,4],
        STUDY = STUDY_CITE
    )
    study_data[[i]] <- study_result
}

result <- bind_rows(study_data)

threshold <- 1/8
result <- result %>% filter(allSNPs >= threshold)

# Scatter plot with facet_wrap (2x2 layout)
p <- ggplot(result, aes(x = allSNPs)) +
  geom_point(aes(y = sc_SNPs, color = "sc-SNPs")) +
  geom_point(aes(y = sc_SNPs_8error, color = "sc-SNPs-8%error")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "European-admixture proportions from all-SNPs genotypes", 
    y = "European-admixture proportions from sc-SNPs and sc-SNPs-8%error genotypes", 
    color = "Genotype"
  ) +
  scale_color_manual(values = c("sc-SNPs" = "grey", "sc-SNPs-8%error" = "blue")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.45)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45)) +
  theme_minimal() +
  facet_wrap(~ STUDY, ncol = 2)

ggsave("Supp_Figure_6.png", plot = p, width = 10, height = 10)



