library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

continents <- c("Africa", "America", "East Asia", "Europe","Middle-East", "South Asia")
study_names <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")

# all-SNPs
allSNPs_pca <- read.table("../../Analysis/02_Power_Analysis/Results/allSNPs/allSNPs_pca_continent_error.txt", header = TRUE, sep = "\t")
n <- allSNPs_pca$reference
allSNPs_pca <- cbind("Continent" = rownames(allSNPs_pca), allSNPs_pca)
allSNPs_pca$pred <- NULL
allSNPs_pca$reference <- NULL
rownames(allSNPs_pca) <- NULL
allSNPs_pca$Method <- "PCA-Distance"
#
allSNPs_rf <- read.table("../../Analysis/02_Power_Analysis/Results/allSNPs/allSNPs_pcaRF_continent_error.txt", header = TRUE, sep = "\t")
allSNPs_rf <- cbind("Continent" = rownames(allSNPs_rf), allSNPs_rf)
allSNPs_rf$pred <- NULL
allSNPs_rf$reference <- NULL
rownames(allSNPs_rf) <- NULL
allSNPs_rf$Method <- "PCA-RandomForest"
#
allSNPs_adm <- read.table("../../Analysis/02_Power_Analysis/Results/allSNPs/allSNPs_admixture_continent_error.txt", header = TRUE, sep = "\t")
allSNPs_adm <- cbind("Continent" = rownames(allSNPs_adm), allSNPs_adm)
allSNPs_adm$pred <- NULL
allSNPs_adm$reference <- NULL
rownames(allSNPs_adm) <- NULL
allSNPs_adm$Method <- "ADMIXTURE"

allSNPs <- rbind(allSNPs_pca, allSNPs_rf, allSNPs_adm)
colnames(allSNPs)[2] <- "error"
allSNPs$Genotype <- "all-SNPs"

# scRNA
scRNA_pca_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis/Results/", study, "/", study, "_pca_continent_error.txt")
    scRNA_pca_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA_pca_study <- cbind("Continent" = rownames(scRNA_pca_study), scRNA_pca_study)
    scRNA_pca_study$pred <- NULL
    scRNA_pca_study$reference <- NULL
    rownames(scRNA_pca_study) <- NULL
    scRNA_pca_result <- rbind(scRNA_pca_result, scRNA_pca_study)
}
colnames(scRNA_pca_result)[2] <- "error"
scRNA_pca_result <- as.data.frame(scRNA_pca_result %>%
  group_by(Continent) %>%
  summarise(error = mean(error, na.rm = TRUE)))
scRNA_pca_result$Method <- "PCA-Distance"
#
scRNA_rf_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis/Results/", study, "/", study, "_pcaRF_continent_error.txt")
    scRNA_rf_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA_rf_study <- cbind("Continent" = rownames(scRNA_rf_study), scRNA_rf_study)
    scRNA_rf_study$pred <- NULL
    scRNA_rf_study$reference <- NULL
    rownames(scRNA_rf_study) <- NULL
    scRNA_rf_result <- rbind(scRNA_rf_result, scRNA_rf_study)
}
colnames(scRNA_rf_result)[2] <- "error"
scRNA_rf_result <- as.data.frame(scRNA_rf_result %>%
  group_by(Continent) %>%
  summarise(error = mean(error, na.rm = TRUE)))
scRNA_rf_result$Method <- "PCA-RandomForest"
#
scRNA_adm_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis/Results/", study, "/", study, "_admixture_continent_error.txt")
    scRNA_adm_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA_adm_study <- cbind("Continent" = rownames(scRNA_adm_study), scRNA_adm_study)
    scRNA_adm_study$pred <- NULL
    scRNA_adm_study$reference <- NULL
    rownames(scRNA_adm_study) <- NULL
    scRNA_adm_result <- rbind(scRNA_adm_result, scRNA_adm_study)
}
colnames(scRNA_adm_result)[2] <- "error"
scRNA_adm_result <- as.data.frame(scRNA_adm_result %>%
  group_by(Continent) %>%
  summarise(error = mean(error, na.rm = TRUE)))
scRNA_adm_result$Method <- "ADMIXTURE"

scRNA <- rbind(scRNA_pca_result, scRNA_rf_result, scRNA_adm_result)
scRNA$Genotype <- "sc-SNPs"


# scRNA with noise
scRNA.noise_pca_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis_noise/Results_0.08/", study, "/", study, "_pca_continent_error.txt")
    scRNA.noise_pca_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA.noise_pca_study <- cbind("Continent" = rownames(scRNA.noise_pca_study), scRNA.noise_pca_study)
    scRNA.noise_pca_study$pred <- NULL
    scRNA.noise_pca_study$reference <- NULL
    rownames(scRNA.noise_pca_study) <- NULL
    scRNA.noise_pca_result <- rbind(scRNA.noise_pca_result, scRNA.noise_pca_study)
}
colnames(scRNA.noise_pca_result)[2] <- "error"
scRNA.noise_pca_result <- as.data.frame(scRNA.noise_pca_result %>%
  group_by(Continent) %>%
  summarise(error = mean(error, na.rm = TRUE)))
scRNA.noise_pca_result$Method <- "PCA-Distance"
#
scRNA.noise_rf_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis_noise/Results_0.08/", study, "/", study, "_pcaRF_continent_error.txt")
    scRNA.noise_rf_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA.noise_rf_study <- cbind("Continent" = rownames(scRNA.noise_rf_study), scRNA.noise_rf_study)
    scRNA.noise_rf_study$pred <- NULL
    scRNA.noise_rf_study$reference <- NULL
    rownames(scRNA.noise_rf_study) <- NULL
    scRNA.noise_rf_result <- rbind(scRNA.noise_rf_result, scRNA.noise_rf_study)
}
colnames(scRNA.noise_rf_result)[2] <- "error"
scRNA.noise_rf_result <- as.data.frame(scRNA.noise_rf_result %>%
  group_by(Continent) %>%
  summarise(error = mean(error, na.rm = TRUE)))
scRNA.noise_rf_result$Method <- "PCA-RandomForest"
#
scRNA.noise_adm_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis_noise/Results_0.08/", study, "/", study, "_admixture_continent_error.txt")
    scRNA.noise_adm_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA.noise_adm_study <- cbind("Continent" = rownames(scRNA.noise_adm_study), scRNA.noise_adm_study)
    scRNA.noise_adm_study$pred <- NULL
    scRNA.noise_adm_study$reference <- NULL
    rownames(scRNA.noise_adm_study) <- NULL
    scRNA.noise_adm_result <- rbind(scRNA.noise_adm_result, scRNA.noise_adm_study)
}
colnames(scRNA.noise_adm_result)[2] <- "error"
scRNA.noise_adm_result <- as.data.frame(scRNA.noise_adm_result %>%
  group_by(Continent) %>%
  summarise(error = mean(error, na.rm = TRUE)))
scRNA.noise_adm_result$Method <- "ADMIXTURE"

scRNA.noise <- rbind(scRNA.noise_pca_result, scRNA.noise_rf_result, scRNA.noise_adm_result)
scRNA.noise$Genotype <- "sc-SNPs-8%error"

result <- rbind(allSNPs, scRNA, scRNA.noise)
mapping = c("eur" = "Europe", 
            "eas" = "East Asia", 
            "amr" = "America", 
            "sas" = "South Asia", 
            "afr" = "Africa", 
            "mid" = "Middle-East")
result$Continent = mapping[result$Continent]
result$error <- result$error / 100
continent_counts <- data.frame(
  Continent = c("Europe", "East Asia", "America", "South Asia", "Africa", "Middle-East"),
  n
)
result <- result %>%
  left_join(continent_counts, by = "Continent") %>%
  mutate(
    se = sqrt((error * (1 - error)) / (n*4))
  )
result$Genotype <- factor(result$Genotype, 
                              levels = c("all-SNPs", "sc-SNPs", "sc-SNPs-8%error"))
result$Method <- factor(result$Method,
                    levels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))


my_labeller <- function(labels) {
  labels$Continent <- sapply(labels$Continent, function(x) {
    n_val <- continent_counts$n[match(x, continent_counts$Continent)]
    paste0(x, " (n = ", n_val, ")")
  })
  return(labels)
}



p <- ggplot(result, aes(x = Genotype, y = 100*error, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = 100*(error - 1.96*se), ymax = 100*(error + 1.96*se)),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  facet_wrap(~ Continent, labeller = my_labeller) +              
  scale_fill_brewer(palette = "Dark2") +  
  labs(x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_bw() +                        
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave("Supp_Figure_1.png", plot = p)

