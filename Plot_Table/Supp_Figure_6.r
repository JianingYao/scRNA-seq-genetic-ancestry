library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)

study_names <- c("GompertsAirwatCfCells", 
                 "HnsccImmuneLandscape", 
                 "SpatialMapSkin",
                 "StemCellsInChronicMyeloidLeukemia",
                 "HumanCellLandscape",
                 "humanHeartFailureCellularLandscape_cell",
                 "humanHeartFailureCellularLandscape_nuclei",
                 "NSCL_lesions_tumor_classification", 
                 "HormoneRegulatedNetworksBreast",
                 "humanPreimplantationEmbryos",
                 "nasalMucosaLifespan")
K <- length(study_names)  
                 
continents <- c("Africa", "East Asia", "Europe", "Middle East", "South Asia", "America")
n <- c(701, 665, 536, 124, 278, 159) * 4

scRNA.noise.adm <- read.table(paste0("../../Analysis/04_power_region_noise/Results_0.08/admixture_error.txt"), header = TRUE, row.names = 1)
scRNA.noise.pca <- read.table(paste0("../../Analysis/04_power_region_noise/Results_0.08/pca_error.txt"), header = TRUE, row.names = 1)
scRNA.noise.rf  <- read.table(paste0("../../Analysis/04_power_region_noise/Results_0.08/rf_error.txt"), header = TRUE, row.names = 1)

scRNA.noise.adm$ADMIXTURE         <- rowMeans(scRNA.noise.adm)
scRNA.noise.pca$`PCA-Distance`    <- rowMeans(scRNA.noise.pca)
scRNA.noise.rf$`PCA-RandomForest` <- rowMeans(scRNA.noise.rf)

noise_adm_se <- apply(scRNA.noise.adm[ , setdiff(colnames(scRNA.noise.adm), "ADMIXTURE")],         1, sd, na.rm = TRUE) / sqrt(K)
noise_pca_se <- apply(scRNA.noise.pca[ , setdiff(colnames(scRNA.noise.pca), "PCA-Distance")],      1, sd, na.rm = TRUE) / sqrt(K)
noise_rf_se  <- apply(scRNA.noise.rf  [ , setdiff(colnames(scRNA.noise.rf),  "PCA-RandomForest")], 1, sd, na.rm = TRUE) / sqrt(K)

data <- as.data.frame(cbind(scRNA.noise.pca$`PCA-Distance`, scRNA.noise.rf$`PCA-RandomForest`, scRNA.noise.adm$ADMIXTURE))
colnames(data) <- c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE")
data$Continent  <- continents
data$Continent  <- factor(data$Continent)
data$SampleSize <- n

data$PCA_Distance_se    <- noise_pca_se
data$PCA_RandomForest_se<- noise_rf_se
data$ADMIXTURE_se       <- noise_adm_se

data_long <- melt(data, id.vars = "Continent", measure.vars = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))
data_se   <- melt(data, id.vars = "Continent", measure.vars = c("PCA_Distance_se", "PCA_RandomForest_se", "ADMIXTURE_se"))
data_long$se <- data_se$value
data_long$variable <- factor(data_long$variable,
                             levels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))
max_value <- max(data_long$value + 1.96*data_long$se, na.rm = TRUE)

p_scRNA.noise <- ggplot(data_long, aes(x = Continent, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (value - 1.96*se), ymax = (value + 1.96*se)), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "C. Genotypes from sc-SNPs-8%error", x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") +
  ylim(0, max_value)

scRNA.adm <- read.table(paste0("../../Analysis/04_power_region/Results/admixture_error.txt"), header = TRUE, row.names = 1)
scRNA.pca <- read.table(paste0("../../Analysis/04_power_region/Results/pca_error.txt"), header = TRUE, row.names = 1)
scRNA.rf  <- read.table(paste0("../../Analysis/04_power_region/Results/rf_error.txt"),  header = TRUE, row.names = 1)

scRNA.adm$ADMIXTURE         <- rowMeans(scRNA.adm[, 2:dim(scRNA.adm)[2]])
scRNA.pca$`PCA-Distance`    <- rowMeans(scRNA.pca[, 2:dim(scRNA.pca)[2]])
scRNA.rf$`PCA-RandomForest` <- rowMeans(scRNA.rf[,  2:dim(scRNA.rf)[2]])

study_cols <- 2:(1 + K)

sc_adm_se <- apply(scRNA.adm[, study_cols, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(K)
sc_pca_se <- apply(scRNA.pca[, study_cols, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(K)
sc_rf_se  <- apply(scRNA.rf[,  study_cols, drop = FALSE], 1, sd, na.rm = TRUE) / sqrt(K)

data <- as.data.frame(cbind(scRNA.pca$`PCA-Distance`, scRNA.rf$`PCA-RandomForest`, scRNA.adm$ADMIXTURE))
colnames(data) <- c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE")
data$Continent  <- continents
data$Continent  <- factor(data$Continent)
data$SampleSize <- n

data$PCA_Distance_se     <- sc_pca_se
data$PCA_RandomForest_se <- sc_rf_se
data$ADMIXTURE_se        <- sc_adm_se

data_long <- melt(data, id.vars = "Continent", measure.vars = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))
data_se   <- melt(data, id.vars = "Continent", measure.vars = c("PCA_Distance_se", "PCA_RandomForest_se", "ADMIXTURE_se"))
data_long$se <- data_se$value
data_long$variable <- factor(data_long$variable,
                             levels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

p_scRNA <- ggplot(data_long, aes(x = Continent, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (value - 1.96*se), ymax = (value + 1.96*se)), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "B. Genotypes from sc-SNPs", x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank()) + 
  theme(legend.position = "none") +
  ylim(0, max_value)

allSNPs.adm <- read.table(paste0("../../Analysis/04_power_region/Results/admixture_error.txt"), header = TRUE, row.names = 1)
allSNPs.pca <- read.table(paste0("../../Analysis/04_power_region/Results/pca_error.txt"), header = TRUE, row.names = 1)
allSNPs.rf  <- read.table(paste0("../../Analysis/04_power_region/Results/rf_error.txt"),  header = TRUE, row.names = 1)

allSNPs.adm$ADMIXTURE         <- allSNPs.adm$allSNPs
allSNPs.pca$`PCA-Distance`    <- allSNPs.pca$allSNPs
allSNPs.rf$`PCA-RandomForest` <- allSNPs.rf$allSNPs

data <- as.data.frame(cbind(allSNPs.pca$`PCA-Distance`, allSNPs.rf$`PCA-RandomForest`, allSNPs.adm$ADMIXTURE))
colnames(data) <- c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE")
data$Continent  <- continents
data$Continent  <- factor(data$Continent)
data$SampleSize <- n

data <- data %>%
  mutate(
    PCA_Distance_se     = 100 * sqrt((`PCA-Distance`    / 100) * (1 - (`PCA-Distance`    / 100)) / SampleSize),
    PCA_RandomForest_se = 100 * sqrt((`PCA-RandomForest`/ 100) * (1 - (`PCA-RandomForest`/ 100)) / SampleSize),
    ADMIXTURE_se        = 100 * sqrt((ADMIXTURE         / 100) * (1 - (ADMIXTURE          / 100)) / SampleSize)
  )

data_long <- melt(data, id.vars = "Continent", measure.vars = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))
data_se   <- melt(data, id.vars = "Continent", measure.vars = c("PCA_Distance_se", "PCA_RandomForest_se", "ADMIXTURE_se"))
data_long$se <- data_se$value
data_long$variable <- factor(data_long$variable,
                             levels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

p_allSNPs <- ggplot(data_long, aes(x = Continent, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (value - 1.96*se), ymax = (value + 1.96*se)), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "A. Genotypes from all-SNPs", x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") +
  ylim(0, max_value)

legend <- get_legend(
  ggplot(data_long, aes(x = Continent, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    labs(fill = "Method") +
    theme_minimal() +
    theme(legend.position = "right", legend.justification = "center")
)

combined_plots <- plot_grid(p_allSNPs, p_scRNA, p_scRNA.noise, ncol = 1, align = "v")
final_plot <- plot_grid(combined_plots, legend, ncol = 2, rel_widths = c(1, 0.2))

ggsave("Supp_Figure_6.png", plot = final_plot, width = 9, height = 12, bg = 'white')
