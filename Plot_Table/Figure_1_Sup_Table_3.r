library(ggplot2)
library(reshape2)
library(dplyr)


# allSNPs
result <- NULL
allSNPs <- read.table("../../Analysis/02_Power_Analysis/Results/allSNPs/allSNPs_error.txt", header = TRUE, sep = "\t")
result <- rbind(result, allSNPs[1,])
# scRNA
study_names <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")
scRNA_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis/Results/", study, "/", study, "_error.txt")
    scRNA_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA_result <- rbind(scRNA_result, scRNA_study[1,])
}
scRNA <- colMeans(scRNA_result)
result <- rbind(result, scRNA)
# scRNA with noise
scRNA.noise_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis_noise/Results_0.08/", study, "/", study, "_error.txt")
    scRNA.noise_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA.noise_result <- rbind(scRNA.noise_result, scRNA.noise_study[1,])
}
scRNA.noise <- colMeans(scRNA.noise_result)
result <- rbind(result, scRNA.noise)
filtered_data <- as.data.frame(result)
colnames(filtered_data) <- c("PCA", "RandomForest", "ADMIXTURE")
filtered_data$Study <- c("all-SNPs", "sc-SNPs", "sc-SNPs-8%error")
filtered_data$Study <- factor(filtered_data$Study, 
                              levels = c("all-SNPs", "sc-SNPs", "sc-SNPs-8%error"))


n <- 3481*4
filtered_data <- filtered_data %>%
  mutate(across(c(PCA, RandomForest, ADMIXTURE), ~ . / 100))
filtered_data <- filtered_data %>%
  mutate(
    PCA_se = sqrt((PCA * (1 - PCA)) / n),
    RandomForest_se = sqrt((RandomForest * (1 - RandomForest)) / n),
    ADMIXTURE_se = sqrt((ADMIXTURE * (1 - ADMIXTURE)) / n)
  )

data_long <- melt(filtered_data, id.vars = "Study", measure.vars = c("PCA", "RandomForest", "ADMIXTURE"))
data_se <- melt(filtered_data, id.vars = "Study", measure.vars = c("PCA_se", "RandomForest_se", "ADMIXTURE_se"))
data_long$se <- data_se$value
data_long$variable <- factor(data_long$variable,
                    levels = c("PCA", "RandomForest", "ADMIXTURE"),
                    labels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

library(RColorBrewer)
p <- ggplot(data_long, aes(x = Study, y = 100*value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = 100*(value - 1.96*se), ymax = 100*(value + 1.96*se)), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank())
ggsave("Figure_1.png", plot = p, width = 10, height = 6)


supp_table_3 <- data_long %>%
  mutate(
    ErrorRate = 100 * value,
    LowerBound = 100 * (value - 1.96 * se),
    UpperBound = 100 * (value + 1.96 * se)
  ) %>%
  select(Study, Method = variable, ErrorRate, LowerBound, UpperBound)

colnames(supp_table_3) <- c("SNP set", "Method", "Error rate (%)", "SE (lower; %)", "SE (upper; %)")

write.csv(supp_table_3, "../Supp_Table_3/supp_table_3.csv", row.names = FALSE)