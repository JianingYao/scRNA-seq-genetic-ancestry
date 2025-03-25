library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(cowplot)
library(RColorBrewer)


# scRNA no error
study_names <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")
scRNA_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis/Results/", study, "/", study, "_error.txt")
    scRNA_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA_result <- rbind(scRNA_result, scRNA_study[1,])
}
data_noError <- as.data.frame(scRNA_result)
colnames(data_noError) <- c("PCA", "RandomForest", "Admixture")
data_noError$Study <- c("Carraro et al.", "Cillo et al.", "Petropoulos et al.", "Giustacchini et al.")
data_noError$Study <- factor(data_noError$Study, 
                              levels = c("Carraro et al.", "Cillo et al.", "Petropoulos et al.", "Giustacchini et al."))
# scRNA with error
scRNA.noise_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis_noise/Results_0.08/", study, "/", study, "_error.txt")
    scRNA.noise_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA.noise_result <- rbind(scRNA.noise_result, scRNA.noise_study[1,])
}
data_error <- as.data.frame(scRNA.noise_result)
colnames(data_error) <- c("PCA", "RandomForest", "Admixture")
data_error$Study <- c("Carraro et al.", "Cillo et al.", "Petropoulos et al.", "Giustacchini et al.")
data_error$Study <- factor(data_error$Study, 
                              levels = c("Carraro et al.", "Cillo et al.", "Petropoulos et al.", "Giustacchini et al."))


n <- 3481

data_noError <- data_noError %>%
  mutate(across(c(PCA, RandomForest, Admixture), ~ . / 100)) %>%
  mutate(
    PCA_se = sqrt((PCA * (1 - PCA)) / n),
    RandomForest_se = sqrt((RandomForest * (1 - RandomForest)) / n),
    Admixture_se = sqrt((Admixture * (1 - Admixture)) / n)
  ) %>%
  mutate(
    PCA_CI = 1.96 * PCA_se,
    RandomForest_CI = 1.96 * RandomForest_se,
    Admixture_CI = 1.96 * Admixture_se
  ) %>%
  mutate(across(c(PCA, RandomForest, Admixture, PCA_CI, RandomForest_CI, Admixture_CI), ~ . * 100)) 

data_error <- data_error %>%
  mutate(across(c(PCA, RandomForest, Admixture), ~ . / 100)) %>%
  mutate(
    PCA_se = sqrt((PCA * (1 - PCA)) / n),
    RandomForest_se = sqrt((RandomForest * (1 - RandomForest)) / n),
    Admixture_se = sqrt((Admixture * (1 - Admixture)) / n)
  ) %>%
  mutate(
    PCA_CI = 1.96 * PCA_se,
    RandomForest_CI = 1.96 * RandomForest_se,
    Admixture_CI = 1.96 * Admixture_se
  ) %>%
  mutate(across(c(PCA, RandomForest, Admixture, PCA_CI, RandomForest_CI, Admixture_CI), ~ . * 100)) 

data_noError_long <- melt(data_noError, id.vars = "Study", measure.vars = c("PCA", "RandomForest", "Admixture"))
data_noError_CI <- melt(data_noError, id.vars = "Study", measure.vars = c("PCA_CI", "RandomForest_CI", "Admixture_CI"))
data_noError_long$CI <- data_noError_CI$value
data_noError_long$variable <- factor(data_noError_long$variable,
                    levels = c("PCA", "RandomForest", "Admixture"),
                    labels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

data_error_long <- melt(data_error, id.vars = "Study", measure.vars = c("PCA", "RandomForest", "Admixture"))
data_error_CI <- melt(data_error, id.vars = "Study", measure.vars = c("PCA_CI", "RandomForest_CI", "Admixture_CI"))
data_error_long$CI <- data_error_CI$value
data_error_long$variable <- factor(data_error_long$variable,
                    levels = c("PCA", "RandomForest", "Admixture"),
                    labels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

max_value <- max(c(data_noError_long$value + data_noError_long$CI, data_error_long$value + data_error_long$CI))

# plot
plot1 <- ggplot(data_noError_long, aes(x = Study, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = value - CI, ymax = value + CI), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "A. Genotypes from sc-SNPs", x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") +
  ylim(0, max_value)

plot2 <- ggplot(data_error_long, aes(x = Study, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = value - CI, ymax = value + CI), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "B. Genotypes from sc-SNPs-8%error", x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") + 
  ylim(0, max_value)

legend <- get_legend(
  ggplot(data_noError_long, aes(x = Study, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    labs(fill = "Method") +
    theme_minimal() +
    theme(legend.position = "right", legend.justification = "center")
)

combined_plots <- plot_grid(plot1, plot2, ncol = 1, align = "v")
final_plot <- plot_grid(combined_plots, legend, ncol = 2, rel_widths = c(1, 0.2))

ggsave("Supp_Figure_2.png", plot = final_plot, width = 8.8, height = 9)



