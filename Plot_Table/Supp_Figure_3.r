library(ggplot2)
library(reshape2)
library(RColorBrewer)

study_names <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")
## results with different error levels
error_rates <- seq(0.01, 0.12, by = 0.01)
all_results <- data.frame() 
for (error_rate in error_rates) {
  result <- NULL
  for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis_noise/Results_", 
                        sprintf("%.2f", error_rate), "/",
                        study, "/", study, "_error.txt")
    study_result <- read.table(file_path, header = TRUE, sep = "\t")
    result <- rbind(result, study_result[1,])
  }
  all_results <- rbind(all_results, colMeans(result))
}
colnames(all_results) <- c("pca", "rf", "admixture")
all_results$noise_level <- error_rates
filtered_data <- all_results
data_long <- melt(filtered_data, id.vars = "noise_level", variable.name = "Method", value.name = "Value")
data_long$Method <- factor(data_long$Method,
                    levels = c("pca", "rf", "admixture"),
                    labels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

## result without error
scRNA_result <- NULL
for (study in study_names) {
    file_path <- paste0("../../Analysis/02_Power_Analysis/Results/", study, "/", study, "_error.txt")
    scRNA_study <- read.table(file_path, header = TRUE, sep = "\t")
    scRNA_result <- rbind(scRNA_result, scRNA_study[1,])
}
scRNA_result <- as.data.frame(scRNA_result)
colnames(scRNA_result) <- c("pca", "rf", "admixture")
data_noError <- colMeans(scRNA_result)

## plot with no noise in the line chart
data_noError <- data.frame(
  noise_level = 0.0,
  Method = factor(c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"),
                  levels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE")),
  Value = data_noError
)

data_long <- rbind(data_noError, data_long)


p <- ggplot(data_long, aes(x = noise_level, y = Value, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Genotyping error rate", y = "Genetic-ancestry inference error rate (%)") +
  scale_x_continuous(limits = c(0.0, 0.12), breaks = seq(0.0, 0.12, by = 0.01), minor_breaks = NULL) +
  theme_minimal()

ggsave("Supp_Figure_3.png", plot = p, width = 10, height = 6)



