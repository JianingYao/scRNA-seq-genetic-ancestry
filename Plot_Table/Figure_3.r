library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyr)

study_names <- c(
  "GompertsAirwatCfCells",
  "HnsccImmuneLandscape",
  "HormoneRegulatedNetworksBreast",
  "HumanCellLandscape",
  "humanHeartFailureCellularLandscape_cell",
  "humanHeartFailureCellularLandscape_nuclei",
  "humanPreimplantationEmbryos",
  "nasalMucosaLifespan",
  "NSCL_lesions_tumor_classification",
  "SpatialMapSkin",
  "StemCellsInChronicMyeloidLeukemia"
)
K <- length(study_names)

error_rates <- c(0.00, seq(0.01, 0.13, by = 0.01))

mapping <- c(eur="Europe", eas="East Asia", amr="America", sas="South Asia", afr="Africa", mid="Middle-East")

read_cont_err <- function(study, suffix, er) {
  base <- if (er == 0) {
    "../../Analysis/02_Power_Analysis/Results_maf05_r2_0p10"
  } else {
    file.path("../../Analysis/02_Power_Analysis_noise/Results_maf05_r2_0p10",
              paste0("Results_", sprintf("%.2f", er)))
  }
  f <- file.path(base, study, paste0(study, "_", suffix, "_continent_error.txt"))
  x <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  x <- cbind(Continent = rownames(x), x)
  x$pred <- NULL; x$reference <- NULL; rownames(x) <- NULL
  colnames(x)[2] <- "error"  
  x$Continent <- mapping[x$Continent]
  x[, c("Continent", "error")]
}

methods <- c("pca","pcaRF","admixture")
method_labels <- c(pca="PCA-Distance", pcaRF="PCA-RandomForest", admixture="ADMIXTURE")

rows <- list()
for (er in error_rates) {
  for (m in methods) {
    df_all <- bind_rows(lapply(study_names, read_cont_err, suffix = m, er = er))
    per_cont <- df_all %>%
      group_by(Continent) %>%
      summarise(mean_err = mean(error, na.rm = TRUE), .groups = "drop")
    overall_pct <- mean(per_cont$mean_err, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      noise_level = er,
      Method = method_labels[[m]],
      Value = overall_pct
    )
  }
}

line_df <- bind_rows(rows)
line_df$Method <- factor(line_df$Method,
                         levels = c("PCA-Distance","PCA-RandomForest","ADMIXTURE"))

p <- ggplot(line_df, aes(x = noise_level, y = Value, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Genotyping error rate",
       y = "Genetic-ancestry inference error rate (%)") +
  scale_x_continuous(limits = c(0.00, 0.13),
                     breaks = seq(0.00, 0.13, by = 0.01),
                     minor_breaks = NULL) +
  theme_minimal()

ggsave("Figure_3.png", plot = p, width = 10, height = 6, dpi = 300, bg = "white")

line_df_wide <- line_df %>%
  mutate(noise_level = sprintf("%.2f", noise_level)) %>% 
  arrange(Method, as.numeric(noise_level)) %>%
  pivot_wider(
    id_cols = Method,
    names_from = noise_level,
    values_from = Value
  )

write.csv(line_df_wide, "../Supp_Table_6/supp_table_6.csv", row.names = FALSE)

