library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(RColorBrewer)

study_names <- c(
  "GompertsAirwatCfCells",
  "SpatialMapSkin",
  "StemCellsInChronicMyeloidLeukemia",
  "NSCL_lesions_tumor_classification",
  "HormoneRegulatedNetworksBreast",
  "humanPreimplantationEmbryos",
  "nasalMucosaLifespan",
  "HnsccImmuneLandscape",
  "HumanCellLandscape",
  "humanHeartFailureCellularLandscape_cell",
  "humanHeartFailureCellularLandscape_nuclei"
)

study_cite <- c(
  "Airway epithelium",
  "Skin",
  "Bone marrow",
  "Lung",
  "Breast",
  "Embryo",
  "Nose",
  "CD45+",
  "Cell landscape",
  "Heart (cell)",
  "Heart (nuclei)"
)

mapping <- c(eur="Europe", eas="East Asia", amr="America", sas="South Asia", afr="Africa", mid="Middle-East")

methods <- c("pca","pcaRF","admixture")
method_labels <- c(pca="PCA-Distance", pcaRF="PCA-RandomForest", admixture="ADMIXTURE")

base_clean <- "../../Analysis/02_Power_Analysis/Results_maf05_r2_0p10"
base_noise <- "../../Analysis/02_Power_Analysis_noise/Results_maf05_r2_0p10/Results_0.08"

read_study_continent_errors <- function(base_dir, study, method_suffix) {
  f <- file.path(base_dir, study, paste0(study, "_", method_suffix, "_continent_error.txt"))
  if (!file.exists(f)) return(numeric(0))
  x <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  errs <- as.numeric(x[[3]])         
  names(errs) <- rownames(x)         
  errs
}

summarise_study <- function(base_dir, study) {
  rows <- lapply(methods, function(m) {
    v <- read_study_continent_errors(base_dir, study, m)
    v <- v[!is.na(v)]
    mu <- if (length(v)) mean(v) else NA_real_
    data.frame(Method = method_labels[[m]], Mean = mu)
  })
  do.call(rbind, rows)
}

clean_list <- lapply(seq_along(study_names), function(i) {
  st <- study_names[i]; cite <- study_cite[i]
  out <- summarise_study(base_clean, st)
  out$Study <- cite
  out$Variant <- "sc-SNPs"
  out
})
data_clean <- bind_rows(clean_list)

noise_list <- lapply(seq_along(study_names), function(i) {
  st <- study_names[i]; cite <- study_cite[i]
  out <- summarise_study(base_noise, st)
  out$Study <- cite
  out$Variant <- "sc-SNPs-8%error"
  out
})
data_noise <- bind_rows(noise_list)

data_clean$Study <- factor(data_clean$Study, levels = study_cite)
data_noise$Study <- factor(data_noise$Study, levels = study_cite)
data_clean$Method <- factor(data_clean$Method, levels = c("PCA-Distance","PCA-RandomForest","ADMIXTURE"))
data_noise$Method <- factor(data_noise$Method, levels = c("PCA-Distance","PCA-RandomForest","ADMIXTURE"))

min_value <- min(c(data_clean$Mean, data_noise$Mean), na.rm = TRUE)
max_value <- max(c(data_clean$Mean, data_noise$Mean), na.rm = TRUE)

plot_panel <- function(df, title_text, x_axis_text = FALSE) {
  ggplot(df, aes(x = Study, y = Mean, fill = Method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = title_text,
         x = if (x_axis_text) "HCA study" else NULL,
         y = "Genetic-ancestry inference error rate (%)",
         fill = "Method") +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      axis.text.x = if (x_axis_text) element_text(angle = 90, vjust = 0.5, hjust = 1) else element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = if (x_axis_text) 0 else 5, r = 5, b = if (x_axis_text) 5 else 0, l = 5)
    )
}

plot1 <- plot_panel(data_clean, "A. Genotypes from sc-SNPs", x_axis_text = FALSE)
plot2 <- plot_panel(data_noise, "B. Genotypes from sc-SNPs-8%error", x_axis_text = TRUE)

legend <- cowplot::get_legend(
  ggplot(data_clean, aes(x = Study, y = Mean, fill = Method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, na.rm = TRUE) +
    scale_fill_brewer(palette = "Set2") +
    labs(fill = "Method") +
    theme_minimal() +
    theme(legend.position = "right", legend.justification = "center")
)

combined_plots <- plot_grid(plot1, plot2, ncol = 1, align = "v", rel_heights = c(1, 1))
final_plot <- plot_grid(combined_plots, legend, ncol = 2, rel_widths = c(1, 0.2))

ggsave("Supp_Figure_3.png",
       plot = final_plot, width = 8.8, height = 10, dpi = 300, bg = "white")



