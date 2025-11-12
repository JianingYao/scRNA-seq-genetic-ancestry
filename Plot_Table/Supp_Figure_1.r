library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

recode_dataset <- function(x) {
  dplyr::recode(
    x,
    "Breast_Mammary_Tissue"      = "Breast Mammary Tissue",
    "Esophagus_Mucosa"           = "Esophagus Mucosa",
    "Esophagus_Muscularis"       = "Esophagus Muscularis",
    "Heart_Left_Ventricle"       = "Heart Left Ventricle",
    "Muscle_Skeletal"            = "Muscle Skeletal",
    "Skin_Sun_Exposed_Lower_leg" = "Skin Sun Exposed Lower leg",
    .default = x
  )
}

strip_x_labels <- function(p) {
  p + theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    plot.margin  = margin(t = 5, r = 5, b = 0, l = 5)
  )
}
keep_x_labels <- function(p) {
  p + theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(),
    plot.margin  = margin(t = 0, r = 5, b = 5, l = 5)
  )
}

#############################################
## PANEL A DATA (Variant-level error rate)
#############################################
summary_df <- read.csv("../Supp_Table_3/genotype_error_summary.csv", header = TRUE)

summary_df$Reference <- NULL
summary_df$n_scSNPs_ref <- NULL
summary_df$n_scSNPs_ref_CDS_UTRs <- NULL

plot_df_A <- summary_df %>%
  filter(Dataset != "Average") %>%
  pivot_longer(-Dataset, names_to = "Variant", values_to = "Error") %>%
  mutate(
    Variant = recode(
      as.character(Variant),
      "error_rate_scSNPs_ref"          = "scSNPs in ref.",
      "error_rate_scSNPs_ref_CDS_UTRs" = "scSNPs in ref. and CDS+UTRs",
      .default = as.character(Variant)
    ),
    Dataset = recode_dataset(as.character(Dataset))
  )

lvl_A <- sort(unique(plot_df_A$Dataset))
lvl_A <- c(setdiff(lvl_A, "PBMCs"), "PBMCs")
plot_df_A <- plot_df_A %>%
  mutate(
    Dataset = factor(Dataset, levels = lvl_A),
    Variant = factor(
      Variant,
      levels = c("scSNPs in ref.", "scSNPs in ref. and CDS+UTRs")
    )
  )

#############################################
## PANEL B DATA (Pair-level error rate)
#############################################
pairs_raw <- read.csv("../Supp_Table_3/genotype_error_pair.csv", header = TRUE)

num_cols_pairs <- setdiff(names(pairs_raw), "Dataset")
avg_row_pairs <- data.frame(
  Dataset = "Average",
  t(colMeans(pairs_raw[, num_cols_pairs, drop = FALSE], na.rm = TRUE)),
  check.names = FALSE
)
pairs_out <- rbind(pairs_raw, avg_row_pairs)
write.csv(pairs_out, "../Supp_Table_3/genotype_error_pair.csv", row.names = FALSE)

plot_df_B <- pairs_raw %>%
  filter(Dataset != "Average") %>%
  pivot_longer(-Dataset, names_to = "Pair", values_to = "Error") %>%
  mutate(
    Dataset = recode_dataset(as.character(Dataset))
  )

lvl_B <- sort(unique(plot_df_B$Dataset))
lvl_B <- c(setdiff(lvl_B, "PBMCs"), "PBMCs")
plot_df_B <- plot_df_B %>%
  mutate(
    Dataset = factor(Dataset, levels = lvl_B)
  )

#############################################
## PANEL C DATA (Allele-level error rate)
#############################################
alleles_raw <- read.csv("../Supp_Table_3/genotype_error_allele.csv", header = TRUE)

num_cols_alleles <- setdiff(names(alleles_raw), "Dataset")
avg_row_alleles <- data.frame(
  Dataset = "Average",
  t(colMeans(alleles_raw[, num_cols_alleles, drop = FALSE], na.rm = TRUE)),
  check.names = FALSE
)
alleles_out <- rbind(alleles_raw, avg_row_alleles)
write.csv(alleles_out, "../Supp_Table_3/genotype_error_allele.csv", row.names = FALSE)

plot_df_C <- alleles_raw %>%
  filter(Dataset != "Average") %>%
  pivot_longer(-Dataset, names_to = "Allele", values_to = "Error") %>%
  mutate(
    Allele = recode(
      as.character(Allele),
      "A_to_C" = "A to C", "A_to_G" = "A to G",
      "C_to_A" = "C to A", "C_to_T" = "C to T",
      "G_to_A" = "G to A", "G_to_T" = "G to T",
      "T_to_C" = "T to C", "T_to_G" = "T to G",
      .default = as.character(Allele)
    ),
    Dataset = recode_dataset(as.character(Dataset))
  )

lvl_C <- sort(unique(plot_df_C$Dataset))
lvl_C <- c(setdiff(lvl_C, "PBMCs"), "PBMCs")
plot_df_C <- plot_df_C %>%
  mutate(
    Dataset = factor(Dataset, levels = lvl_C)
  )

## PANEL A with "Average" group
avg_rows_A <- plot_df_A %>%
  group_by(Variant) %>%
  summarize(Error = mean(Error, na.rm = TRUE), .groups = "drop") %>%
  mutate(Dataset = "Average")

plot_df_A_bar <- bind_rows(
  plot_df_A %>% mutate(Dataset = as.character(Dataset)),
  avg_rows_A
) %>%
  mutate(
    Dataset = recode_dataset(Dataset),
    Dataset = factor(
      Dataset,
      levels = c(levels(plot_df_A$Dataset), "Average")
    ),
    Variant = factor(Variant, levels = levels(plot_df_A$Variant))
  )

p1B <- ggplot(
  plot_df_A_bar,
  aes(x = Dataset, y = Error, fill = Variant)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  labs(
    x = NULL,
    y = "Error rate",
    fill = "Variant",
    subtitle = "Genotype error rate summary (Variants)"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )
p1B <- strip_x_labels(p1B)

## PANEL B with "Average" group
avg_rows_B <- plot_df_B %>%
  group_by(Pair) %>%
  summarize(Error = mean(Error, na.rm = TRUE), .groups = "drop") %>%
  mutate(Dataset = "Average")

plot_df_B_bar <- bind_rows(
  plot_df_B %>% mutate(Dataset = as.character(Dataset)),
  avg_rows_B
) %>%
  mutate(
    Dataset = recode_dataset(Dataset),
    Dataset = factor(
      Dataset,
      levels = c(levels(plot_df_B$Dataset), "Average")
    )
  )

p2B <- ggplot(
  plot_df_B_bar,
  aes(x = Dataset, y = Error, fill = Pair)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  labs(
    x = NULL,
    y = "Error rate",
    fill = "Pair",
    subtitle = "Genotype error rate per polymorphism pair"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )
p2B <- strip_x_labels(p2B)

## PANEL C with "Average" group
avg_rows_C <- plot_df_C %>%
  group_by(Allele) %>%
  summarize(Error = mean(Error, na.rm = TRUE), .groups = "drop") %>%
  mutate(Dataset = "Average")

plot_df_C_bar <- bind_rows(
  plot_df_C %>% mutate(Dataset = as.character(Dataset)),
  avg_rows_C
) %>%
  mutate(
    Dataset = recode_dataset(Dataset),
    Dataset = factor(
      Dataset,
      levels = c(levels(plot_df_C$Dataset), "Average")
    )
  )

p3B <- ggplot(
  plot_df_C_bar,
  aes(x = Dataset, y = Error, fill = Allele)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  labs(
    x = "Tissue",
    y = "Error rate",
    fill = "Allele",
    subtitle = "Genotype error rate per allele"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )
p3B <- keep_x_labels(p3B)

fig_avgbars <- (p1B / p2B / p3B) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "right")

ggsave(
  "Supp_Figure_1.png",
  plot = fig_avgbars,
  width = 10,
  height = 14,
  dpi = 300
)
