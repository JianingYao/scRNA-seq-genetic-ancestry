library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

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

result <- tidyr::crossing(
  maf = c(0.01, 0.02, 0.05),
  r2  = c(0.10, 0.15, 0.20)
)

mapping <- c(eur="Europe", eas="East Asia", amr="America", sas="South Asia", afr="Africa", mid="Middle-East")


parse_maf_r2 <- function(res_dir) {
  folder <- basename(normalizePath(res_dir, winslash="/", mustWork = FALSE))
  maf_str <- sub(".*maf([0-9]+)(?:$|_).*", "\\1", folder)
  r2_str  <- sub(".*r2_0p([0-9]+).*", "\\1", folder)
  maf <- as.numeric(paste0("0.", maf_str))   # "05" -> 0.05
  r2  <- as.numeric(paste0("0.", r2_str))    # "10" -> 0.10
  data.frame(maf = maf, r2 = r2)
}

read_allSNPs_across_continents <- function(base_dir, method_suffix) {
  f <- file.path(base_dir, "allSNPs", paste0("allSNPs_", method_suffix, "_continent_error.txt"))
  if (!file.exists(f)) return(NA_real_)
  x <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  x <- cbind(Continent = rownames(x), x)
  x$pred <- NULL; x$reference <- NULL; rownames(x) <- NULL
  colnames(x)[2] <- "error"
  mean(x$error, na.rm = TRUE)
}

read_sc_study_continent <- function(base_dir, study, method_suffix) {
  f <- file.path(base_dir, study, paste0(study, "_", method_suffix, "_continent_error.txt"))
  if (!file.exists(f)) return(NULL)
  x <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  x <- cbind(Continent = rownames(x), x)
  x$pred <- NULL; x$reference <- NULL; rownames(x) <- NULL
  colnames(x)[2] <- "error"
  x$Continent <- mapping[x$Continent]
  x[, c("Continent", "error")]
}

sc_avg_over_continents <- function(base_dir, method_suffix) {
  df_all <- bind_rows(lapply(study_names, read_sc_study_continent,
                             base_dir = base_dir, method_suffix = method_suffix))
  if (is.null(df_all) || nrow(df_all) == 0) return(NA_real_)
  per_cont <- df_all %>%
    group_by(Continent) %>%
    summarise(mean_err = mean(error, na.rm = TRUE), .groups = "drop")
  mean(per_cont$mean_err, na.rm = TRUE)
}

sc_noise_avg_over_continents <- function(res_dir_nonnoise, method_suffix) {
  results_folder <- basename(normalizePath(res_dir_nonnoise, winslash="/", mustWork = FALSE))
  noise_dir <- file.path("../../Analysis/02_Power_Analysis_noise", results_folder, "Results_0.08")
  df_all <- bind_rows(lapply(study_names, read_sc_study_continent,
                             base_dir = noise_dir, method_suffix = method_suffix))
  if (is.null(df_all) || nrow(df_all) == 0) return(NA_real_)
  per_cont <- df_all %>%
    group_by(Continent) %>%
    summarise(mean_err = mean(error, na.rm = TRUE), .groups = "drop")
  mean(per_cont$mean_err, na.rm = TRUE)
}

results_dirs <- Sys.glob("../../Analysis/02_Power_Analysis/Results*/")

rows <- list()
for (res_dir in results_dirs) {
  key <- parse_maf_r2(res_dir)

  all_pca <- read_allSNPs_across_continents(res_dir, "pca")
  all_rf  <- read_allSNPs_across_continents(res_dir, "pcaRF")
  all_adm <- read_allSNPs_across_continents(res_dir, "admixture")

  sc_pca  <- sc_avg_over_continents(res_dir, "pca")
  sc_rf   <- sc_avg_over_continents(res_dir, "pcaRF")
  sc_adm  <- sc_avg_over_continents(res_dir, "admixture")

  scn_pca <- sc_noise_avg_over_continents(res_dir, "pca")
  scn_rf  <- sc_noise_avg_over_continents(res_dir, "pcaRF")
  scn_adm <- sc_noise_avg_over_continents(res_dir, "admixture")

  rows[[length(rows)+1]] <- data.frame(
    maf = key$maf, r2 = key$r2,
    allSNPs_PCA = all_pca, allSNPs_RandomForest = all_rf, allSNPs_ADMIXTURE = all_adm,
    scSNPs_PCA  = sc_pca,  scSNPs_RandomForest  = sc_rf,  scSNPs_ADMIXTURE  = sc_adm,
    scSNPs8_PCA = scn_pca, scSNPs8_RandomForest = scn_rf, scSNPs8_ADMIXTURE = scn_adm,
    check.names = FALSE
  )
}

combined <- bind_rows(rows)

result <- result %>%
  inner_join(combined, by = c("maf", "r2")) %>%
  arrange(maf, r2)

write.csv(result, "power_MAF_LD.csv", row.names = FALSE)

df <- result

long <- df %>%
  pivot_longer(-c(maf, r2), names_to = "metric", values_to = "error") %>%
  mutate(
    Method  = sub(".*_", "", metric),                  # PCA / RandomForest / ADMIXTURE
    Variant = sub("_[^_]+$", "", metric)               # allSNPs / scSNPs / scSNPs8
  ) %>%
  select(-metric)

long <- long %>%
  mutate(
    Variant = recode(Variant,
      "allSNPs"  = "all-SNPs",
      "scSNPs"   = "sc-SNPs",
      "scSNPs8"  = "sc-SNPs-8%error"
    ),
    Method = recode(Method,
      "PCA"          = "PCA-Distance",
      "RandomForest" = "PCA-RandomForest",
      "ADMIXTURE"    = "ADMIXTURE"
    ),
    combo = sprintf("MAF=%.2f, r²=%.2f", maf, r2)
  )

combo_levels <- long %>%
  distinct(maf, r2) %>%
  arrange(maf, r2) %>%
  transmute(combo = sprintf("MAF=%.2f, r²=%.2f", maf, r2)) %>%
  pull(combo)

long$combo <- factor(long$combo, levels = combo_levels)
long$Method <- factor(long$Method,
                      levels = c("PCA-Distance", "PCA-RandomForest", "ADMIXTURE"))

library(grid) 

p <- ggplot(long, aes(x = combo, y = error, fill = Method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~ Variant, nrow = 1) +
  labs(x = "MAF–LD combination", y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 10, 20),      
    panel.spacing = unit(8, "pt")               
  ) +
  coord_cartesian(clip = "off")                 


png_cairo <- function(filename, width, height, ...) {
  png(filename, width = width, height = height, units = "in",
      res = 300, type = "cairo", ...)
}

ggsave("Supp_Figure_4.png", p, width = 8.5, height = 5, bg = "white",
       device = png_cairo)

