### change pop to asw or acb or balochi

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

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
  "Airway epithelium","Skin","Bone marrow","Lung","Breast",
  "Embryo","Nose","CD45+","Cell landscape","Heart (cell)","Heart (nuclei)"
)
stopifnot(length(study_names) == length(study_cite))
names(study_cite) <- study_names 

RES_ROOT   <- "/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis/Results_maf05_r2_0p10"
NOISE_ROOT <- "/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/02_Power_Analysis_noise/Results_maf05_r2_0p10/Results_0.08"

target_ancestries <- c("eur", "eas", "amr", "sas", "afr", "mid")

sc_cols    <- c(eur = 1, eas = 2, amr = 3, sas = 4, afr = 5, mid = 6)
noise_cols <- c(afr = 1, amr = 2, eas = 3, eur = 4, mid = 5, sas = 6) 

## Threshold on allSNPs to keep rows (set as needed; 0 == keep all)
threshold <- 0

avg_out_csv      <- file.path("../Supp_Table_7/supp_table_7_acb_avg.csv")
by_study_out_csv <- file.path("../Supp_Table_8/supp_table_8_acb_by_study.csv")
summary_txt      <- file.path("../Supp_Table_7/supp_table_7_acb_summary.txt")

acb.allSNPs <- as.data.frame(read.table(
  file.path(RES_ROOT, "allSNPs/admixture_results/allSNPs.noacb.txt"),
  header = FALSE
))
colnames(acb.allSNPs) <- c("eur", "eas", "amr", "sas", "afr", "mid")

INFO_DIR <- "/project2/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/00_TGP/"
info <- read.table(file.path(INFO_DIR, "info.txt"), header = TRUE, sep = "\t")
POP  <- "acb"
acb_iids <- info[info$POP == POP, "IID"]

read_sc_matrix <- function(studies, root_dir, ancestry_col_index){
  mat_list <- lapply(studies, function(STUDY){
    f <- file.path(root_dir, STUDY, "admixture_results", "scSNPs.noacb.txt")
    tab <- read.table(f, header = FALSE)
    tab[, ancestry_col_index, drop = FALSE]
  })
  do.call(cbind, mat_list)
}

combined_list <- list()

for (anc in target_ancestries){
  sc_mat    <- read_sc_matrix(study_names, RES_ROOT,  ancestry_col_index = sc_cols[[anc]])
  sc_mean   <- rowMeans(sc_mat)

  sc_noise_mat  <- read_sc_matrix(study_names, NOISE_ROOT, ancestry_col_index = noise_cols[[anc]])
  sc_noise_mean <- rowMeans(sc_noise_mat)

  all_prop <- acb.allSNPs[[anc]]

  df <- data.frame(
    IID = acb_iids,
    ancestry = anc,
    allSNPs = all_prop,
    `sc-SNPs` = sc_mean,
    `sc-SNPs-8%error` = sc_noise_mean,
    stringsAsFactors = FALSE
  ) %>%
    filter(allSNPs >= threshold)

  combined_list[[anc]] <- df
}

combined_result <- bind_rows(combined_list)
colnames(combined_result) <- c("acb IID", "Ancestry", "allSNPs", "sc-SNPs", "sc-SNPs-8%error")

combined_result <- combined_result %>%
  mutate(Population = "acb") %>%
  relocate(Population, .before = 1) %>%
  relocate(`acb IID`, .after = Population)

write_csv(combined_result, avg_out_csv)

per_study_list <- list()

for (anc in target_ancestries) {
  sc_mat       <- read_sc_matrix(study_names, RES_ROOT,  ancestry_col_index = sc_cols[[anc]])
  sc_noise_mat <- read_sc_matrix(study_names, NOISE_ROOT, ancestry_col_index = noise_cols[[anc]])

  sc_df     <- as.data.frame(sc_mat);       colnames(sc_df)     <- study_names
  sc_err_df <- as.data.frame(sc_noise_mat); colnames(sc_err_df) <- study_names

  sc_long <- sc_df %>%
    mutate(IID = acb_iids) %>%
    pivot_longer(-IID, names_to = "Study", values_to = "sc-SNPs")

  sc_err_long <- sc_err_df %>%
    mutate(IID = acb_iids) %>%
    pivot_longer(-IID, names_to = "Study", values_to = "sc-SNPs-8%error")

  all_prop <- acb.allSNPs[[anc]]
  keep_iids <- acb_iids[all_prop >= threshold]

  per_anc <- sc_long %>%
    inner_join(sc_err_long, by = c("IID", "Study")) %>%
    mutate(
      Population = "acb",
      Ancestry   = anc,
      allSNPs    = all_prop[match(IID, acb_iids)],
      Study_cite = unname(study_cite[Study])
    ) %>%
    filter(IID %in% keep_iids) %>%
    relocate(Population, IID, Ancestry, Study, Study_cite, allSNPs, `sc-SNPs`, `sc-SNPs-8%error`)

  per_study_list[[anc]] <- per_anc
}

per_study_result <- bind_rows(per_study_list)
write_csv(per_study_result, by_study_out_csv)

plot_df <- combined_result %>%
  mutate(Ancestry = factor(Ancestry, levels = target_ancestries)) %>%
  pivot_longer(cols = c(`sc-SNPs`, `sc-SNPs-8%error`),
               names_to = "Genotype",
               values_to = "y_prop")

sink(summary_txt)
cat("Per-ancestry correlations and slopes (averaged sc-SNPs vs all-SNPs)\n")
cat(sprintf("Threshold on allSNPs: %g\n\n", threshold))

for (anc in target_ancestries){
  sub <- combined_result %>% filter(Ancestry == anc)
  x  <- sub$allSNPs
  y1 <- sub$`sc-SNPs`
  y2 <- sub$`sc-SNPs-8%error`

  corr_sc     <- suppressWarnings(cor(x, y1, use = "complete.obs", method = "pearson"))
  corr_sc_err <- suppressWarnings(cor(x, y2, use = "complete.obs", method = "pearson"))

  slope_sc     <- tryCatch(coef(lm(y1 ~ x))[2], error = function(e) NA_real_)
  slope_sc_err <- tryCatch(coef(lm(y2 ~ x))[2], error = function(e) NA_real_)

  msg <- sprintf("[%s] N=%d  r(sc)=%.4f  slope(sc)=%.4f  r(sc_err)=%.4f  slope(sc_err)=%.4f",
                  toupper(anc), nrow(sub), corr_sc, slope_sc, corr_sc_err, slope_sc_err)
  cat(msg, "\n")
  message(msg)
}
sink()

message("Done.")
message("Saved:")
message(paste0("  - Averaged CSV: ", avg_out_csv))
message(paste0("  - Per-study CSV: ", by_study_out_csv, " (includes Study_cite mapping)"))
message(paste0("  - Summary: ", summary_txt))



