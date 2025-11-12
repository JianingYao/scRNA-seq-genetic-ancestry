library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

continents_order <- c("Africa","America","East Asia","Europe","Middle East","South Asia")
study_names <- c(
  "GompertsAirwatCfCells","HnsccImmuneLandscape","HormoneRegulatedNetworksBreast",
  "HumanCellLandscape","humanHeartFailureCellularLandscape_cell",
  "humanHeartFailureCellularLandscape_nuclei","humanPreimplantationEmbryos",
  "nasalMucosaLifespan","NSCL_lesions_tumor_classification","SpatialMapSkin",
  "StemCellsInChronicMyeloidLeukemia"
)
K <- length(study_names)

# ---------- all-SNPs ----------
allSNPs_pca <- read.table("../../Analysis/02_Power_Analysis/Results_maf05_r2_0p10/allSNPs/allSNPs_pca_continent_error.txt",
                          header=TRUE, sep="\t", check.names=FALSE)
n_vec <- allSNPs_pca$reference
allSNPs_pca <- cbind(Continent = rownames(allSNPs_pca), allSNPs_pca)
allSNPs_pca$pred <- NULL; allSNPs_pca$reference <- NULL; rownames(allSNPs_pca) <- NULL
allSNPs_pca$Method <- "PCA-Distance"

allSNPs_rf <- read.table("../../Analysis/02_Power_Analysis/Results_maf05_r2_0p10/allSNPs/allSNPs_pcaRF_continent_error.txt",
                         header=TRUE, sep="\t", check.names=FALSE)
allSNPs_rf <- cbind(Continent = rownames(allSNPs_rf), allSNPs_rf)
allSNPs_rf$pred <- NULL; allSNPs_rf$reference <- NULL; rownames(allSNPs_rf) <- NULL
allSNPs_rf$Method <- "PCA-RandomForest"

allSNPs_adm <- read.table("../../Analysis/02_Power_Analysis/Results_maf05_r2_0p10/allSNPs/allSNPs_admixture_continent_error.txt",
                          header=TRUE, sep="\t", check.names=FALSE)
allSNPs_adm <- cbind(Continent = rownames(allSNPs_adm), allSNPs_adm)
allSNPs_adm$pred <- NULL; allSNPs_adm$reference <- NULL; rownames(allSNPs_adm) <- NULL
allSNPs_adm$Method <- "ADMIXTURE"

allSNPs <- bind_rows(allSNPs_pca, allSNPs_rf, allSNPs_adm)
colnames(allSNPs)[2] <- "error"   
allSNPs$Genotype <- "all-SNPs"

mapping <- c(eur="Europe", eas="East Asia", amr="America", sas="South Asia", afr="Africa", mid="Middle East")
allSNPs$Continent <- mapping[allSNPs$Continent]

continent_counts <- data.frame(
  Continent = c("Europe","East Asia","America","South Asia","Africa","Middle East"),
  n = n_vec,
  check.names = FALSE
)

allSNPs <- allSNPs %>%
  left_join(continent_counts, by="Continent") %>%
  mutate(
    error_prop = error/100,
    se = sqrt(error_prop * (1 - error_prop) / n)
  ) %>%
  select(Continent, Method, Genotype, error_prop, se)

read_cont_err <- function(study, suffix, noise=FALSE) {
  base <- if (noise) "../../Analysis/02_Power_Analysis_noise/Results_maf05_r2_0p10/Results_0.08/"
          else        "../../Analysis/02_Power_Analysis/Results_maf05_r2_0p10/"
  f <- paste0(base, study, "/", study, "_", suffix, "_continent_error.txt")
  x <- read.table(f, header=TRUE, sep="\t", check.names=FALSE)
  x <- cbind(Continent = rownames(x), x)
  x$pred <- NULL; x$reference <- NULL; rownames(x) <- NULL
  colnames(x)[2] <- "error"
  x
}

agg_sc <- function(df, method_name, genotype_label) {
  df %>%
    group_by(Continent) %>%
    summarise(
      mean_err = mean(error, na.rm=TRUE),   
      sd_err   = sd(error,   na.rm=TRUE),   
      .groups="drop"
    ) %>%
    transmute(
      Continent = mapping[Continent],
      Method    = method_name,
      Genotype  = genotype_label,
      error_prop = mean_err/100,
      se = (sd_err/100)/sqrt(K)            
    )
}

# ---------- sc-SNPs (average across studies) ----------
sc_pca <- bind_rows(lapply(study_names, read_cont_err, suffix="pca",   noise=FALSE))
sc_rf  <- bind_rows(lapply(study_names, read_cont_err, suffix="pcaRF", noise=FALSE))
sc_adm <- bind_rows(lapply(study_names, read_cont_err, suffix="admixture", noise=FALSE))

scRNA <- bind_rows(
  agg_sc(sc_pca, "PCA-Distance",   "sc-SNPs"),
  agg_sc(sc_rf,  "PCA-RandomForest","sc-SNPs"),
  agg_sc(sc_adm, "ADMIXTURE",      "sc-SNPs")
)

# ---------- sc-SNPs-8%error (average across studies) ----------
scn_pca <- bind_rows(lapply(study_names, read_cont_err, suffix="pca",   noise=TRUE))
scn_rf  <- bind_rows(lapply(study_names, read_cont_err, suffix="pcaRF", noise=TRUE))
scn_adm <- bind_rows(lapply(study_names, read_cont_err, suffix="admixture", noise=TRUE))

scRNA_noise <- bind_rows(
  agg_sc(scn_pca, "PCA-Distance",   "sc-SNPs-8%error"),
  agg_sc(scn_rf,  "PCA-RandomForest","sc-SNPs-8%error"),
  agg_sc(scn_adm, "ADMIXTURE",      "sc-SNPs-8%error")
)

# ---------- combine ----------
result <- bind_rows(allSNPs, scRNA, scRNA_noise)

result$Genotype  <- factor(result$Genotype,  levels=c("all-SNPs","sc-SNPs","sc-SNPs-8%error"))
result$Method    <- factor(result$Method,    levels=c("PCA-Distance","PCA-RandomForest","ADMIXTURE"))
result$Continent <- factor(result$Continent, levels=continents_order)

result <- result %>%
  mutate(
    y      = 100*error_prop,
    y_lwr  = 100*(error_prop - 1.96*se),
    y_upr  = 100*(error_prop + 1.96*se)
  )

# ----------------------- PLOTTING -----------------------

avg_over_cont_simple <- result %>%
  group_by(Genotype, Method) %>%
  summarise(
    mu      = mean(y/100, na.rm = TRUE),               
    se_comb = sd(y/100, na.rm = TRUE) / sqrt(n()),     
    .groups = "drop"
  ) %>%
  mutate(
    Continent = "Average",
    y     = 100 * mu,
    y_lwr = 100 * (mu - 1.96 * se_comb),
    y_upr = 100 * (mu + 1.96 * se_comb)
  ) %>%
  select(Continent, Method, Genotype, y, y_lwr, y_upr)

plot_df <- result %>%
  select(Continent, Method, Genotype, y, y_lwr, y_upr) %>%
  bind_rows(avg_over_cont_simple)


facet_levels <- c("Africa","America","East Asia","[EMPTY]",
                  "Europe","Middle East","South Asia","Average")

plot_df$Continent <- factor(as.character(plot_df$Continent),
                            levels = facet_levels)

total_n <- sum(continent_counts$n, na.rm = TRUE)
my_labeller <- function(labels) {
  labels$Continent <- sapply(labels$Continent, function(x) {
    n_val <- continent_counts$n[match(x, continent_counts$Continent)]
    if (!is.na(n_val)) {
      paste0(x, " (n = ", n_val, ")")
    } else if (x == "Average") {
      paste0("Average")
    } else {
      ""  # for the empty facet
    }
  })
  labels
}

p1 <- ggplot(plot_df, aes(x = Genotype, y = y, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~ Continent, ncol = 4, labeller = my_labeller, drop = FALSE) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = NULL, y = "Genetic-ancestry inference error rate (%)", fill = "Method") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("Figure_2.png", p1, bg = "white", width = 10, height = 6, units = "in", dpi = 300)


plot_df$y_lwr <- NULL
plot_df$y_upr <- NULL
colnames(plot_df)[4] <- "Error rate (%)"
order_vec <- c("Africa","America","East Asia","Europe","Middle East","South Asia")
ord <- ifelse(plot_df$Continent %in% order_vec,
              match(plot_df$Continent, order_vec),
              ifelse(grepl("^Average", plot_df$Continent), length(order_vec)+1, length(order_vec)+2))
plot_df <- plot_df[order(ord), ]

write.csv(plot_df, "../Supp_Table_4/supp_table_4.csv", row.names = FALSE)