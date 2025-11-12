library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)

scRNA_studies <- c("GompertsAirwatCfCells", 
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

study_cite <- c("Carraro et al.", 
                 "Cillo et al.", 
                 "Ganier et al.",
                 "Giustacchini et al.",
                 "Han et al.",
                 "Koenig et al. (cell)",
                 "Koenig et al. (nuclei)",
                 "Leader et al.", 
                 "Murrow et al.",
                 "Petropoulos et al.",
                 "Winkley et al.")

names(study_cite) <- scRNA_studies

info <- read.table(
  "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region/Results/myinfo.split.txt",
  header = TRUE, sep = "\t"
)
info$POP <- as.character(info$POP)

display_label <- function(anc) {
  switch(anc,
         "Africa"     = "Africa",
         "American"   = "America",
         "EastAsia"   = "East Asia",
         "European"   = "Europe",
         "MiddleEast" = "Middle East",
         "SouthAsia"  = "South Asia",
         anc)
}

get_pca_data <- function(panel, ANC, study = NULL) {
  if (panel == "allSNPs") {
    path <- paste0(
      "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region/Results/allSNPs/pca_results/",
      "allSNPs.", ANC, ".eigenvec"
    )
  } else if (panel == "scSNPs") {
    path <- paste0(
      "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region/Results/",
      study, "/pca_results/", study, ".", ANC, ".eigenvec"
    )
  } else stop("panel must be 'allSNPs' or 'scSNPs'")
  read.table(path, header = FALSE, stringsAsFactors = FALSE)
}

plot_pca_single <- function(panel, ANC, study = NULL) {
  myinfo_full <- subset(info, info$myREG == ANC)
  pca <- get_pca_data(panel, ANC, study)

  stopifnot(ncol(pca) >= 4)
  colnames(pca)[1:4] <- c("FID","IID","PC1","PC2")

  common <- intersect(myinfo_full$IID, pca$IID)
  myinfo <- myinfo_full[match(common, myinfo_full$IID), , drop = FALSE]
  pca    <- pca[match(common, pca$IID), , drop = FALSE]

  stopifnot(nrow(myinfo) == nrow(pca))

  mypop   <- unique(myinfo$POP)
  ngroups <- length(mypop)
  cols <- if (ngroups <= 8) brewer.pal(ngroups, "Dark2")
          else colorRampPalette(brewer.pal(8, "Dark2"))(ngroups)

  pt_col <- cols[match(myinfo$POP, mypop)]

  main_title <- display_label(ANC)
  plot(pca$PC1, pca$PC2, col = pt_col, pch = 16,
       xlab = "PC1", ylab = "PC2", main = main_title)
  legend("topleft", legend = mypop, pch = 16, col = cols, cex = 0.8, bty = "n")
}

draw_title_row <- function(label) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.02, 0.5, label, adj = 0, cex = 1.6, font = 2)
}

ancestries <- c("Africa", "American", "EastAsia", "European", "MiddleEast", "SouthAsia")

for (study in scRNA_studies) {
  panelB_title <- if (!is.na(study_cite[study])) study_cite[study] else study

  pdf_filename <- paste0(study, "_Supp_Figure_5.pdf")
  pdf(pdf_filename, width = 12, height = 12)

  layout_mat <- rbind(
    c( 1,  1,  1),
    c( 2,  3,  4),
    c( 5,  6,  7),
    c( 8,  8,  8),
    c( 9, 10, 11),
    c(12, 13, 14)
  )
  layout(layout_mat, heights = c(0.5, 1, 1, 0.5, 1, 1), widths = c(1, 1, 1))

  draw_title_row("A. Genotypes from all-SNPs")

  par(mar = c(4, 4, 2, 1))
  for (ANC in ancestries) plot_pca_single("allSNPs", ANC)

  draw_title_row(paste0("B. Genotypes from sc-SNPs: ", panelB_title))

  par(mar = c(4, 4, 2, 1))
  for (ANC in ancestries) plot_pca_single("scSNPs", ANC, study = study)

  dev.off()
  cat("Saved:", pdf_filename, "\n")
}




