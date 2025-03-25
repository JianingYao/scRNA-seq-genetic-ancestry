library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)

scRNA_studies <- c("GompertsAirwatCfCells", "HnsccImmuneLandscape", 
                   "humanPreimplantationEmbryos", "StemCellsInChronicMyeloidLeukemia")

info <- read.table("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region/Results/myinfo.split.txt",
                   header = TRUE, sep = "\t")
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
    path <- paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region/Results/allSNPs/pca_results/",
                   "allSNPs.", ANC, ".eigenvec")
  } else if (panel == "scSNPs") {
    path <- paste0("/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/04_power_region/Results/",
                   study, "/pca_results/", study, ".", ANC, ".eigenvec")
  }
  read.table(path, header = FALSE)
}

plot_pca_single <- function(panel, ANC, study = NULL) {
  myinfo <- subset(info, info$myREG == ANC)
  pca    <- get_pca_data(panel, ANC, study)
  mypop  <- unique(myinfo$POP)
  n_groups <- length(mypop)
  
  if (n_groups <= 8) {
    colors <- brewer.pal(n_groups, "Dark2")
  } else {
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_groups)
  }
  
  mycol <- rep("grey", nrow(pca))
  for (i in seq_along(mypop)) {
    idx <- which(myinfo$POP == mypop[i])
    if (length(idx) > 0) {
      mycol[idx] <- colors[i]
    }
  }
  
  disp_title <- display_label(ANC)
  
  plot(pca$V3, pca$V4, col = mycol, pch = 16,
       xlab = "PC1", ylab = "PC2", main = disp_title)
  legend("topleft", legend = mypop, pch = 16, col = colors, cex = 0.8)
}

ancestries <- c("Africa", "American", "EastAsia", "European", "MiddleEast", "SouthAsia")

for (study in scRNA_studies) {
  pdf_filename <- paste0(study, "_Supp_Figure_4.pdf")
  pdf(pdf_filename, width = 12, height = 12)
  
  # We create a 4×3 layout: 
  #  - top 2 rows (6 cells) for Panel A (allSNPs)
  #  - bottom 2 rows (6 cells) for Panel B (scSNPs)
  layout_matrix <- rbind(
    matrix(1:3,   nrow = 1, byrow = TRUE),  # Row 1 of Panel A
    matrix(4:6,   nrow = 1, byrow = TRUE),  # Row 2 of Panel A
    matrix(7:9,   nrow = 1, byrow = TRUE),  # Row 1 of Panel B
    matrix(10:12, nrow = 1, byrow = TRUE)   # Row 2 of Panel B
  )
  layout(layout_matrix)
  
  par(oma = c(0, 0, 3, 0), mar = c(4, 4, 2, 1))
  
  ## -- Panel A: AllSNPs --
  for (ANC in ancestries) {
    plot_pca_single("allSNPs", ANC)
  }
  
  ## -- Panel B: scSNPs --
  for (ANC in ancestries) {
    plot_pca_single("scSNPs", ANC, study = study)
  }
  
  dev.off()
  cat("Saved plot for", study, "as", pdf_filename, "\n")
}
