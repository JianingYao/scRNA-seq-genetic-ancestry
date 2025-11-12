study_names <- c("GompertsAirwatCfCells", 
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

DATA_DIR="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/Analysis/01_Ancestry/Results_maf05_r2_0p10/"
result = NULL

for (i in 1:length(study_names)){
  STUDY = study_names[i]
  ADM_DIR=paste0(DATA_DIR, STUDY, "/ADMIXTURE_results/")
  admixture = read.table(paste0(ADM_DIR, STUDY,".TGP_HGDP.pca.6.Q"),h=F)
  admixture.ind = read.table(paste0(DATA_DIR, STUDY,"/",STUDY,".TGP_HGDP.pca.fam"))[,2]

  colnames(admixture)=c("ADM Europe","ADM East Asia","ADM America","ADM South Asia","ADM Africa","ADM Middle-East")
  admixture = admixture[, order(colnames(admixture))]
  rownames(admixture)=admixture.ind

  myadmixture  = admixture[as.character(read.table(paste0(DATA_DIR, STUDY,"/",STUDY,".plink.fam"))[,2]),] 
  myadmixture = cbind("Donor ID" = rownames(myadmixture), myadmixture)
  myadmixture$Study = study_cite[i]
  study_result = myadmixture[, c("Study", setdiff(names(myadmixture), "Study"))]
  study_result$`ADMIXTURE prediction` <- apply(
  study_result[, c("ADM Europe", "ADM East Asia", "ADM America", "ADM South Asia", "ADM Africa", "ADM Middle-East")],1,
  function(x) {
    adm_pred <- names(x)[which.max(x)]
    gsub("ADM ", "", adm_pred)
  }
)

  PCA_DIR = paste0(DATA_DIR, STUDY, "/PCA_results/")
  pca_result = read.table(paste0(PCA_DIR, STUDY,".pca_center.txt"),h=T)
  pca_mapping = c("eur" = "Europe", 
                    "eas" = "East Asia", 
                    "amr" = "America", 
                    "sas" = "South Asia", 
                    "afr" = "Africa", 
                    "mea" = "Middle-East")
  colnames(pca_result) = pca_mapping[colnames(pca_result)]
  pca_prediction = apply(pca_result, 1, function(row_vals) {
    cols_with_1 = colnames(pca_result)[which(row_vals == 1)]
    paste(cols_with_1, collapse = ", ")
  })
  study_result$`PCA prediction` = pca_prediction

  rf_DIR = paste0(DATA_DIR, STUDY, "/rf_results/")
  rf_result = read.table(paste0(rf_DIR, STUDY,".pca_rf.txt"),header=TRUE, row.names=1, sep="\t")
  colnames(rf_result) = c("RF Africa", "Rf America", "RF East Asia", "RF Europe","RF Middle-East", "RF South Asia", "Random Forest prediction")
  rf_mapping = c("European" = "Europe", 
                  "East Asia" = "East Asia", 
                  "American" = "America", 
                  "South Asia" = "South Asia", 
                  "Africa" = "Africa", 
                  "Middle-East" = "Middle-East")
  rf_result$`Random Forest prediction` = rf_mapping[rf_result$`Random Forest prediction`]
  study_result = merge(study_result, rf_result, by = "row.names", all = TRUE)
  study_result$Row.names = NULL
  pred_col = c("PCA prediction", "ADMIXTURE prediction", "Random Forest prediction")
  cols = names(study_result)
  cols = cols[!(cols %in% pred_col)]
  new_order = append(cols, pred_col, after = 2)
  study_result = study_result[, new_order]

  result = rbind(result, study_result)
}

write.csv(result, "supp_table_5.csv", row.names = FALSE)

admixture_mean = mean(result$`ADM Europe`)
admixture_median = median(result$`ADM Europe`) 
admixture_50 = sum(result$`ADM Europe` > 0.5)/nrow(result) 
admixture_80 = sum(result$`ADM Europe` > 0.8)/nrow(result) 
admixture_europe = sum(result$`ADM Europe` > 0.5) 
admixture_africa = sum(result$`ADM Africa` > 0.5) 
admixture_america = sum(result$`ADM America` > 0.5)
admixture_eastasia = sum(result$`ADM East Asia` > 0.5) 
admixture_mideast = sum(result$`ADM Middle-East` > 0.5) 
admixture_southasia = sum(result$`ADM South Asia` > 0.5) 
pca_europe = sum(result$`PCA prediction` == "Europe")/nrow(result) 
rf_50 = sum(result$`RF Europe` > 0.5)/nrow(result) 
rf_80 = sum(result$`RF Europe` > 0.8)/nrow(result) 
rf_90 = sum(result$`RF Europe` > 0.9)/nrow(result) 