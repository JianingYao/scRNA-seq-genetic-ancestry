info="/project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/02_AddRG/humanPreimplantationEmbryos_RG.tsv"
humanPreimplantationEmbryos = read.table(info, sep = "\t", header = TRUE) 

# submit a job for each bam file
for (i in 1:dim(humanPreimplantationEmbryos)[1]) {
    FILE =  humanPreimplantationEmbryos[i, "RGID"]
    RGLB = humanPreimplantationEmbryos[i, "RGLB"]
    RGPL = humanPreimplantationEmbryos[i, "RGPL"]
    RGPU = humanPreimplantationEmbryos[i, "RGPU"]
    RGSM = humanPreimplantationEmbryos[i, "RGSM"]
    RGDS = paste0(humanPreimplantationEmbryos[i, "RGDS1"], '_', humanPreimplantationEmbryos[i, "RGDS2"])

    # Set environment variables
    Sys.setenv(FILE = FILE)
    Sys.setenv(RGLB = RGLB)
    Sys.setenv(RGPL = RGPL)
    Sys.setenv(RGPU = RGPU)
    Sys.setenv(RGSM = RGSM)
    Sys.setenv(RGDS = RGDS)

    system("sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/humanPreimplantationEmbryos/02_AddRG/02_AddRG_humanPreimplantationEmbryos.sh $FILE $RGLB $RGPL $RGPU $RGSM $RGDS")
}

