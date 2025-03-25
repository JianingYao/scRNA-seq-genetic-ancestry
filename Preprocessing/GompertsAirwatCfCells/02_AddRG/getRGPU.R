info = "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/02_AddRG/GompertsAirwatCfCells_RG.tsv"
GompertsAirwatCfCells = read.table(info, sep = "\t", header = TRUE)

GompertsAirwatCfCells$RGPU <- rep(0, dim(GompertsAirwatCfCells)[1])

DIR = "/scratch1/yaojiani/GompertsAirwatCfCells"

for (i in 1:dim(GompertsAirwatCfCells)[1]) {
    FILE = GompertsAirwatCfCells[i, "RGID"]
    FOLDER = GompertsAirwatCfCells[i, "folder_name"]
    PROTOCOL = GompertsAirwatCfCells[i, "RGDS1"]

    if (PROTOCOL == "Drop-seq") {
        cmd <- paste0("zcat ", DIR, "/dropseq/r1/data/", FOLDER, "/", FILE, "_R1.fastq.gz | head -n 1 | awk 'NR==1 {split($2,a,\":\"); print a[1] \":\" a[2] \":\" a[3]}'")
    }
    if (PROTOCOL == "10X3v2") {
        cmd <- paste0("zcat ", DIR, "/v2/r1/data/", FOLDER, "/", FILE, "_R1.fastq.gz | head -n 1 | awk 'NR==1 {split($2,a,\":\"); print a[1] \":\" a[2] \":\" a[3] \":\" a[4]}'")
    } 
    if (PROTOCOL == "10X3v3") {
        cmd <- paste0("zcat ", DIR, "/v3/r1/data/", FOLDER, "/", FILE, "_R1.fastq.gz | head -n 1 | awk 'NR==1 {split($2,a,\":\"); print a[1] \":\" a[2] \":\" a[3] \":\" a[4]}'")
    }

    GompertsAirwatCfCells[i, "RGPU"] <- system(cmd, intern = TRUE)
}

write.table(GompertsAirwatCfCells, "/project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/02_AddRG/GompertsAirwatCfCells_RG1.tsv", row.names = FALSE, sep="\t", quote=FALSE)