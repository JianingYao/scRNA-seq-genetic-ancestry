#!/bin/bash

export TMPDIR=/scratch2/yaojiani/job_tmp

skip_first_line=true
while IFS= read -r line
do
    if [ "$skip_first_line" = true ]; then
        skip_first_line=false
        continue
    fi
    RGID=$(echo "$line" | awk -F'\t' '{print $2}')
    PROTOCOL=$(echo "$line" | awk -F'\t' '{print $5}')
    FILE_R1=$(find /scratch1/yaojiani/GompertsAirwatCfCells/data/ -name "*${RGID}_R1.fastq.gz") 
    FILE_R2=$(find /scratch1/yaojiani/GompertsAirwatCfCells/data/ -name "*${RGID}_R2.fastq.gz")

    if [ -z "$FILE_R1" ] || [ -z "$FILE_R2" ]; then
        echo "Missing files for RGID: $RGID"
        continue
    fi

    DIR_R1=$(dirname "$FILE_R1")
    DIR_R2=$(dirname "$FILE_R2")
    if [ "$PROTOCOL" = "Drop-seq" ]; then
        sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/01_STARSolo/01_STARsolo_GompertsAirwatCfCells_dropseq.sh $FILE_R1 $FILE_R2 $RGID $DIR_R1 $DIR_R2
        echo "Process $RGID using $PROTOCOL"
    fi
    if [ "$PROTOCOL" = "10X3v2" ]; then
        sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/01_STARSolo/01_STARsolo_GompertsAirwatCfCells_v2.sh $FILE_R1 $FILE_R2 $RGID $DIR_R1 $DIR_R2
        echo "Process $RGID using $PROTOCOL"
    fi
    if [ "$PROTOCOL" = "10X3v3" ]; then
        sbatch /project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/01_STARSolo/01_STARsolo_GompertsAirwatCfCells_v3.sh $FILE_R1 $FILE_R2 $RGID $DIR_R1 $DIR_R2
        echo "Process $RGID using $PROTOCOL"
    fi
done < /project/gazal_569/jianing/sc-RNA_ancestry/scripts/GompertsAirwatCfCells/02_AddRG/GompertsAirwatCfCells_RG.tsv








