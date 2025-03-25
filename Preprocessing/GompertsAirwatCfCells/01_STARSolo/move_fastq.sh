#!/bin/bash

# Function to move directories
move_directories() {
    local folder_list=$1
    local dest_dir=$2

    # Check if the destination directory exists
    if [ ! -d "$dest_dir" ]; then
        echo "Destination directory does not exist. Creating it now..."
        mkdir -p "$dest_dir"
    fi

    # Loop through each line in the file
    while IFS= read -r folder; do
        # Move the folder to the destination directory
        if [ -d "/scratch1/yaojiani/GompertsAirwatCfCells/$folder" ]; then
            mv "/scratch1/yaojiani/GompertsAirwatCfCells/$folder" "$dest_dir/"
        else
            echo "Directory $folder does not exist."
        fi
    done < "$folder_list"
}


move_directories "dropseq_r1.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/dropseq/r1/data"
move_directories "dropseq_r2.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/dropseq/r2/data"
move_directories "v2_i1.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/v2/i1/data"
move_directories "v2_r1.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/v2/r1/data"
move_directories "v2_r2.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/v2/r2/data"
move_directories "v3_i1.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/v3/i1/data"
move_directories "v3_r1.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/v3/r1/data"
move_directories "v3_r2.txt" "/scratch1/yaojiani/GompertsAirwatCfCells/v3/r2/data"