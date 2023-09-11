#!/usr/bin/env bash
GENUS_NAMES=$1
NUM_OF_THREADS=$2
OUTPUT_PATH=$3


#Remove plasmids
: > "${OUTPUT_PATH}/missing_assemblies.txt"

python "scripts/data_processing/cluster_1_prep.py" $OUTPUT_PATH 

while true; do
    count=$(wc -l "${OUTPUT_PATH}/missing_assemblies.txt" | awk '{print $1}')
    echo "You are missing this many assemblies: $count"
    read -p "Do you want to retry downloading assemblies (y/n): " choice
    case $choice in
        [Yy]* )
            "scripts/data_processing/download_redo.sh" $GENUS_NAMES $OUTPUT_PATH
            echo "Exiting the program."
            exit 0
            ;;
        [Nn]* )
            echo "Continuing the program."
            break
            ;;
        * )
            echo "Please enter 'y' or 'n'."
            ;;
    esac
done

#Cluster similiar assemblies using ANI
"scripts/data_processing/compute_ANI.sh" $GENUS_NAMES $NUM_OF_THREADS $OUTPUT_PATH 
python "scripts/data_processing/cluster_2_prep.py" $OUTPUT_PATH

#Cluster coding sequences with 99% similiarity
"scripts/data_processing/cluster_cds.sh" $GENUS_NAMES $NUM_OF_THREADS $OUTPUT_PATH 
