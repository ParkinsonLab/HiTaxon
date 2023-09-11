#!/usr/bin/env bash
ASSEMBLY_SUMMARY=$1
GENUS_NAMES=$2
OUTPUT_PATH=$3

#Create directory to store RefSeq Sequences
if [ ! -d "$OUTPUT_PATH" ]; then
    mkdir $OUTPUT_PATH
fi

#Download RefSeq Assembly Summary File if file path is not provided

if [ "$ASSEMBLY_SUMMARY" = "default" ] && [ ! -e "$OUTPUT_PATH/assembly_summary.txt" ]; then
    ASSEMBLY_SUMMARY="$OUTPUT_PATH/assembly_summary.txt"
    wget "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" -O $ASSEMBLY_SUMMARY
fi

#Filter for high-quality assemblies, pertaining to genera of interest
python "scripts/data_collection/collection.py" $ASSEMBLY_SUMMARY $GENUS_NAMES $OUTPUT_PATH

#Download assemblies
"scripts/data_collection/download.sh" $GENUS_NAMES $OUTPUT_PATH 

