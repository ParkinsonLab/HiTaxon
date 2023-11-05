#!/usr/bin/env bash
GENUS_NAMES=$1
BWA_PATH=$2
OUTPUT_PATH=$3

#Create directory if not exist
if [ ! -d "$BWA_PATH" ]; then
    echo "BWA directory not found. Creating Directory..."
    mkdir $BWA_PATH
fi

#Prepare sequences for BWA
python "scripts/align/bwa_build.py" $BWA_PATH $OUTPUT_PATH $GENUS_NAMES

#Create BWA indices
"scripts/align/bwa_build.sh" $GENUS_NAMES $BWA_PATH
