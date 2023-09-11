#!/usr/bin/env bash
OUTPUT_PATH=$1 
KAIJU_PATH=$2 
NUM_OF_THREADS=$3 

#Combine all Kaiju formamtted FASTA files for each species together into a single file
find $OUTPUT_PATH -name kaiju.fa |while read fname; do cat $fname >> "$OUTPUT_PATH/total_kaiju.fa"; done
#Create Kaiju database
kaiju-mkbwt -a ACDEFGHIKLMNPQRSTVWY -o "$KAIJU_PATH/proteins" "$OUTPUT_PATH/total_kaiju.fa" 
kaiju-mkfmi "$KAIJU_PATH/proteins"
