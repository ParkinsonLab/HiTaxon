#!/usr/bin/env bash
GENUS_NAMES=$1
NUM_OF_THREADS=$2
OUTPUT_PATH=$3

#Cluster CDS derived from represenative assemblues with 99% similiarity on a per species basis
for genus in $(cut -f1 $GENUS_NAMES); do
    for species in $(cut -f1 $OUTPUT_PATH/$genus/output_path.txt); do
        if [ -e $species/non_redundant.fa ]; then
            continue
        else
            cd-hit-est -i $species/agg_CDS.fa -o $species/non_redundant.fa -c 0.99 -M 0 -T $NUM_OF_THREADS
        fi
    done
done
