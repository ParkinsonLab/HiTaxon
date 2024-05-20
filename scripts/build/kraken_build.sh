#!/usr/bin/env bash
GENUS_NAMES=$1
KRAKEN_NAME=$2
KRAKEN_PATH=$3
NUM_OF_THREADS=$4
OUTPUT_PATH=$5

kraken2-build --add-to-library $OUTPUT_PATH/"human_kraken.fa" --db $KRAKEN_PATH/$KRAKEN_NAME --threads $NUM_OF_THREADS

#Add sequences to library
for genus in $(cut -f1 $GENUS_NAMES); do
    find $OUTPUT_PATH/$genus -name kraken.fa | while read fname; do
       kraken2-build --add-to-library $fname --db $KRAKEN_PATH/$KRAKEN_NAME --threads $NUM_OF_THREADS;
    done
done

#Build database
kraken2-build --build --db $KRAKEN_PATH/$KRAKEN_NAME --threads $NUM_OF_THREADS --fast-build
