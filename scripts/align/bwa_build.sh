#!/usr/bin/env bash
GENUS_NAMES=$1
BWA_PATH=$2

#Add sequences to library
for genus in $(cut -f1 $BWA_PATH/genus2add.txt); do
    bwa index $BWA_PATH/${genus}.fa
done

