#!/usr/bin/env bash
GENUS_NAMES=$1 
OUTPUT_PATH=$2 

for genus in $(cut -f1 $GENUS_NAMES); do
    for species in $(cut -f1 "$OUTPUT_PATH/$genus/output_path.txt"); do
        #Use Prodigal to find protein-coding genes from HiTaxon's non-redundant coding sequences and convert to protein
        prodigal -i $species/non_redundant.fa -o $species/gene.coords.gbk -a $species/protein.fa
    done
done
