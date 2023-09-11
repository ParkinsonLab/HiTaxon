#!/usr/bin/env bash
GENUS_NAMES=$1
OUTPUT_PATH=$2


while IFS=',' read -r genus species genome; do
    rm -rf $OUTPUT_PATH/$genus/$species/$genome
    species=$OUTPUT_PATH/$genus/$species
    echo $species
    echo $genome
    datasets download genome accession $genome --include cds,seq-report --filename $species/$genome.zip --dehydrated
    unzip $species/$genome.zip -d $species/$genome
    datasets rehydrate --directory $species/$genome/
done < "$OUTPUT_PATH/missing_assemblies.txt"

#Create RefSeq metadata files
find $OUTPUT_PATH -type f -name "assembly_data_report.jsonl" | while read fname; do
    dataformat tsv genome --inputfile $fname > "${fname%/*}/genome.tsv"
done

find $OUTPUT_PATH -type f -name "sequence_report.jsonl"|while read fname; do 
    dataformat tsv genome-seq --inputfile $fname > "${fname%/*}/genome_seq.tsv"
done

