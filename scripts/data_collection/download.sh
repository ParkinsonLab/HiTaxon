#!/usr/bin/env bash
GENUS_NAMES=$1
OUTPUT_PATH=$2


for genus in $(cut -f1 $GENUS_NAMES); do
    #Create a file listing all species subdirectories
    printf '%s\n' $OUTPUT_PATH/$genus/* > $OUTPUT_PATH/$genus/output_path.txt
    sed -i '/output_path/d' $OUTPUT_PATH/$genus/output_path.txt
    #Download assemblies for each species in genus
    for species in $(cut -f1 $OUTPUT_PATH/$genus/output_path.txt); do
        while read genome; do
            if [ -e "$species/$genome" ]; then
                continue
            else
                datasets download genome accession $genome --include cds,seq-report --filename $species/$genome.zip --dehydrated
                unzip $species/$genome.zip -d $species/$genome
                datasets rehydrate --directory $species/$genome/
            fi
        done < $species/download_assemblies.txt
    done
done 

#Create RefSeq metadata files
find $OUTPUT_PATH -type f -name "assembly_data_report.jsonl" | while read fname; do
    dataformat tsv genome --inputfile $fname > "${fname%/*}/genome.tsv"
done

find $OUTPUT_PATH -type f -name "sequence_report.jsonl"|while read fname; do 
    dataformat tsv genome-seq --inputfile $fname > "${fname%/*}/genome_seq.tsv"
done

