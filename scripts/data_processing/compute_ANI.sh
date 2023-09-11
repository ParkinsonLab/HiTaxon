#!/usr/bin/env bash
GENUS_NAMES=$1
NUM_OF_THREADS=$2
OUTPUT_PATH=$3

#Compute ANI between all same-species assemblies 
for genus in $(cut -f1 $GENUS_NAMES); do
    for species in $(cut -f1 $OUTPUT_PATH/$genus/output_path.txt); do
        if [ -e $species/ani_output.matrix ]; then
            continue
        else
            #If a large number of comparisons (> 400) use alternative fastANI command (--matrix) to avoid downstream memory contstraints
            if [ $(wc -l < $species/query.txt) -gt 400 ]; then
                fastANI --ql $species/query.txt --rl $species/query.txt -o $species/ani_output --matrix -t $NUM_OF_THREADS
            else
                fastANI --ql $species/query.txt --rl $species/query.txt -o $species/ani_output.matrix -t $NUM_OF_THREADS
            fi
        fi
    done
done
