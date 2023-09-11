#!/usr/bin/env bash
GENUS_NAMES=$1
OUTPUT_PATH=$2

SIM_NUM=500000
#Simulate paired-end reads for species
for genus in $(cut -f1 "$GENUS_NAMES"); do
    for species in $(cut -f1 "$OUTPUT_PATH/$genus/output_path.txt"); do
        #Simulate a minimum of SIM_NUM reads
        TOTAL=$(seqkit stats "$species/non_redundant.fa" -T | cut -f5 | tail -n 1)
        if [[ $TOTAL -lt $SIM_NUM*150 ]]; then
            SEQCOUNT=$(echo "600000 * 150 * 2 / $TOTAL" | bc)
            SEQCOUNT=$(($SEQCOUNT + 1))
        else
            SEQCOUNT=1
        fi
        species_name=$(basename "$species")
        if [ -e "$species/${species_name}_paired2.fq" ]; then
            continue
        else
            echo $species_name
            art_illumina -i $file $species/non_redundant.fa  -m 200 -s 10 -l 150 -ss HS25 -f $SEQCOUNT -q -na -rs 7 -o  "$species/${species_name}_paired"
        fi
    done
done

#Alter FASTQ headers to include species name
for genus in $(cut -f1 "$GENUS_NAMES"); do
    for species in $(cut -f1 "$OUTPUT_PATH/$genus/output_path.txt"); do
        species_name=$(basename "$species")
        if [ -e "$species/${species_name}_paired2.fq_renamed_headers" ]; then
           continue
        else:
            for file in "$species"/*fq; do
                var=$(basename "$file")
                var+="^"
                var+=$(basename "$species")
                awk -v var="$var" '/@/ {
                    sub("@", "&"var"<")
                    sub(/\fq/, x)
                }1' "$file" > "${file%}_renamed_headers"
            done
         fi
    done
done

#Sample SIM_NUM reads from FASTQ files 
for i in {1..2}; do
    for genus in $(cut -f1 "$GENUS_NAMES"); do
        for species in $(cut -f1 "$OUTPUT_PATH/$genus/output_path.txt"); do
            species_name=$(basename "$species")
            if [ -e "$species/${species_name}_sample_2.fq" ]; then
               continue
            else
               seqtk sample -s100 $species/*_paired${i}.fq_renamed_headers $SIM_NUM  > "$species/${species_name}_sample_${i}.fq"
            fi
        done
    done
done
