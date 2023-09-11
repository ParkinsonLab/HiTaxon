#!/usr/bin/env bash
GENUS_NAMES=$1
OUTPUT_PATH=$2



while read -r genus; do
  if [ -e "$OUTPUT_PATH/$genus/${genus}_shuffled.fa" ]; then
      continue
  else
      while read -r species; do
         # Use fastq-join to merge the sample fastq files
         species_name=$(basename "$species")
         fastq-join "$species/${species_name}_sample_1.fq" "$species/${species_name}_sample_2.fq" -o "$species/${species_name}_merged.fq"

         # Convert the merged fastq to fasta format using seqtk
         seqtk seq -a "$species/${species_name}_merged.fqjoin" > "$species/${species_name}_merged.fa"

         # Concatenate the species merged fasta files into a genus-level merged file
         cat "$species/${species_name}_merged.fa" >> "$OUTPUT_PATH/$genus/${genus}_merged.fa"
         done < "$OUTPUT_PATH/$genus/output_path.txt"
  fi
  cat "$OUTPUT_PATH/$genus/${genus}_merged.fa" | perl scripts/train/seq-shuf > "$OUTPUT_PATH/$genus/${genus}_shuffled.fa"
done < "$GENUS_NAMES"

#shuffling script can be found at https://github.com/thackl/seq-scripts"
