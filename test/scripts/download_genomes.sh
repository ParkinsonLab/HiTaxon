#!/usr/bin/env bash
TESTDIR=$1 

# Download coding seqeunces for each assembly in the Test Assembly CSV file
for genome in $(cut -d ',' -f 2 "$TESTDIR/test_genomes.csv"); do
    datasets download genome accession $genome --include cds --filename "$TESTDIR/${genome}.zip"
    unzip "$TESTDIR/$genome.zip" -d "$TESTDIR/$genome"
done
