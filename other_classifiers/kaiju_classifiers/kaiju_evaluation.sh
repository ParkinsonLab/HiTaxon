#!/usr/bin/env bash
SEQUENCE_FILE=$1 #Path to FASTA file to be analyzed by Kaiju
KAIJU_PATH=$2 
REPORT_PATH=$3 
REPORT_NAME=$4 
NUM_OF_THREADS=$5 

#Run Kaiju
kaiju -z $NUM_OF_THREADS -t $KAIJU_PATH/nodes.dmp -f $KAIJU_PATH/proteins.fmi -i ${SEQUENCE_FILE} -o $REPORT_PATH/${REPORT_NAME}_kaiju.out 
#Add taxa names to Kaiju Output
kaiju-addTaxonNames -r phylum,class,order,family,genus,species -t $KAIJU_PATH/nodes.dmp -n $KAIJU_PATH/names.dmp -i $REPORT_PATH/${REPORT_NAME}_kaiju.out -o $REPORT_PATH/${REPORT_NAME}_kaiju_taxa.out
