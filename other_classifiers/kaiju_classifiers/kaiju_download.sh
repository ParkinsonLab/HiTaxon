#!/usr/bin/env bash
KAIJU_PATH=$1 #Path in which kaiju database is to built

mkdir $KAIJU_PATH
#Download and unzip necessary files to create custom Kaiju database
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip -P $KAIJU_PATH
unzip $KAIJU_PATH/new_taxdump.zip -d $KAIJU_PATH


