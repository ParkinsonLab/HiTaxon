#!/usr/bin/env bash
KRAKEN_NAME=$1
KRAKEN_PATH=$2 


mkdir $KRAKEN_PATH
#Create taxonomy
kraken2-build --download-taxonomy --db $KRAKEN_PATH/$KRAKEN_NAME
