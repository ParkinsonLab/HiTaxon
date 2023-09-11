#!/usr/bin/env bash
GENUS_NAMES=$1
MODEL_PATH=$2
OUTPUT_PATH=$3

#Create directory to store model path if not createed
if [ ! -d "$MODEL_PATH" ]; then
  mkdir -p "$MODEL_PATH"
  echo "Directory created: $MODEL_PATH"
else
  echo "Directory already exists: $MODEL_PATH"
fi

#Simulate paired-end reads from non-redundant sequences
"scripts/train/simulator.sh" $GENUS_NAMES $OUTPUT_PATH 

#Process simulated reads to prepare for training
"scripts/train/prepare.sh" $GENUS_NAMES $OUTPUT_PATH 

#Train specie-level ML classifiers
python "scripts/train/train.py" $GENUS_NAMES $MODEL_PATH $OUTPUT_PATH
