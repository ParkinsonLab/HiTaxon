#!/usr/bin/env bash
KRAKEN_NAME=$1
KRAKEN_PATH=$2
MODE=$3
MODEL_PATH=$4
NUM_OF_THREADS=$5
ASSEMBLY_SUMMARY=$6
REPORT_NAME=$7
REPORT_PATH=$8
SEQUENCE_FILE=$9

echo $REPORT_PATH
echo $SEQUENCE_FILE

#Create directory to store taxonomic predictions if not createed
if [ ! -d "$REPORT_PATH" ]; then
  mkdir -p "$REPORT_PATH"
  echo "Directory created: $REPORT_PATH"
else
  echo "Directory already exists: $REPORT_PATH"
fi


#Generate predictions with Kraken2
if [ -e "${REPORT_PATH}/${REPORT_NAME}_lineage_kraken.csv" ]; then
    echo "Kraken Processed File Exist"
else
    echo $SEQUENCE_FILE
    scripts/evaluation/kraken_evaluation.sh $KRAKEN_NAME $KRAKEN_PATH $NUM_OF_THREADS $REPORT_NAME $REPORT_PATH $SEQUENCE_FILE 
    #Get entire lineage of predictions
    python "scripts/evaluation/kraken_evaluation.py" $ASSEMBLY_SUMMARY $REPORT_NAME $REPORT_PATH 
fi

if [ "$MODE" = "Kraken2" ]; then
    echo "MODE is set to Kraken2. Exiting the script."
    exit 1
fi

#Ensemble Kraken2's output with ML classifiers
python "scripts/evaluation/fasttext_evaluation.py" $MODEL_PATH $REPORT_NAME $REPORT_PATH $SEQUENCE_FILE
