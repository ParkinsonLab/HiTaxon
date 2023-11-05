#!/usr/bin/env bash
KRAKEN_NAME=$1
KRAKEN_PATH=$2
MODE=$3
SPECIALIZED_PATH=$4
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

if [ ! -e $ASSEMBLY_SUMMARY ]; then
    wget "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" -O $ASSEMBLY_SUMMARY
    while true; do
        read -p "Do you want to quit? (yes/no): " answer
        case $answer in
            [Yy]|[Yy][Ee][Ss])
                echo "Goodbye!"
                exit
                ;;
            [Nn]|[Nn][Oo])
                break
                ;;
            *)
                echo "Invalid input. Please enter 'yes' or 'no'."
                ;;
        esac
    done
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
#Ensemble Kraken2's output with ML classifiers
elif [ "$MODE" = "Kraken2_ML" ]; then
    echo "MODE is set to Kraken2_ML"
    python "scripts/evaluation/fasttext_evaluation.py" $SPECIALIZED_PATH $REPORT_NAME $REPORT_PATH $SEQUENCE_FILE

else
    echo "MODE is set to Kraken2_BWA"
    python "scripts/evaluation/bwa_evaluation.py" $SPECIALIZED_PATH $REPORT_NAME $REPORT_PATH $SEQUENCE_FILE $MODE
fi


