#!/usr/bin/env bash
ASSEMBLY_SUMMARY=$1
GENUS_NAMES=$2
KRAKEN_NAME=$3
KRAKEN_PATH=$4
NUM_OF_THREADS=$5
OUTPUT_PATH=$6

#Create directory and install taxonomy
if [ ! -d "$KRAKEN_PATH/$KRAKEN_NAME" ]; then
    echo "Kraken DB directory not found. Creating Directory and Downloading Taxonomy..."
    "scripts/build/kraken_download.sh" $KRAKEN_NAME $KRAKEN_PATH 
    echo "Taxonomy installed"
    # Ask if want to quit, as remaining scripts can be used offline (I.E SciNet Niagara compute nodes)
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

#If no RefSeq assembly provided, use newest assembly installed during collection
if [ "$ASSEMBLY_SUMMARY" = "default" ]; then
    ASSEMBLY_SUMMARY="$OUTPUT_PATH/assembly_summary.txt"
fi

#Prepare sequences for Kraken2
python "scripts/build/kraken_build.py" $ASSEMBLY_SUMMARY $KRAKEN_PATH $OUTPUT_PATH 

#Add sequences to Kraken2 library and build database
"scripts/build/kraken_build.sh" $GENUS_NAMES $KRAKEN_NAME $KRAKEN_PATH $NUM_OF_THREADS $OUTPUT_PATH
