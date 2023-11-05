source config.file

# Initialize variables to false
ACTION_C=false
ACTION_P=false
ACTION_B=false
ACTION_A=false
ACTION_T=false
ACTION_E=false
ACTION_F=false
ACTION_O=false
ACTION_M=false
HELP=false

# Listing HiTaxon's command-line arguments
usage() {
    echo "Usage: $0 [-c|--collect] [-p|--process] [-b|--build] [-t|--train] [-e|--evaluate] [--help]"
    echo "Options:"
    echo "  -c, --collect     Download assemblies from RefSeq"
    echo "  -p, --process     Perform data reduction of RefSeq sequences"
    echo "  -b, --build       Construct Kraken2 Database"
    echo "  -a, --align       Generate indices for BWA"
    echo "  -t, --train       Train FastText machine learning model"
    echo "  -e, --evaluate    Evaluate FASTA file with either Kraken2 or Ensemble"
    echo "  -f, --fasta       Path to FASTA file"
    echo "  -o, --output      Name of output report"
    echo "  -m, --mode        Evaluation mode, options being Kraken2, Kraken2_ML, Kraken2_BWA"
    echo "  --help            Display this help message"
    exit 1
}

# Parse command-line arguments 
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--collect) 
            if [ "$ACTION_P" = true ] || [ "$ACTION_B" = true ] || [ "$ACTION_A" = true ] || [ "$ACTION_T" = true ] || [ "$ACTION_E" = true ] || [ "$ACTION_F" = true ] || [ "$ACTION_O" = true ] || [ "$ACTION_M" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            ACTION_C=true
            shift
            ;;
        -p|--process)
            if [ "$ACTION_C" = true ] || [ "$ACTION_B" = true ] || [ "$ACTION_A" = true ] || [ "$ACTION_T" = true ] || [ "$ACTION_E" = true ] || [ "$ACTION_F" = true ] || [ "$ACTION_O" = true ] || [ "$ACTION_M" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            ACTION_P=true
            shift
            ;;
        -b|--build)
            if [ "$ACTION_C" = true ] || [ "$ACTION_P" = true ] || [ "$ACTION_A" = true ] || [ "$ACTION_T" = true ] || [ "$ACTION_E" = true ] || [ "$ACTION_F" = true ] || [ "$ACTION_O" = true ] || [ "$ACTION_M" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            ACTION_B=true
            shift
            ;;
        -a|--align)
            if [ "$ACTION_C" = true ] || [ "$ACTION_P" = true ] || [ "$ACTION_B" = true ] || [ "$ACTION_T" = true ] || [ "$ACTION_E" = true ] || [ "$ACTION_F" = true ] || [ "$ACTION_O" = true ] || [ "$ACTION_M" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            ACTION_A=true
            shift
            ;;
        -t|--train)
            if [ "$ACTION_P" = true ] || [ "$ACTION_B" = true ] || [ "$ACTION_A" = true ] || [ "$ACTION_C" = true ] || [ "$ACTION_E" = true ] || [ "$ACTION_F" = true ] || [ "$ACTION_O" = true ] || [ "$ACTION_M" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            ACTION_T=true
            shift
            ;;
        -e|--evaluate)
            if [ "$ACTION_P" = true ] || [ "$ACTION_B" = true ] || [ "$ACTION_A" = true ] || [ "$ACTION_C" = true ] || [ "$ACTION_T" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            ACTION_E=true
            shift
            ;;
        -f|--fasta)
            if [ "$ACTION_E" = true ]; then
                SEQUENCE_FILE="$2"
                shift 2
            else
                echo "Option -f can only be used with -e."
                usage
            fi
            ;;
        -o|--output)
            if [ "$ACTION_E" = true ]; then
                REPORT_NAME="$2"
                shift 2
            else
                echo "Option -o can only be used with -e."
                usage
            fi
            ;;
        -m|--mode)
            if [ "$ACTION_E" = true ]; then
                MODE="$2"
                if [ "$MODE" != "Kraken2" ] && [ "$MODE" != "Kraken2_ML" ] && [ "$MODE" != "Kraken2_BWA" ]; then
                    echo "Invalid option for the -m flag. Please use 'Kraken2', 'Kraken2_ML', or Kraken2_BWA."
                    usage
                    exit 1
                fi
                shift 2
            else
                echo "Option -m can only be used with -e."
                usage
            fi
            ;;
        --help)
            if [ "$ACTION_P" = true ] || [ "$ACTION_B" = true ] || [ "$ACTION_C" = true ] || [ "$ACTION_A" = true ] || [ "$ACTION_T" = true ] || [ "$ACTION_E" = true ] || [ "$ACTION_F" = true ] || [ "$ACTION_O" = true ] || [ "$ACTION_M" = true ] || [ "$HELP" = true ]; then
                echo "Options -c, -p, -b, -a, -t, -e and --help cannot be used together."
                usage
            fi
            HELP=true
            usage
            ;;
        *)
            usage
            ;;
    esac
done

# If command-line argument given, run appropriate script
if [ "$ACTION_C" = true ]; then
    echo "Collecting assemblies"
    ./scripts/collection.sh $ASSEMBLY_SUMMARY $GENUS_NAMES $OUTPUT_PATH
fi

if [ "$ACTION_P" = true ]; then
    echo "Processing Sequences"
    ./scripts/processing.sh $GENUS_NAMES $NUM_OF_THREADS $OUTPUT_PATH 
fi

if [ "$ACTION_B" = true ]; then
    echo "Building Kraken2 Database"
    ./scripts/build.sh $ASSEMBLY_SUMMARY $GENUS_NAMES $KRAKEN_NAME $KRAKEN_PATH $NUM_OF_THREADS $OUTPUT_PATH
fi

if [ "$ACTION_A" = true ]; then
    echo "Generating BWA indices"
    ./scripts/align.sh $GENUS_NAMES $BWA_PATH $OUTPUT_PATH 
fi

if [ "$ACTION_T" = true ]; then
    echo "Training machine learning classifiers"
    ./scripts/train.sh $GENUS_NAMES $MODEL_PATH $OUTPUT_PATH  
fi

if [ "$ACTION_E" == true ]; then
    if [ -n "$SEQUENCE_FILE" ] && [ -n "$MODE" ] && [ -n "$REPORT_NAME" ]; then
        echo "Evaluating sequence data"
        if [ "$ASSEMBLY_SUMMARY" = "default" ]; then
            ASSEMBLY_SUMMARY="$OUTPUT_PATH/assembly_summary.txt"
        fi
        if [ "$MODE" = "Kraken2_ML" ]; then
            ./scripts/evaluation.sh $KRAKEN_NAME $KRAKEN_PATH $MODE $MODEL_PATH $NUM_OF_THREADS $ASSEMBLY_SUMMARY $REPORT_NAME $REPORT_PATH $SEQUENCE_FILE 
        else
            ./scripts/evaluation.sh $KRAKEN_NAME $KRAKEN_PATH $MODE $BWA_PATH $NUM_OF_THREADS $ASSEMBLY_SUMMARY $REPORT_NAME $REPORT_PATH $SEQUENCE_FILE
        fi
    else
        echo "-f, -m, -o are required."
    fi
fi

if [ "$ACTION_C" != true ] && [ "$ACTION_P" != true ] && [ "$ACTION_B" != true ] && [ "$ACTION_T" != true ]  && [ "$ACTION_E" != true ]  && [ "$HELP" != true ]; then
    usage
fi
                                                                      
