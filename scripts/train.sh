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

# Initialize the mode variable
mode=""

# Function to show prompt and read user input
ask_mode() {
  echo "Please enter the mode you want to be in:"
  echo "1) data_creation"
  echo "2) model_creation"
  read -p "Enter your choice (1 or 2): " choice
}

# Loop until a valid input is provided
while true; do
  ask_mode

  case $choice in
    1|data_creation)
      mode="data_creation"
      echo "You have chosen 'data_creation' mode."
      break
      ;;
    2|model_creation)
      mode="model_creation"
      echo "You have chosen 'model_creation' mode."
      break
      ;;
    *)
      echo "Invalid input. Please enter 1 for 'data_creation' or 2 for 'model_creation'."
      ;;
  esac
done

#Train species-level ML classifiers
python "scripts/train/train.py" $GENUS_NAMES $MODEL_PATH $OUTPUT_PATH $mode
