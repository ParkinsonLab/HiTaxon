import argparse
import numpy as npi
import fasttext as ft
import os
from ete3 import NCBITaxa
from train_utils import lineage_creator, genus_select, create_lcpn_data, create_lcl_data, model_train

""" 
Train LCL Classifier
"""
def main():
    ncbi = NCBITaxa()
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which RefSeq sequences are downloaded and processed")
    parser.add_argument("genus_names", type = str, help = "text file with all genera")
    parser.add_argument("model_path", type = str, help = "path in which models are saved")

    args = parser.parse_args()

    output_path = args.output_path
    model_path = args.model_path
    genus_names = open(f"{output_path}/{args.genus_names}").read().splitlines()

    #Create DataFrame with lineage for each genus
    lineages = lineage_creator(genus_names, ncbi)
    ranks = ["species","genus", "family", "class", "order", "phylum"]
    for rank in ranks:
        #Create training data
        create_lcl_data(taxa_of_interest, rank, lineages, output_path, model_path)
        #Shuffle training data
        os.system(f"shuf {model_path}/{rank}_train.txt > {model_path}/{rank}_train_shuffled.txt")
        
        model = model_train(rank, model_path)
        model.save_model(f"{model_path}/{rank}_model.bin")
if __name__ == '__main__':
    main()
