import argparse
import numpy as np
import fasttext as ft
import os

from HiTaxon.train_utils import build_kmers, data_creation, negative_set, model_train

""" 
Train species-level machine learnring classifiers 
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genus_names", type = str, help = "path to text file which lists the genera of interest, seperated by newline characters")
    parser.add_argument("model_path", type = str, help = "path in which machine learning classifiers are saved")
    parser.add_argument("output_path", type = str, help = "path in which RefSeq seqeunces are collected and processed" )
    args = parser.parse_args()

    output_path = args.output_path
    model_path = args.model_path

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    genus_names = open(args.genus_names).read().splitlines()
    sorted_genus = []
    #Train models from smallest to largest
    for genus in genus_names:
        if not(os.path.exists(f"{output_path}/{genus}/{genus}_shuffled.fa")):
            continue
        else:
            sorted_genus.append((genus, os.path.getsize(f"{output_path}/{genus}/{genus}_shuffled.fa")))
    sorted_genus = sorted(sorted_genus, key=lambda x: x[1])

    for genus, size in sorted_genus:
        if os.path.exists(f"{model_path}/{genus}_model.bin"):
            continue
        #K-merize and shuffle training data
        if not(os.path.exists(f"{model_path}/{genus}_train_shuffled.txt")):
            data_creation(genus, 13, genus_names, ncbi, output_path, model_path)
            os.system(f"shuf {model_path}/{genus}_train.txt > {model_path}/{genus}_train_shuffled.txt")
        #Train models
        model = model_train(genus, model_path)
        model.save_model(f"{model_path}/{genus}_model.bin")
        os.remove(f"{model_path}/{genus}_train.txt")
        os.remove(f"{model_path}/{genus}_train_shuffled.txt")




if __name__ == '__main__':
    main()
