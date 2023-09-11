import argparse
import numpy as np
import pandas as pd
from ete3 import NCBITaxa

from eval_utils import evaluate_taxa, lcpn_prediction, lcpn_evaluation
from train_utils import lineage_creator

""" 
Evaluate FASTA with LCPN Method 
"""

def main():
    ncbi = NCBITaxa()

    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which downloaded and processed assembly data is stored")
    parser.add_argument("species_present", type = str, help ="text file of species")
    parser.add_argument("model_path", type = str, help = "path in which models are saved")
    parser.add_argument("report_path", type = str, help = "path in which outputs are saved")
    parser.add_argument("report_name", type = str, help = "file name of output")
    parser.add_argument("sequence_file", type = str, help = "file path for K-merized sequence file to analyze")

    args = parser.parse_args()
    output_path = args.output_path
    model_path = args.model_path
    report_path = args.report_path
    report_name = args.report_name
    sequence_file = args.sequence_file

    #Determine lineage of each species
    species_present = open(f"{output_path}/{args.species_present}").read().splitlines()
    genus_present = [species.split(" ")[0] for species in species_present]
    lineages = lineage_creator(genus_present, ncbi)
    lineages["species"] = species_present

    #Generate predictions 
    output = lcpn_prediction(lineages, sequence_file, report_path, report_name, model_path)
    output.to_csv(f"{report_path}/{report_name}")

    #Run LCPN inference scheme
    lcpn_evaluation(f"{report_name}", report_path, report_name)

if __name__ == "__main__":
    main()
