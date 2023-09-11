import argparse
import numpy as np
import pandas as pd
from ete3 import NCBITaxa

from eval_utils import hierarchical_lcl_prediction, hierarchical_lcl_evaluation
from train_utils import lineage_creator

""" 
Evaluate FASTA with Hierarchy informed LCL Method 
"""

def main():
    ncbi = NCBITaxa()
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which RefSeq sequences are downloaded and processed")
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

    #Generate Predictions
    hierarchical_lcl_prediction(sequence_file, report_path, report_name, model_path)

    #Run Heirarchy-informed LCL inference scheme 
    species_record = open(f"{output_path}/{args.species_present}").read().splitlines()
    genus_record = [species.split(" ")[0] for species in species_record]
    lineages = lineage_creator(genus_record, ncbi)
    hierarchical_lcl_evaluation(lineages, output_path, species_record, report_path, report_name)

if __name__ == "__main__":
    main()
