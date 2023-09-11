import argparse
import numpy as np
import pandas as pd
from ete3 import NCBITaxa

from eval_utils import standard_lcl_prediction,standard_lcl_evaluation

""" 
Evaluate FASTA with LCL Method 
"""

def main():
    ncbi = NCBITaxa()
    parser = argparse.ArgumentParser()
    parser.add_argument("model_path", type = str, help = "path in which models are saved")
    parser.add_argument("report_path", type = str, help = "path in which outputs are saved")
    parser.add_argument("report_name", type = str, help = "file name of output")
    parser.add_argument("sequence_file", type = str, help = "file path for K-merized sequence file to analyze")

    args = parser.parse_args()
    model_path = args.model_path
    report_path = args.report_path
    report_name = args.report_name
    sequence_file = args.sequence_file

    #Generate predictions 
    output = standard_lcl_prediction(sequence_file, report_path, report_name, model_path)
    output.to_csv(f"{report_path}/{report_name}_lcl.csv")

    #Run LCL inference scheme
    standard_lcl_evaluation(f"{report_name}", report_path, report_name)

if __name__ == "__main__":
    main()
