import argparse
import numpy as np
import pandas as pd
from ete3 import NCBITaxa

from eval_utils import lcl_lcpn_prediction, lcl_lcpn_evaluation

""" 
Evaluate FASTA with LCL-LCPN Method 
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

    #Generate predictions
    output = lcl_lcpn_prediction(args.sequence_file, report_path, report_name, model_path)
    output.to_csv(f"{report_path}/{report_name}_lcl_lcpn.csv")

    #Run LCL-LCPN inference scheme  
    lcl_lcpn_evaluation(f"{report_name}", report_path, report_name)

if __name__ == "__main__":
    main()
