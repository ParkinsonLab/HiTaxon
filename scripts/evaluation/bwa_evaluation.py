import argparse
import pandas as pd
import os

from ete3 import NCBITaxa

from HiTaxon.evaluation_utils import fasta2bwa, evaluation_bwa, ensemble

"""
Generate Ensemble Predictions
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("specialized_path", type = str, help = "path in which bwa indices are stored")
    parser.add_argument("report_name", type = str, help = "file name of output")
    parser.add_argument("report_path", type = str, help = "path to store classifer output")
    parser.add_argument("sequence_file", type = str, help = "file path of FASTA file to analyze")
    parser.add_argument("mode", type = str, help = "The ensemble mode")
    args = parser.parse_args()

    report_path = args.report_path
    report_name = args.report_name
    sequence_file = args.sequence_file
    specialized_path = args.specialized_path
    mode = args.mode
     
    ncbi = NCBITaxa()
    #K-merize FASTA file to be analyzed
    if not(os.path.exists("{report_path}/{report_name}_bwa.fa")):
        fasta2bwa(sequence_file, report_path, report_name)
    #Generate predictions using ML classifiers
    bwa_output = evaluation_bwa(report_path, report_name, specialized_path)
    bwa_output.to_csv(f"{report_path}/{report_name}_bwa.csv")
    #Ensemble ML predictions with Kraken2
    ensemble_output = ensemble(report_path, report_name, mode)
    ensemble_output.to_csv(f"{report_path}/{report_name}_ensemble_bwa.csv")

if __name__ == "__main__":
    main()
         
