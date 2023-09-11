import argparse
import pandas as pd
import os

from ete3 import NCBITaxa

from HiTaxon.evaluation_utils import fasta2kmer, evaluation, ensemble

"""
Generate Ensemble Predictions
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_path", type = str, help = "path in which models are stored")
    parser.add_argument("report_name", type = str, help = "file name of output")
    parser.add_argument("report_path", type = str, help = "path to store classifer output")
    parser.add_argument("sequence_file", type = str, help = "file path of FASTA file to analyze")
    args = parser.parse_args()

    report_path = args.report_path
    report_name = args.report_name
    sequence_file = args.sequence_file
    model_path = args.model_path

    ncbi = NCBITaxa()
    #K-merize FASTA file to be analyzed
    if not(os.path.exists("{report_path}/{report_name}_kmer.txt")):
        fasta2kmer(sequence_file, report_path, report_name)
    #Generate predictions using ML classifiers
    ml_output = evaluation(report_path, report_name, model_path, 0.5)
    ml_output.to_csv(f"{report_path}/{report_name}_ml.csv")
    #Ensemble ML predictions with Kraken2
    ensemble_output = ensemble(report_path, report_name)
    ensemble_output.to_csv(f"{report_path}/{report_name}_ensemble.csv")

if __name__ == "__main__":
    main()
         
