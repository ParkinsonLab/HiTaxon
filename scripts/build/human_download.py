import argparse
import pandas as pd
import os
from HiTaxon.build_utils import human_genome_adder
from ete3 import NCBITaxa

"""
Create FASTA compatible with Kraken2

"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which RefSeq data is collected and stored")
    args = parser.parse_args()

    output_path = args.output_path

    human_genome_adder(output_path)

if __name__ == '__main__':
    main()