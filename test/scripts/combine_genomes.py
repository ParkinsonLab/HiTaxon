import argparse
import numpy as np
import os
import pandas as pd


from Bio import SeqIO

def combine_reads(test_instance, test_path):
    """
    Create FASTA files composed of simulated reads from different species
    Args:
        test_instance: simulated test instance (1...10)
        test_path: path of test directory
    """
    #List individual species directories
    directories = (os.listdir(f"{test_path}/simulated_reads/{test_instance}"))
    #Default files Polyester simulations
    file1 = "sample_02_1.fasta"
    file2 = "sample_02_2.fasta"
    #file paths for FASTA files to be created
    fasta1 = f"{test_path}/simulated_reads/{test_instance}/merged1.fa"
    fasta2 = f"{test_path}/simulated_reads/{test_instance}/merged2.fa"

    #Aggregrate paired-end reads from each simulated species to appropriate FASTA file
    fasta_files = [(fasta1, file1), (fasta2, file2)]
    for fasta in fasta_files:
            f = open(f"{fasta[0]}", "w")
            #indiv_directory == indiv_species_directory
            for indiv_directory in directories:
                       	if not(os.path.exists(f"{test_path}/simulated_reads/{test_instance}/{indiv_directory}/{fasta[1]}")):
                                continue
                        for record in SeqIO.parse(f"{test_path}/simulated_reads/{test_instance}/{indiv_directory}/{fasta[1]}", "fasta"):
                                f.write(f">{indiv_directory}" + "\n" +  str(record.seq) + "\n")
            f.close()

"""
Aggregate simulated paired-ends from multiple organisms to create metagenomics/metatranscriptomic test set

"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("test_path", type = str, help = "path of test directory")
    args = parser.parse_args()

    test_path = args.test_path

    #For test instance 1 to 10, aggregate reads into independent FASTA files
    for sample in range(1, 11):
        combine_reads(sample, test_path,)

if __name__ == "__main__":
    main()
