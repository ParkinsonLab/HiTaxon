import argparse
import numpy as np
import pandas as pd


from Bio import SeqIO
from ete3 import NCBITaxa

def species_ground_truth(sample, test_path):
    """
    Return species classification for each simulated read
    Args:
        sample: simulated test instance (1...10)
        test_path: path of test directory
    """
    ground_truth = []
    #append the FASTA header as the have the species name
    for record in SeqIO.parse(f"{test_path}/simulated_reads/{sample}/merged.fa", "fasta"):
            ground_truth.append(str(record.id))
    return ground_truth


def expand_lineage(taxa_of_interest, ncbi):
    """
    Determine the phylum to genus taxonomic lineage of a species
    Args:
         taxa_of_interest: species of interest
         ncbi: NCBITaxa()
    """
    taxa_of_interest = taxa_of_interest.replace("_"," ")

    #Determine taxid and return entire list of lineage in taxids
    name2taxid = ncbi.get_name_translator([taxa_of_interest])[taxa_of_interest][0]
    lineage = ncbi.get_lineage(name2taxid) 

    #Subset for specific ranks  
    named_lineage = {"species": "NA", "genus": "NA", "family": "NA", "order": "NA", "class": "NA", "phylum": "NA"}
    if lineage is None:
        return named_lineage
    for taxa in lineage:
        rank = ncbi.get_rank([taxa])[taxa]
        if rank in named_lineage.keys():
            #Convert taxid to name
            named_lineage[rank] =  ncbi.get_taxid_translator([taxa])[taxa]
    return named_lineage


def expand_truth(ground_truth, ncbi):
    """
    Return entire taxonomic lineage (Phylum to Species) of each read in simulated test set
    Args:
        ground_truth: list of species pertaining to each read in the simulated test set
        ncbi: NCBITaxa()

    """
    #Determine lineage for each unique species in the test set
    unique_taxa = np.unique(np.array(ground_truth))
    lineages = {}
    for taxa in unique_taxa:
        lineages[taxa] = expand_lineage(taxa, ncbi)
    #For each read, define its entire lineage
    truth_expanded = [lineages[taxa] for taxa in ground_truth]
    truth_expanded = pd.DataFrame(truth_expanded)
    return truth_expanded

"""
Ground truth values of test reads for entire lineage (Species to Phylum)
"""

def main():

    ncbi = NCBITaxa()
    parser = argparse.ArgumentParser()
    parser.add_argument("test_path", type = str, help = "path of test directory")
    args = parser.parse_args()

    test_path = args.test_path

    #For test instance 1 to 10, determine the taxonomic lineage for each read
    for sample in range(1, 11):
        ground_truth = species_ground_truth(sample, test_path)
        expanded_ground_truth = expand_truth(ground_truth, ncbi)
        expanded_ground_truth.to_csv(f"{test_path}/ground_truth_{sample}.csv")

if __name__ == "__main__":
    main()
