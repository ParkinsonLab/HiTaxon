import argparse
import numpy as np
import pandas as pd
import os
from HiTaxon.align_utils import bwa_FASTA_generator


"""
Create FASTA for BWA indices

"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bwa_path", type = str, help = "path in which bwa indices are to be stored")
    parser.add_argument("output_path", type = str, help = "path in which RefSeq data is collected and stored")
    parser.add_argument("taxa_path", type = str, help = "path in pertaining to relevant genera")
    args = parser.parse_args()

    output_path = args.output_path
    bwa_path = args.bwa_path
    all_species = open(f'{output_path}/species_record.txt').read().splitlines()
    missing_species =  open(f'{output_path}/missing_species.txt').read().splitlines()
    
    #Track which species were unable to be added to BWA indices
    f = open(f"{bwa_path}/species_not_added.txt","w")
    
    #Track which BWA indices need to be made
    f2 = open(f"{bwa_path}/genus2add.txt","w")
    
    relevant_genus = open(args.taxa_path).read().splitlines()
    genus2add = []
    #Write BWA FASTA files with appropriate headers for each species
    for species in all_species:
        genus = species.split(" ")[0]
        if species in missing_species or genus not in relevant_genus:
            continue
        species_name = species.replace(" ","_")
        genus = species_name.split("_")[0]
        if os.path.exists(f"{bwa_path}/{genus}.fa.ann"):
            continue
        genus2add.append(genus)
        bwa_FASTA_generator(species_name, output_path, bwa_path)
    f.close()
    for genus in np.unique(genus2add):
        f2.write(genus + "\n")
    f2.close()

if __name__ == '__main__':
    main()
