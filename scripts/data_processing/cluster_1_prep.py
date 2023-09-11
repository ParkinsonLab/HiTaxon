import argparse
import os
from HiTaxon.processing_utils import clean, query

"""
ANI Preparation
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which RefSeq sequence data is to be collected and processed")
    args = parser.parse_args()

    output_path = args.output_path
    f1 = open(f"{output_path}/missing_species.txt", "w")
    all_species = open(f'{output_path}/species_record.txt').read().splitlines()
    for species in all_species:
        species_name = species.replace(" ","_")
        genus = species_name.split("_")[0]
        #if metadata available, remove plasmid sequences from assemblies
        clean(species_name, output_path)
        #create FastANI query files composed of downloaded assemblies
        query(species_name, output_path)
        if os.stat(f'{output_path}/{genus}/{species_name}/query.txt').st_size == 0:
            f1.write(species + "\n")
    f1.close()
    

if __name__ == '__main__':
    main()
