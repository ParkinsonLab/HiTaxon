import argparse
import pandas as pd
import os
from HiTaxon.build_utils import id_generator, kraken_FASTA_generator
from ete3 import NCBITaxa

"""
Create FASTA compatible with Kraken2

"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_summary", type = str, help = "File name of NCBI Refseq assembly")
    parser.add_argument("kraken_path", type = str, help = "path in which Kraken2 database is to be stored")
    parser.add_argument("output_path", type = str, help = "path in which RefSeq data is collected and stored")
    args = parser.parse_args()

    output_path = args.output_path
    kraken_path = args.kraken_path
    assembly_summary = args.assembly_summary
    all_species = open(f'{output_path}/species_record.txt').read().splitlines()
    missing_species =  open(f'{output_path}/missing_species.txt').read().splitlines()
    ncbi = NCBITaxa()

    refseq = pd.read_csv(assembly_summary, delimiter = "\t", skiprows = 1)
    refseq = refseq[["#assembly_accession", "refseq_category", "species_taxid", "organism_name", "infraspecific_name", "assembly_level", "genome_rep"]]
    refseq["species"] = refseq["organism_name"].apply(lambda x: " ".join(x.split(" ")[:2]))
    
    #Track which species were unable to be added to Kraken2
    f = open(f"{kraken_path}/species_not_added.txt","w")

    #Write Kraken2 FASTA files with appropriate headers for each species
    for species in all_species:
        if species in missing_species:
            continue
        name2taxid = id_generator(species, refseq, ncbi)
        if (name2taxid == None):
            f.write(missing + "\n")
            continue
        else:
            species_name = species.replace(" ","_")
            genus = species_name.split("_")[0]
            if os.path.exists(f"{output_path}/{genus}/{species_name}/kraken.fa"):
                continue
            kraken_FASTA_generator(species_name, name2taxid, output_path)
    f.close()

if __name__ == '__main__':
    main()
