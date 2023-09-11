import argparse
from Bio import SeqIO
from ete3 import NCBITaxa

def kaiju_FASTA_generator(species, output_path, ncbi):
    """
    Reformat Prodigal FASTA for a single species to appropriate format for Kaiju custom database
    Args:
        species = species of interest
        output_path = Path in which RefSeq sequences were downloaded and processed is stored; including prodigal sequences 
        ncbi = NCBITaxa()

    """
    genus = species.split("_")[0]
    species_name = species.replace("_"," ")
    name2taxid = ncbi.get_name_translator([species_name])[species_name][0] 
    kaiju_fasta = open(f"{output_path}/{genus}/{species}/kaiju.fa", "w")
    counter = 0
    for record in SeqIO.parse(f"{output_path}/{genus}/{species}/protein.fa", "fasta"):
        kaiju_fasta.write(f">{counter}_{name2taxid}" + "\n" + str(record.seq).replace("*","") + "\n")
        counter +=1
    kaiju_fasta.close()

"""
Create Kaiju compatible FASTA files from prodigal outputs
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "Path in which RefSeq sequences were downloaded and processed is stored; including prodigal sequences")
    args = parser.parse_args()

    output_path = args.output_path

    #Get list of all species data was collected and processed for
    all_species = open(f'{output_path}/species_record.txt').read().splitlines()
    missing_species = open(f'{output_path}/species_record.txt').read().splitlines()
    ncbi = NCBITaxa()
    for species in all_species:
        if species in missing_species:
            continue
        species = species.replace(" ","_")
        kaiju_FASTA_generator(species, output_path, ncbi)

if __name__ == '__main__':
    main()
