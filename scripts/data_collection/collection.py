import argparse
import os
from HiTaxon.collection_utils import assembly_collection

"""
Collect high-quality assemblies from RefSeq
"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_summary", type = str, help = "path to NCBI RefSeq Assembly File")
    parser.add_argument("genus_names", type = str, help = "path to text file which list all genera of interest")
    parser.add_argument("output_path", type = str, help = "path in which RefSeq sequences will be collected and processed")
    args = parser.parse_args()
    
    assembly_summary = args.assembly_summary
    genus_names = open(args.genus_names).read().splitlines()
    output_path = args.output_path
    
    #Collect assemblies per species 
    collected_species, missing_genera = assembly_collection(assembly_summary, genus_names, output_path)
    
    if os.path.exists(f"{output_path}/species_record.txt"):
        species_record = open(f"{output_path}/species_record.txt").read().splitlines()
    else:
        species_record = []

    #Create record of all species assemblies have been collected for
    species_output= open(f"{output_path}/species_record.txt", "a")

    for species, assemblies in collected_species.items():
        genus = species.split(" ")[0]
        species_name = species.replace(" ","_")
        if not(species in species_record):
            species_output.write(species + "\n")
            os.makedirs(f"{output_path}/{genus}/{species_name}")
        #Create list of assemblies to download for each species
        download_output = open(f"{output_path}/{genus}/{species_name}/download_assemblies.txt", "w")
        for assembly in assemblies:
            download_output.write(assembly + "\n")
        download_output.close()
    species_output.close()

    #Create record of genus in which no assemblies were found 
    missing_genera_output = open(f"{output_path}/missing_genera.txt", "w")
    for genus in missing_genera:
        missing_genera_output.write(genus + "\n")
        print(f"Formal named species for {genus} is not present in RefSeq Assembly")
    missing_genera_output.close()


if __name__ == '__main__':
    main()

