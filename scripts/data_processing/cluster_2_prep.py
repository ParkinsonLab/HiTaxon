import argparse
import os
import pandas as pd
import shutil

from HiTaxon.processing_utils import cluster_default, cluster_alternative, find_representative, aggregator
"""
Create species-specific sets of represenative coding sequences
"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which sequence data is collected and processed")
    args = parser.parse_args()

    output_path = args.output_path

    all_species = open(f'{output_path}/species_record.txt').read().splitlines()
    missing_species = open(f'{output_path}/missing_species.txt').read().splitlines()
    for species in all_species:
        if species in missing_species:
            continue
        species = species.replace(" ","_")
        genus = species.split("_")[0]
        if os.path.exists(f"{output_path}/{genus}/{species}/rep_assemblies.txt"):
            continue
        if os.path.exists(f"{output_path}/{genus}/{species}/ani_output.matrix"):
            #Generate Clusters of similiar assemblies
            cluster_output, assemblies = cluster_default(species, output_path, 400)
            
            #For instances in which species had one assembly in RefSeq
            if(cluster_output == None): 
                f = open(f"{output_path}/{genus}/{species}/rep_assemblies.txt", "w")
                f.write(assemblies + "\n")
                f.close()
            #Find represenative assembly from each cluster
            else:
                clustering_results = pd.DataFrame({'cluster': cluster_output.labels_,'assemblies': assemblies})
                clustering_results.to_csv(f"{output_path}/{genus}/{species}/clustered.csv")
                representative_assemblies = find_representative(species, cluster_output, assemblies, output_path)
                f = open(f"{output_path}/{genus}/{species}/rep_assemblies.txt", "w")
                for assembly_path in representative_assemblies:
                    f.write(assembly_path + "\n")
                f.close()
        else:
            #Instances in which no ANI output was generated (i.e all assemblies had less then 80% ANI)
            shutil.copyfile(f"{output_path}/{genus}/{species}/query.txt", f"{output_path}/{genus}/{species}/rep_assemblies.txt")
        #Combine coding sequences from represenative assemblies into single file
        aggregator(species, output_path)


if __name__ == '__main__':
    main()
