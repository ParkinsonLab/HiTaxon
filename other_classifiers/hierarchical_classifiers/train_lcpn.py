import argparse
import numpy as np
import os
from ete3 import NCBITaxa
from train_utils import lineage_creator, genus_select, create_lcpn_data, create_lcl_data, model_train

""" 
Train LCPN Classifier
"""

def main():
    ncbi = NCBITaxa()
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", type = str, help = "path in which RefSeq sequences are downloaded and processed ")
    parser.add_argument("genus_names", type = str, help = "text file with all genera")
    parser.add_argument("model_path", type = str, help = "path in which models are saved")

    args = parser.parse_args()

    output_path = args.output_path
    model_path = args.model_path
    genus_names = open(f"{output_path}/{args.genus_names}").read().splitlines()

    #Create dataframe with lineage for each genus
    lineages = lineage_creator(genus_names, ncbi)
    ranks = ["genus", "family", "class", "order", "phylum"]
    #Exclude creating Genus -> Species classifiers, as can be done using ./HiTaxon -t
    for rank in ranks[1:]:
        current_idx = ranks.index(rank)
        for taxa in np.unique(lineages[rank].values):
            lineage_subset = lineages[lineages[rank] == taxa]
            #Check if more then one unique descedent for this taxa
            if len(np.unique(lineage_subset[ranks[current_idx-1]].values)) > 1:
                mode = "multi"
            else:
                mode = "binary"
            selected_genera = genus_select(taxa, rank, lineages, mode)
            #This will occur if an outgroup is needed
            if selected_genera == None:
               exit()
            #Create training data
            create_lcpn_data(taxa, selected_genera, ranks[current_idx-1], lineages, output_path, model_path)
            #Shuffle training data
            os.system(f"shuf {model_path}/{taxa}_train.txt > {model_path}/{taxa}_train_shuffled.txt")
            
            model = model_train(taxa, model_path)
            model.save_model(f"{model_path}/{taxa}_model.bin")
    #Need to create a multi-class Phylum model to initilize the LCPN architecture, is the same classifier used in LCL method
    create_lcl_data(taxa_of_interest, "phylum", lineages, output_path, model_path)
    os.system(f"shuf {model_path}/phylum_train.txt > {model_path}/phylum_train_shuffled.txt")
    model = model_train("phylum", model_path)
    model.save_model(f"{model_path}/phylum_model.bin")

if __name__ == '__main__':
    main()
