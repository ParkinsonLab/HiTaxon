import numpy as np
import fasttext as ft
import pandas as pd

from Bio import SeqIO
from ete3 import NCBITaxa
from HiTaxon.train_utils import build_kmers


def lineage_creator(genus_names, ncbi):
    """
    Determine the entire lineage of all genera of interest
    Args:
        genus_names: list of all genera
        ncbi = NCBITaxa()
    """
    ranks = ["genus", "family", "class", "order", "phylum"]
    #Get taxids of lineage for each genus
    name2taxid = [ncbi.get_name_translator([taxa])[taxa][0] for taxa in genus_names]
    lineages = [ncbi.get_lineage(taxid) for taxid in name2taxid]
    taxid2rank_lineages = [ncbi.get_rank(lineage) for lineage in lineages]
    rank2taxid_lineages = []
    #Invert dict from taxid:rank to rank:taxid for each genus lineage
    for lineage in taxid2rank_lineages:
        rank2taxid_lineages.append({rank:taxid for taxid,rank in lineage.items()})
    #Convert taxid to name across all ranks for each genus
    lineages_subset = []
    for lineage in rank2taxid_lineages:
        lineage_subset = []
        for rank in ranks:
            if rank in lineage.keys():
                lineage_subset.append(ncbi.get_taxid_translator([lineage[rank]])[lineage[rank]])
            else:
                lineage_subset.append("unclassified")
        lineages_subset.append(lineage_subset) 
    lineages_df = pd.DataFrame(lineages_subset, columns = ranks) 
    return lineages_df

def genus_select(taxa, rank, lineages, mode):
    """
    Determine which sequence data (which is grouped according to genus) to use to train classifiers for a specific taxa at a particular rank
    Args:
        taxa: taxa of interest
        rank: rank of interest
        lineages: dataframe containing entire lineage of all genera of interest
        mode: If mode is "binary", need to select random negative example to create binary classifer for that taxa else use all descendent examples
    """
    genus = []
    negative = []
    #For taxa with multiple descedents, create list of a subset of genera corresponding to these descendents
    if mode == "multi":
        members = lineages[lineages[rank] == taxa]
        genus.append(np.unique(members["genus"].values))
    #For taxa with one descendent, need to select negative example
    else:
        genus.append(lineages[lineages[rank] == taxa]["genus"].values)
        phylum_of_taxa = lineages[lineages[rank] == taxa]["phylum"].values
        remaining_subset = lineages[~(lineages[rank] == taxa)]
        #If all trained species are from the same taxa at a particular rank, need to add atleast one outgroup to create binary classifiers
        if len(remaining_subset) == 0:
            print(f"All organisms are in {taxa}, please incorporate an organism not apart of this taxa")
            return None
        #Preferentially select negative examples from the same phylum
        related = remaining_subset[remaining_subset["phylum"].isin(phylum_of_taxa)]
        if len(related) == 0:
            #If none exist, then use datasets from random genera as the negative example
            negative_example = remaining_subset["genus"]
        else:
            negative_example = related["genus"]
        negative.append(np.unique(negative_example.values))
    return [genus, negative]

def create_lcpn_data(taxa, selected_genera, next_rank, lineages, output_path, model_path):
    """
    Create LCPN training data for specific taxa at a particular rank
    Args:
        taxa: taxa  of interest
        selected_genera: genera to source sequence data for training 
        next_rank: rank corresponding to descendent taxa being predicted
        lineages: dataframe containing entire lineage of all trained species
        output_path: path in which RefSeq sequences are downloaded and processed 
        model_path: path in which models are saved
    """
    #List of all genera corresponding to descendents of the taxa of interest
    genera = selected_genera[0][0]
    total_positive = 0
    #Create training file
    f = open(f"{model_path}/{taxa}_train.txt", "w")
    for genus in genera:
        print(genus)
        #Create label for descendent of taxa
        label = lineages[lineages["genus"] == genus][next_rank].values
        #Add k-merized sequences to training data
        for record in SeqIO.parse(f"{output_path}/{genus}/{genus}_shuffled.fa", "fasta"):
            seq2kmer = build_kmers(str(record.seq), 13)
            f.write("__label__" + str(label[0])  +  " " + seq2kmer  + "\n")
            total_positive +=1
    counter = 0
    #If training binary classifier, add negative genus
    if len(selected_genera[1]) > 0:
        negatives = selected_genera[1][0]
        for negative in negatives:
            #To create a balanced training set for the binary classifier, use an equal amount of positive and negative instances
            if total_positive == counter:
                break       
            print(negative)
            for record in SeqIO.parse(f"{output_path}/{negative[0]}/{negative[0]}_shuffled.fa", "fasta"):
                label = "different"
                seq2kmer = build_kmers(str(record.seq), 13)
                f.write("__label__" + str(label) +  " " + seq2kmer  + "\n")
                counter +=1
                if total_positive == counter:
                    break       
    f.close()
    
def create_lcl_data(genus_names, rank, lineages, output_path, model_path):
    """
    Create LCL training data for specific taxonomic rank
    Args:
        genus_names: list of all genera of interest
        rank: rank of interest
        lineages: dataframe containing entire lineage of all genera
        output_path: path in which RefSeq sequences are downloaded and processed 
        model_path: path in which models are saved
    """
    f = open(f"{model_path}/{rank}_train.txt", "w")
    if rank == "species":
        for genus in taxa_of_interest:
            for record in SeqIO.parse(f"{output_path}/{genus}/{genus}_shuffled.fa", "fasta"):
                #Get species name from header of FASTA files
                label = (record.id).split("^")[1].split("<")[0]

                seq2kmer = build_kmers(str(record.seq), 13)
                f.write("__label__" + str(label)  +  " " + seq2kmer  + "\n")
    else:
        #Create a dictionary mapping genus (k) to higher-order taxa (v) at the rank of interest
        genus_mapped2_rank = dict(zip(lineages["genus"], lineages[rank])) 

        for k,v in genus_mapped2_rank.items():
            for record in SeqIO.parse(f"{output_path}/{k}/{k}_shuffled.fa", "fasta"):
                label = v
                seq2kmer = build_kmers(str(record.seq), 13)
                f.write("__label__" + str(label)  +  " " + seq2kmer  + "\n")
    f.close()


def model_train(to_train, model_path, lr = 0.5, wordNgrams = 1, epoch = 15, dim = 200, thread = 80):
    """
    Train ML model for specific taxa/rank
    Args:
        to_train: Can either be taxa of interest for LCPN models or ranks for LCL models; serves as the pre-fix for training file
        model_path: path in which models are saved
        lr: learning rate
        wordNgrams: max length of word ngram
        epoch: epochs
        dim: Embedding dimension
        thread: Number of threads to use
    """
    model = ft.train_supervised(f"{model_path}/{to_train}_shuffled.txt", lr = lr, wordNgrams = wordNgrams, epoch = epoch,  dim = dim, thread = thread)
    return model

 
