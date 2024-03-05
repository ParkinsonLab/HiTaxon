from Bio import SeqIO
import fasttext as ft
import numpy as np
import pandas as pd

from ete3 import NCBITaxa
from tqdm import tqdm

def build_kmers(sequence, ksize):
    """ 
    Convert single sequence into kmers 
    Args:
        Sequence: sequence/read in fasta file
        ksize: k-mer size to generate
    """
    kmers = []
    ksize = int(ksize)
    n_kmers = len(sequence) - ksize + 1
    for num in range(n_kmers):
        kmer = sequence[num:num + ksize]
        kmers.append(kmer)
    kmers = " ".join(kmers)
    return kmers



def data_creation(genus, ksize, seen_genera, ncbi, output_path, model_path):
    """
    Create k-merized training set from simulated reads 
    Args:
        genus: genus of interest
        ksize: size of k-mers to be created from sequence data
        seen_genera: all the genera in which training reads have been simulated for; necessary for creating negative examples
        ncbi: NCBITaxa()
        output_path: path in which sequences are collected and processed
        model_path: path in which models are stored

    """
    #Keep track of numnber of k-merized reads for each species
    instances = 0
    species_and_sequences = {}

    #K-merize sequences until 500000 reads are k-merized for each species
    f = open(f"{model_path}/{genus}_train.txt","w")
    for record in tqdm(SeqIO.parse(f"{output_path}/{genus}/{genus}_shuffled.fa", "fasta")):
        label = (record.id).split("^")[1].split("<")[0]
        if label not in species_and_sequences.keys():
            species_and_sequences[label] = 0
        if species_and_sequences[label] < 500000:
            species_and_sequences[label] +=1
            kmer = build_kmers(str(record.seq), ksize)
            f.write("__label__" + label + " " + kmer  + "\n")
        instances +=1

    #In instances where genus leads to single species, need negative example to create binary classifier
    if len(species_and_sequences.keys()) == 1:
        negative_training_set = negative_set(genus, seen_genera, ncbi)
        species_and_sequences["different"] = []
        #For the negative set of examples, 500000 reads are k-merized in total (i.e 5 negative species, the 100K reads used from each species)
        negative_ancestors = [negative_example[-1][0] for negative_example in negative_training_set]
        instances2add = int(500000/len(negative_ancestors))
        for negative_genus in negative_ancestors:
            instances = 0
            for record in tqdm(SeqIO.parse(f"{output_path}/{negative_genus}/{negative_genus}_shuffled.fa", "fasta")):
                label = "different"
                kmer = build_kmers(str(record.seq), ksize)
                f.write("__label__" + label + " " + kmer  + "\n")
                instances +=1
                if instances == instances2add:
                    break
    f.close()
    return None


def negative_set(binary_taxa, seen_genera, ncbi):
    """ 
    Determine negative examples for single species binary classifiers 
    Args:
        binary_taxa: genus-level taxa with only one species
        seen_genera: list of all genera in which reads have been simulated for
        ncbi = NCBITaxa()

    """
    #Find taxonomic lineage of each genus in seen_genera
    genera_taxa_id = []
    for genus in seen_genera:
        genera_taxa_id.append(ncbi.get_name_translator([genus])[genus][0])
    genera_lineages = [ncbi.get_lineage(genus_id) for genus_id in genera_taxa_id]
    genera_lineages = [ncbi.get_rank(lineage) for lineage in genera_lineages]
    rank2taxid_lineages = []
    for lineage in genera_lineages:
            rank2taxid_lineages.append({rank:taxa_id for taxa_id,rank in lineage.items()})
    fam2phylum_lineages = []
    ranks = ["family", "order", "class", "phylum"]
    for lineage in rank2taxid_lineages:
        single_fam2phylum = []
        for rank in ranks:
            if rank in lineage.keys():
                single_fam2phylum.append(lineage[rank])
            else:
                single_fam2phylum.append("none")
        fam2phylum_lineages.append(single_fam2phylum)
    fam2phylum_lineages_df = pd.DataFrame(fam2phylum_lineages, columns = ranks)
    gen2phylum_lineages_df =  pd.concat([pd.DataFrame(seen_genera), fam2phylum_lineages_df], axis = 1)

    #Keep lineages of all genera besides the binary genera
    remaining_lineages_df = gen2phylum_lineages_df[~(gen2phylum_lineages_df[0] == binary_taxa)]
    negative_set_ancestor_data = []

    #Preferentially use more related organisms as negative examples, terminating at the most descriptive rank in which atleast one negative example is found
    for rank in ranks:
        ancestors_from_diff_genera = remaining_lineages_df[remaining_lineages_df[rank].isin([gen2phylum_lineages_df[gen2phylum_lineages_df[0] == binary_taxa][rank].values[0]])]
        if len(ancestors_from_diff_genera) == 0 and rank != "phylum":
            continue
        elif len(ancestors_from_diff_genera) >=1:
            negative_set_ancestor_data.append((binary_taxa, rank, list(ancestors_from_diff_genera[0].values)))
            break
        else:
            #If not organisms are related, use a random species from a random genera as the negative example
            random_genera_idx = np.random.randint(len(remaining_lineages_df))
            random_genera = list(remaining_lineages_df[0].values)[random_genera_idx]
            negative_set_ancestor_data.append((binary_taxa, rank, [random_genera]))
    return negative_set_ancestor_data



def model_train(genus, model_path, lr = 0.5, wordNgrams = 1, epochs = 15, dim = 200, thread = 80):
    """ 
    Train FastText Model
    Args:
        genus: genus of interest
        model_path: path in which models are stored
        lr: learning rate
        wordNgrams: max length of word ngram
        epochs: number of epochs to train
        dim: embedding dimensions
        thread: CPU threads to use for training
    """
    model = ft.train_supervised(f"{model_path}/{genus}_train_shuffled.txt", lr = lr, wordNgrams = wordNgrams, epoch = epochs,  dim = dim, thread = thread)
    return model
