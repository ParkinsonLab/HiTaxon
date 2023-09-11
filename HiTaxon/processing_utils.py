import numpy as np
import pandas as pd
import os

from Bio import SeqIO
from sklearn.cluster import OPTICS


def clean(species, output_path):
    """ 
    Remove plasmid sequences from assemblies 
    Args:
        species: species of interest
        output_path: path in which sequence data is collected and processed
    """
    genus = species.split("_")[0]
    downloaded_assemblies = pd.read_csv(f"{output_path}/{genus}/{species}/download_assemblies.txt", header = None)

    for assembly in downloaded_assemblies[0]:
        if not(os.path.exists(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/cds_from_genomic.fna")):
            continue
        assembly_path = f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/genome_seq.tsv"

        #For assemblies without metadata, write the sequences to a new destination
        if not(os.path.exists(assembly_path)) or not(os.path.getsize(assembly_path) > 0):
            retain(species, output_path)
            continue

        assembly_metadata = pd.read_csv(assembly_path, delimiter = "\t")
        assembly_metadata = assembly_metadata[assembly_metadata["Molecule type"] == "Plasmid"]
        plasmid_accession = pd.unique(assembly_metadata["RefSeq seq accession"])
        if os.path.exists(f"{output_path}/{genus}/{species}/{assembly}.fa") and os.path.getsize(f"{output_path}/{genus}/{species}/{assembly}.fa") > 0:
           continue
        
        #Write all non-plasmid sequences from assembly to new FASTA file 
        fasta = open(f"{output_path}/{genus}/{species}/{assembly}.fa", "w")
        for record in SeqIO.parse(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/cds_from_genomic.fna", "fasta"):
            for accession in plasmid_accession:
                if accession in str(record.id):
                    continue
            else:
                fasta.write(">" + record.id + "\n" + str(record.seq) + "\n")
        fasta.close()


def retain(species, output_path):
    """
    Copy all sequences from assembly to new FASTA; used for assemblies without Plasmid Metadata
    Args:
        species: species of interest
        output_path: path in which sequence data is collected and processed
    """
    genus = species.split("_")[0]
    downloaded_assemblies = pd.read_csv(f"{output_path}/{genus}/{species}/download_assemblies.txt", header = None)
    for assembly in downloaded_assemblies[0]:
        if (os.path.exists(f"{output_path}/{genus}/{species}/{assembly}.fa") and os.path.getsize(f"{output_path}/{genus}/{species}/{assembly}.fa") > 0) or not(os.path.exists(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/cds_from_genomic.fna")):
           continue
        fasta = open(f"{output_path}/{genus}/{species}/{assembly}.fa", "w")
        for record in SeqIO.parse(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/cds_from_genomic.fna", "fasta"):
            fasta.write(">" + record.id + "\n" + str(record.seq) + "\n")
        fasta.close()




def query(species, output_path):
    """
    Create ANI query files
        Args:
        species: species of interest
        output_path: path in which sequence data is collected and processed
    """
    genus = species.split("_")[0]
    downloaded_assemblies = pd.read_csv(f"{output_path}/{genus}/{species}/download_assemblies.txt", header = None)
    query = open(f"{output_path}/{genus}/{species}/query.txt", "w")
    f = open(f"{output_path}/missing_assemblies.txt", "a")
    #Write list of file paths for all cleaned/retained FASTA files
    for assembly in downloaded_assemblies[0]:
        #Ignore assemblies that were not downloaded or have no sequence data
        if not(os.path.exists(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/cds_from_genomic.fna")) or os.stat(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/{assembly}/cds_from_genomic.fna").st_size == 0:
            f.write(str(genus) + "," + str(species) + "," + str(assembly) + "," + "\n") 
            continue
        query.write(f"{output_path}/{genus}/{species}/{assembly}.fa" + "\n")
    query.close()
    f.close()




def cluster_default(species, output_path, large = 400):
    """
    Cluster assemblies using ANI
    Args:
        species: species of interest
        output_path: path in which sequence data is collected and processed
        large: Maximum number of assemblies in a query set that can be clustered using the default approach
    """

    genus = species.split("_")[0]
    query = pd.read_csv(f"{output_path}/{genus}/{species}/query.txt", header = None, delimiter = "\t")

    #Cluster with altenative memory efficient approach if too many assemblies, otherwise use default
    if len(query) >= large:
        clustering,assemblies = cluster_alternative(species, output_path)
    else:
        ani_output = pd.read_csv(f"{output_path}/{genus}/{species}/ani_output.matrix", header = None, delimiter = "\t")
        assemblies = list(pd.unique(query[0]))

        #Create assembly x assembly matrix of ANI values
        mtx = np.zeros((len(assemblies), len(assemblies)))
        mtx = pd.DataFrame(mtx, columns = assemblies, index = assemblies)
        
        #If one assembly for a species, no need to create matrix and cluster downstream
        if len(mtx) == 1:
            return None, assemblies[0]

        #If ANI > 80% between assembly i and assembly j use fastANI value, else use 0 (as no output below 80%)
        for assembly_i in assemblies:
            if assembly_i not in list(ani_output[0].values):
                continue
            for assembly_j in assemblies:
                assembly_i_ani = ani_output[ani_output[0] == assembly_i]
                if len(assembly_i_ani[assembly_i_ani[1] == assembly_j][2]) != 0:
                    mtx.loc[assembly_i][assembly_j] = assembly_i_ani[assembly_i_ani[1] == assembly_j][2].values[0]
                else:
                    mtx.loc[assembly_i][assembly_j] = 0
        mtx = np.array(mtx)

        #Generate clusters
        clustering = OPTICS(min_samples=2).fit(mtx)
    return clustering,assemblies




def cluster_alternative(species, output_path):
    """
    Cluster assemblies corresponding to species with > 400 assemblies using an alternative FastANI output to avoid OOM errors 
    Args:
        species: species of interest
        output_path: path in which sequence data is collected and processed
    Note: Alternative matrix output computes ANI slightly different then conventional output: https://github.com/ParBLiSS/FastANI/issues/36
    """
    genus = species.split("_")[0]
    ani_output = open(f"{output_path}/{genus}/{species}/ani_output.matrix").read().splitlines()
    mtx_raw = [row.split("\t") for row in ani_output]
    mtx_raw = mtx_raw[1:]
    mtx_processed = []

    #Create list of query comparisons
    for assemblies in mtx_raw:
        add = []
        for items in assemblies:
            #Impute value of 0 in instances no ANI is provided
             if str(items) == "NA":
                     add.append(0)
             else:
                     add.append(items)
        mtx_processed.append(add)

    #Create diagonal matrix
    diag_mtx = []
    for assembly in mtx_processed:
        total_assemblies = len(assembly) - 1
        #Add trailing 0's to keep length of lists for each assembly consistent
        add_placeholder_zero = len(mtx_processed) - total_assemblies
        ani_values = assembly[1:].copy()
        ani_values.extend([0]*add_placeholder_zero)
        diag_mtx.append([float(ani) for ani in ani_values])
    mtx_unfilled = np.zeros(shape = (len(mtx_processed), len(mtx_processed)))
    for i in range(0, len(mtx_processed)):
        mtx_unfilled[i] = diag_mtx[i]

    #Create assembly x assembly matrix
    mtx = mtx_unfilled + mtx_unfilled.T
    diag_100 = np.repeat(100, len(mtx_processed))
    mtx = mtx + np.diag(diag_100)

    #Generate clusters
    clustering = OPTICS(min_samples=2).fit(mtx)
    assemblies = [i[0] for i in mtx_processed]
    return clustering,assemblies




def find_representative(species, clustering, assemblies, output_path):
    """ 
    Select represenative assemblies from clustering output
    Args:
        species: species of interest
        clustering: OPTICS clustering output
        assemblies: list of assemblies compared corresponding to the species of interest
        output_path: path in which sequence data is collected and processed

    """
    genus = species.split("_")[0]
    rep_assemblies = []
    clustering_results = pd.DataFrame({'cluster': clustering.labels_,'assemblies': assemblies})
    #For each unique cluster, Select assembly with highest N50 value 
    for cluster in np.unique(clustering.labels_):
         cluster_subset = clustering_results[clustering_results["cluster"] == cluster]
         list_of_assemblies = cluster_subset["assemblies"].values
         best_n50 = 0
         for assembly_path in list_of_assemblies:
            if cluster == -1:
                #keep all assemblies that are not a part of a cluster (denoted as -1)
                rep_assemblies.append(assembly_path)
            else:
                assembly = assembly_path.split("/")[-1][:-3]
                if best_n50 == 0:
                   #placeholder in case no N50 metadata is available
                   best_assembly = assembly_path
                assembly_metadata = pd.read_csv(f"{output_path}/{genus}/{species}/{assembly}/ncbi_dataset/data/genome.tsv", delimiter = "\t")
                if assembly_metadata["Assembly Stats Scaffold N50"].values[0] > best_n50:
                   best_n50 = assembly_metadata["Assembly Stats Scaffold N50"].values[0]
                   best_assembly = assembly_path
         #keep assembly with highest N50 score in cluster, excluding assemblies not in cluster 
         if cluster != -1:
               rep_assemblies.append(best_assembly)
    return rep_assemblies


def aggregator(species, output_path):
    """
    Combine coding sequences from different assemblies into single FASTA
    Args:
        species: species of interest
        output_path: path in which sequence data is collected and processed
    """
    genus = species.split("_")[0]
    rep_assemblies = pd.read_csv(f"{output_path}/{genus}/{species}/rep_assemblies.txt", header = None)

    #Write sequences from represenative assemblues in agg_CDS fasta
    f = open(f"{output_path}/{genus}/{species}/agg_CDS.fa", "w")
    for assemblies in rep_assemblies[0]:
        for record in SeqIO.parse(assemblies, "fasta"):
            f.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
    f.close()
