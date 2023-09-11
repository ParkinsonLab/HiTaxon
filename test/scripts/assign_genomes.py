import argparse
import numpy as np
import os
import pandas as pd
import shutil


def reorganize_genomes(test_genomes, species, test_path):
    """
    Move test assemblies to their species folder, and determine which assemblies were actually downloaded

    Args:
        test_genomes: file path for test_genomes.csv
        species: list of species interested in aquiring assemblies for
        test_path: path of test directory
    """

    test_genomes = pd.read_csv(test_genomes)
    accessions = test_genomes["0"].values
    downloaded_genomes = []
    #Move downloaded assemblies from assembly_accession/ncbi_dataset/data/assembly_accession to assembly_accession/
    for accession in accessions:
        accession = accession.split("\n")[0]
        path_of_genomes = []
        for root, dirs, files in os.walk(f"{test_path}/{accession}"):
            for file in files:
                    if file.endswith(".fna"):
                             path_of_genomes.append(os.path.join(root, file))
                             name = test_genomes["1"][test_genomes["0"] == accession]
        for genome in path_of_genomes:
            new_path = test_path + "/" + str(accession) + "/" + str(accession) + ".fna"
            shutil.move(genome,new_path)

    #Create directories for each species
    for taxa in species:
        os.makedirs(f"{test_path}/{taxa.replace(' ','_')}", exist_ok=True)

    #Move downloaded assemblies to their corresponding species directory
    for index, row in test_genomes.iterrows():
        name = row["1"].replace(" ","_")
        accession = row["0"]
        old = f"{test_path}/{accession}/{accession}.fna"
        new = f"{test_path}/{name}/{accession}.fna"
        if os.path.exists(new):
                downloaded_genomes.append((accession, name))
                continue
        elif os.path.exists(old):
                downloaded_genomes.append((accession, name))
                shutil.move(old, new)
    return downloaded_genomes



def mock_abundance(downloaded_genomes, refseq_assembly, state, test_path):
     """
     For each species, assign an assembly alongside a relative abundance

     Args:
         downloaded_genomes: list of tuples pertaining to GenBank assemblies which were downloaded and their corresponding species
         refseq_assembly: file path for NCBI RefSeq Assembly
         state: Random seed used in this simulation
         test_path: path of test directory
     """

    downloaded_genomes = pd.DataFrame(downloaded_genomes)
    #A secondary filter used to remove RefSeq assemblies; extra precaution to prevent data leakage
    refseq = pd.read_csv(refseq_assembly, delimiter = "\t", skiprows = 1)
    refseq = refseq[["# assembly_accession", "refseq_category", "species_taxid", "organism_name", "infraspecific_name", "assembly_level", "genome_rep"]]
    refseq["# assembly_accession"] = refseq["# assembly_accession"].apply(lambda x: x.split(".")[0])
    refseq["# assembly_accession"] = refseq["# assembly_accession"].apply(lambda x: x.replace("GCF","GCA"))
    accessions_no_versions = downloaded_genomes[0].apply(lambda x: x.split(".")[0])
    refseq_accessions = pd.unique(refseq["# assembly_accession"])
    downloaded_genomes = downloaded_genomes[~accessions_no_versions.isin(refseq_accessions)]

    #Sample one assembly per species
    downloaded_sample = downloaded_genomes.groupby(1, group_keys=False).apply(lambda df: df.sample(1, random_state = state))

    #Assign each species a relative abundance
    mu, sigma = 0, 1.5
    np.random.seed(state)
    dist = np.random.lognormal(mu, sigma, len(downloaded_sample))
    log_dist = (dist  - dist.min()) / (dist .max() - dist .min())
    log_dist = log_dist / log_dist.sum()
    np.random.shuffle(log_dist)
    downloaded_sample["dist"] = log_dist

    #Count the number of coding sequences for each species
    transcript_counter = []
    downloaded_sample[1] = downloaded_sample[1].apply(lambda species: species.replace(" ","_"))
    for index, row in downloaded_sample .iterrows():
            genome = row[0] + ".fna"
            name = row[1]
            counter = 0
            for record in SeqIO.parse(f"{test_path}/{name}/{genome}", "fasta"):
                    counter +=1
            transcript_counter.append(counter)
    downloaded_sample["num_transcript"] =  transcript_counter

    #With a baseline of 1000000, assign each species a baseline number of reads per CDS (minimum 1 read)
    test_reads = 1000000
    baselines = []
    for index, row in downloaded_sample.iterrows():
            total = test_reads * row["dist"]
            baseline = int(total/row["num_transcript"])
            if baseline < 1:
                    baseline = 1
            baselines.append(baseline)
    downloaded_sample["baselines"] =  baselines
    return downloaded_sample


"""
Define simulated environments 
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("refseq_assembly", type = str, help = "file path for NCBI RefSeq Assembly")
    parser.add_argument("species", type = str, help = "file path for .txt file with the species of interest to be tested") 
    parser.add_argument("test_path", type = str, help = "path of test directory")
    parser.add_argument("state", type = str, help = "seed used for this simulation")
    args = parser.parse_args()

    refseq_assembly = args.refseq_assembly
    species = open(args.species).read().splitlines()
    test_path = args.test_path
    test_state = int(args.state)

    #Determine which assemblies were downloaded and move them to their correct species directory
    downloaded_genomes = reorganize_genomes(f"{test_path}/test_genomes.csv", species, test_path)
    
    #Define 10 simulated environments
    for state in range(1, 11):
        test_community = mock_abundance(downloaded_genomes, refseq_assembly, test_state, test_path)
        test_community.to_csv(f"{test_path}/sample_{state}.csv")


if __name__ == "__main__":
    main()
