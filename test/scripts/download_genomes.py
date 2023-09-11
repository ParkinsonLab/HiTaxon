import argparse
import pandas as pd

def genbank_selection(refseq_summary, genbank_summary, species, random_seed):
    """
    Create a list of genomes from Genbank that were unseen (I.E RefSeq-Excluded) during data collection & processing + training

    Args:
        refseq_summary: file path for NCBI RefSeq Assembly
        genbank_summary: file path for NCBI Genbank Assembly
        species: list of species interested in aquiring assemblies for
        random_seed: Random seed used to shuffle GenBank assembly for this simulation
    """
    #Remove all RefSeq assemblies from Genbank Assembly Summary
    refseq = pd.read_csv(refseq_summary, delimiter = "\t", skiprows = 1)
    genbank = pd.read_csv(genbank_summary, delimiter = "\t", skiprows = 1)
    refseq["accession_reformatted"] = refseq["# assembly_accession"].apply(lambda x: x.split("_")[1].split(".")[0])
    genbank["accession_reformatted"] = genbank["# assembly_accession"].apply(lambda x: x.split("_")[1].split(".")[0])
    shared = list(set(genbank["accession_reformatted"].values).intersection(set(refseq["accession_reformatted"].values)))
    genbank = genbank[~(genbank["accession_reformatted"].isin(shared))]
    genbank[~(genbank["refseq_category"] == "representative genome")]

    #Subset Genbank dataframe for list of species of interest
    genbank["organism_name"] = genbank["organism_name"].apply(lambda x: " ".join(x.split(" ")[:2]))
    genbank_subset = genbank[genbank["organism_name"].isin(species)]

    #Shuffle Genbank dataframe and select a maximum of 10 assemblies per species to use for training
    genbank_subset = genbank_subset.sample(frac=1, random_state = random_seed)
    genome_counter = {species_name:0 for species_name in pd.unique(genbank_subset["organism_name"])}
    genome_list = []
    for index, row in genbank_subset.iterrows():
        if genome_counter[row["organism_name"]] > 10:
            continue
        else:
             genome_list.append((row["# assembly_accession"], row["organism_name"]))
             genome_counter[row["organism_name"]] +=1
    return genome_list



"""
Test Assembly Selection
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("refseq_summary", type = str, help = "file path for NCBI RefSeq Assembly")
    parser.add_argument("genbank_summary", type = str, help = "file path for NCBI Genbank Assembly")
    parser.add_argument("species", type = str, help = "file path for .txt file with the species of interest to be tested") 
    parser.add_argument("test_path", type = str, help = "path of test directory")
    parser.add_argument("random_seed", type = str, help = "seed used for this simulation")
    args = parser.parse_args()

    refseq_summary = args.refseq_summary
    genbank_summary = args.genbank_summary
    species = open(args.species).read().splitlines()
    test_path = args.test_path
    test_seed = int(args.random_seed)

    #Select Genbank assemblies and save selection to test path
    test_genomes = genbank_selection(refseq_summary, genbank_summary, species, test_seed)
    test_genomes = pd.DataFrame(test_genomes)
    test_genomes.to_csv(f"{test_path}/test_genomes.csv")

if __name__ == "__main__":
    main()
