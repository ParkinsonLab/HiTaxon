import pandas as pd

def assembly_collection(assembly_summary, genus_names, output_path):
    """
    Filters RefSeq for high-quality formal-named species assemblies encompassed in the genera of interest

    Args:
        assembly_summary: NCBI RefSeq assembly summary
        genus_names: list of genera to collect assemblies for
        output_path: path in which sequence data is collected and processed
    """
    #Subset RefSeq assembly for formal-named species in specified genera
    refseq = pd.read_csv(assembly_summary, delimiter = "\t", skiprows = 1)
    refseq = refseq[["#assembly_accession", "refseq_category", "species_taxid", "organism_name", "infraspecific_name", "assembly_level", "genome_rep"]]
    refseq["species"] = refseq["organism_name"].apply(lambda x: " ".join(x.split(" ")[:2]))
    refseq["genus"] = refseq.organism_name.apply(lambda x: x.split(" ")[0])
    refseq_of_interest = refseq[refseq.genus.isin(genus_names)]
    formal_named_species = []
    for species in pd.unique(refseq_of_interest["species"]):
        if "sp." in species:
                continue
        elif "phage" in species:
                continue
        else:
             	formal_named_species.append(species)
    refseq_of_interest = refseq_of_interest[refseq_of_interest["species"].isin(formal_named_species)]
    
    #Create list of genera not found in RefSeq
    missing_genera = set(genus_names) - set(refseq_of_interest["genus"].values)

    #For each species, progressively decrease quality threshold until have found atleast 10 assemblies
    maximum = 10
    refseq_of_interest = refseq_of_interest[refseq_of_interest["genome_rep"] == "Full"]
    species_files = {}
    for species in formal_named_species:
        assemblies = []
        species_subset = refseq_of_interest[refseq_of_interest['species'] == species]
        species_name = species.replace(" ","_")
        genus = species.split(" ")[0]
        assembly_level = {"Complete Genome":[], "Chromosome":[], "Scaffold":[], "Contig":[]}
        assembly_quality = ["reference genome", "representative genome", "na"]
        keep_assemblies = []
        strains_encountered = []
        maximum_met = False
        for level in assembly_level.keys():
            if maximum_met == True:
                break
            level_subset = species_subset[species_subset["assembly_level"] == level]
            for quality in assembly_quality:
                added = level_subset[level_subset["refseq_category"] == quality]["#assembly_accession"].values
                if len(keep_assemblies) + len(added) >= maximum:
                    assembly_level[level].extend(added)
                    maximum_met = True
                    break
                else:
                    keep_assemblies.extend(added)
                    strains_encountered.extend(level_subset[level_subset["refseq_category"] == quality]["infraspecific_name"].values)

        #When > 10 assemblies found, only add additional assemblies at the lowest quality threshold which correspond to unseen strains.
        strains_encountered_subset = species_subset[species_subset["infraspecific_name"].isin(set(strains_encountered))]
        strains_encountered_subset = strains_encountered_subset["#assembly_accession"].values
        for level in assembly_level.keys():
            if len(assembly_level[level]) > 0:
                for assembly in assembly_level[level]:
                    if assembly not in strains_encountered_subset or maximum_met == False:
                        keep_assemblies.append(assembly)
        species_files[species] = keep_assemblies
    return species_files, missing_genera

