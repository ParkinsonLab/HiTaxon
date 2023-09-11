from Bio import SeqIO
from ete3 import NCBITaxa

def id_generator(species, reference_assembly, ncbi):
    """
    Get species taxid
    Args:
        species: species of interest
        reference_assembly = DataFrame of NCBI RefSeq Summary
        ncbi: NCBITaxa()
    """
    genus = species.split(" ")[0]
    #Translate name to taxid
    name2taxid = ncbi.get_name_translator([species])

    #In instances NCBI ete3 fails, use taxid used in RefSeq Assembly
    if len(name2taxid) == 0:
        name2taxid = reference_assembly[reference_assembly["species"] == species]["species_taxid"].values
        if len(name2taxid) == 0:
            return None
        else:
            name2taxid = name2taxid[0]
    else:
        name2taxid = name2taxid[species][0]
    return name2taxid 



def kraken_FASTA_generator(species, name2taxid, output_path):
    """
    Create FASTA files from non-redundant sequences that are compatible with Kraken2
    Args:
        species: species of interest
        name2taxid: taxid for species of interest
        output_path: path in which sequence data is stored and collected
    """
    genus = species.split("_")[0]

    #Writing sequences with Kraken2-compatible headers
    kraken_fasta = open(f"{output_path}/{genus}/{species}/kraken.fa", "w")
    counter = 0
    for record in SeqIO.parse(f"{output_path}/{genus}/{species}/non_redundant.fa", "fasta"):
        kraken_fasta.write(f">sequence{counter}|kraken:taxid|{name2taxid}" + "\n" + str(record.seq) + "\n")
        counter +=1
    kraken_fasta.close()
