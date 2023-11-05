from Bio import SeqIO

def bwa_FASTA_generator(species, output_path, bwa_path):
    """
    Create FASTA files from non-redundant sequences that are compatible for building BWA genus to set of species pair indices
    Args:
        species: species of interest
        output_path: path in which sequence data is stored and collected
        bwa_path: path in which bwa fasta are stored
    """
    genus = species.split("_")[0]

    #Writing sequences with Kraken2-compatible headers
    bwa_fasta = open(f"{bwa_path}/{genus}.fa", "a")
    counter = 0
    for record in SeqIO.parse(f"{output_path}/{genus}/{species}/non_redundant.fa", "fasta"):
        bwa_fasta.write(f">sequence{counter}|{species}" + "\n" + str(record.seq) + "\n")
        counter +=1
    bwa_fasta.close()
