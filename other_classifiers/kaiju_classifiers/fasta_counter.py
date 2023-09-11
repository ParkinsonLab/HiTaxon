import argparse
from Bio import SeqIO

"""
Add positional counter to FASTA header to be analyzed; allows for easier analysis of Kaiju results downstream
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence_file", type = str, help = "File path for Fasta to be analyzed")
    args = parser.parse_args()
    sequence_file = args.sequence_file
    #Loop through reads in FASTA file and append counter value to header
    counter = 0
    fasta_counted = open(sequence_file.split(".")[0] + "_sorted.fa", "w")
    for record in SeqIO.parse(sequence_file, "fasta"):
        fasta_counted.write(">" + record.id + str(counter) + "\n" + str(record.seq) + "\n")
        counter +=1
    fasta_counted.close()


if __name__ == "__main__":
     main()
