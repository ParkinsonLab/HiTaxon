import argparse
import pandas as pd

from ete3 import NCBITaxa
from HiTaxon.evaluation_utils import expand_lineage, expand_predictions

"""
Reformat Kraken2's prediction to include entire lineage
"""

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_summary", type = str, help = "file path to assembly summary")
    parser.add_argument("report_name", type = str, help = "file name of output")
    parser.add_argument("report_path", type = str, help = "path to store classifier output")
    args = parser.parse_args()

    ncbi = NCBITaxa()
    assembly_summary = args.assembly_summary
    report_name = args.report_name
    report_path = args.report_path
    predictions = report_path + "/" + report_name + ".kraken"
    
    refseq = pd.read_csv(assembly_summary, delimiter = "\t", skiprows = 1)
    refseq = refseq[["#assembly_accession", "refseq_category", "species_taxid", "organism_name", "infraspecific_name", "assembly_level", "genome_rep"]]
    refseq["species"] = refseq["organism_name"].apply(lambda x: " ".join(x.replace("_"," ").split(" ")[:2]))

    #Convert taxids to names, alongside expanding lineage
    kraken_predictions = expand_predictions(predictions, ncbi, refseq)
    kraken_predictions.to_csv(f"{report_path}/{report_name}_lineage_kraken.csv")



if __name__ == "__main__":
    main()
