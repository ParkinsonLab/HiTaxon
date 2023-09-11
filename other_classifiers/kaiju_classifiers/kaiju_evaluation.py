import argparse
import numpy as np
import pandas as pd
import re

from ete3 import NCBITaxa


def kaiju_expand_lineage(prediction, ncbi):
    """
    Determine the taxonomic lineage of a single prediction (phylum to lowest possible rank, with species as the maximum)
    Args:
        prediction: predicted taxa
        ncbi: NCBITaxa()
    """
    #Get entire lineage of taxa, in taxids
    lineage = ncbi.get_lineage(prediction)

    #Convert taxids to name for relevant ranks
    named_lineage = {"species": "NA", "genus": "NA", "family": "NA", "order": "NA", "class": "NA", "phylum": "NA"}
    if lineage is None:
        return named_lineage
    for taxa in lineage:
        rank = ncbi.get_rank([taxa])[taxa]
        if rank in named_lineage.keys():
            named_lineage[rank] =  ncbi.get_taxid_translator([taxa])[taxa]
    return named_lineage

def kaiju_expand_predictions(predictions, ncbi):
    """
    Determine the taxonomic lineage of all predictions (phylum to lowest possible rank, with species as the maximum)
    Args:
        predictions: path to kaiju taxonomic predictions
        ncbi: NCBITaxa()
    """
    kaiju_predictions = open(predictions, "r")
    kaiju_predictions_formatted = []
    kaiju_predictions_position = []
    for prediction in kaiju_predictions:
        prediction = prediction.split("\t")
        #Get the taxid of prediction
        kaiju_predictions_formatted.append(prediction[2])
        #Get the positional counter on header of Kaiju test FASTA
        r = re.compile(r'(\D+)(\d+)')
        x = int(r.match(prediction[1]).groups()[1])
        kaiju_predictions_position.append(x)
    unique_predictions = np.unique(np.array(kaiju_predictions_formatted))
    lineages = {}
    #Get lineage for each unique predictions
    for prediction in unique_predictions:
        lineages[prediction] = kaiju_expand_lineage(int(prediction), ncbi)
    #Get lineage for each prediction
    kaiju_expanded = [lineages[prediction] for prediction in kaiju_predictions_formatted]
    kaiju_expanded = pd.DataFrame(kaiju_expanded)
    #Sort predictions in appropriate order using position counter
    kaiju_expanded["proper_order"] = kaiju_predictions_position
    kaiju_expanded = kaiju_expanded.sort_values(by=["proper_order"])
    return kaiju_expanded

"""
Expand lineage of Kaiju's predictions
"""

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("report_path", type = str, help = "path in which Kaiju output is saved")
    parser.add_argument("report_name", type = str, help = "name of Kaiju output")
    args = parser.parse_args()

    ncbi = NCBITaxa()
    report_name = args.report_name
    report_path = args.report_path

    predictions = report_path + "/" + report_name + "_kaiju_taxa.out"
    kaiju_predictions = kaiju_expand_predictions(predictions, ncbi)
    kaiju_predictions.to_csv(f"{report_path}/{report_name}_lineage_kaiju.csv")



if __name__ == "__main__":
    main()

