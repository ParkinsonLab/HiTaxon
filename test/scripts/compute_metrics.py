import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef

""" Compute precision, recall, F1 Score, MCC across all taxonomic levels """
def compute_metrics(truth, predictions, rank):
    tested_taxa = np.unique(truth[rank].values)
    per_recall = {i:0 for i in tested_taxa}
    per_precision = {i:0 for i in tested_taxa}
    per_f1 = {i:0 for i in tested_taxa}
    for taxa in tested_taxa:
        subset_truth = truth[truth[rank] == taxa].index.values
        subset_pred = predictions[predictions[rank] == taxa].index.values
        tp = np.sum(predictions[rank].iloc[subset_truth] == taxa)
        fp = len(subset_pred) - tp
        if tp == 0:
            per_precision[taxa] = 0
        else:
            per_precision[taxa] = tp/(tp + fp)
        per_recall[taxa] = tp/(len(subset_truth))
        if per_recall[taxa] == 0 and per_precision[taxa] == 0:
            per_f1[taxa] = 0
        else:
            per_f1[taxa] = (2*per_recall[taxa]*per_precision[taxa])/(per_recall[taxa]+per_precision[taxa])
    f1 = np.mean([v for v in per_f1.values()])
    avg_p = np.mean([v for v in per_precision.values()])
    avg_r = np.mean([v for v in per_recall.values()])
    mcc = matthews_corrcoef(truth[rank].values, predictions[rank].values)
    return [mcc, f1, avg_p, avg_r]


""" 
Compute metrics for classifier relative to ground truth
"""
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("classifier_output", type = str, help = "file path classifier output")
    parser.add_argument("ground_truth", type = str, help = "file path to ground truth for test set")
    parser.add_argument("report_path", type = str, help = "path in which outputs are saved")
    parser.add_argument("report_name", type = str, help = "file name of output")
    args = parser.parse_args()


    truth = pd.read_csv(args.ground_truth)
    classifier = pd.read_csv(args.classifier_output)
    report_path = args.report_path
    report_name = args.report_name
    ranks = ["species", "genus", "family", "order", "class", "phylum"]
    for rank in ranks:
        truth[rank] = truth[rank].apply(lambda x: str(x).replace("_", " "))
        if rank == "species":
            classifier[rank] = classifier[rank].apply(lambda x: " ".join(str(x).replace("_", " ").split(" ")[:2]))
        else:
            classifier[rank] = classifier[rank].apply(lambda x: " ".join(str(x).replace("_", " ").split(" ")[:1]))
        results = compute_metrics(truth, classifier, rank)
    results_df = pd.DataFrame(results, columns = ["f1", "precision", "recall", "mcc"])
    results_df.to_csv(f"{report_path}/{report_name}_metrics.csv")



if __name__ == "__main__":
    main()