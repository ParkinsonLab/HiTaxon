import argparse
import pandas as pd

def kaiju_kraken_ensemble(report_path, report_name):
    """
    Ensemble Kaiju and Kraken2
    Args:
        report_path: path in which Kraken2 and Kaiju outputs are saved
        report_name: names of Kraken2 and Kaiju outputs
    """
    kaiju_preds = pd.read_csv(f"{report_path}/{report_name}_lineage_kaiju.csv")
    kraken_preds = pd.read_csv(f"{report_path}/{report_name}_lineage_kraken.csv")
    ensemble = []
    ranks = ["species", "genus", "family", "order", "class", "phylum"]
    for i in range(0, len(kraken_preds)):
        #If no outputs are provided at the species level, use Kraken2's predictions
        if str(kraken_preds["species"].iloc[i]) == "nan" and str(kaiju_preds["species"].iloc[i]) == "nan":
            ensemble.append(kraken_preds.iloc[i])
        #else if no outputs are provided by Kraken2 at the species level, use Kaiju's predictions
        elif str(kraken_preds["species"].iloc[i]) == "nan":
            ensemble.append(kaiju_preds.iloc[i])
        #by default, use Kraken2's predictions
        else:
            ensemble.append(kraken_preds.iloc[i])
    ensemble_df = pd.DataFrame(ensemble)
    return ensemble_df[ranks]

"""
Create Kaiju-Kraken2 Ensemble
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("report_path", type = str, help = "path in which Kraken2 and Kaiju outputs are saved")
    parser.add_argument("report_name", type = str, help = "names of Kraken2 and Kaiju outputs")
    args = parser.parse_args()

    report_path = args.report_path
    report_name = args.report_name

    ensemble_output = kaiju_kraken_ensemble(report_path, report_name)
    ensemble_output.to_csv(f"{report_path}/{report_name}_kaiju_kraken.csv")

if __name__ == "__main__":
    main()

