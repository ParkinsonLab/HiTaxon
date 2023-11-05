import numpy as np
import pandas as pd
import fasttext as ft
import os

from Bio import SeqIO
from ete3 import NCBITaxa
from HiTaxon.train_utils import build_kmers


def expand_lineage(prediction, ncbi, reference_assembly):
    """
    Determine taxonomic lineage of prediction (Phylum to prediciton)
    Args:
        prediction: predicted taxa
        ncbi: NCBITaxa()
    """
    #Get taxid at each rank for prediction
    named_lineage = {"species": "NA", "genus": "NA", "family": "NA", "order": "NA", "class": "NA", "phylum": "NA"}
    try:
        lineage = ncbi.get_lineage(prediction) 
    except:
        taxid2name = reference_assembly[reference_assembly["species_taxid"] == prediction]["species"].values[0]
        named_lineage = {"species": taxid2name, "genus": "other", "family": "other", "order": "other", "class": "other", "phylum": "other"}
        return named_lineage
    if lineage is None:
        return named_lineage
    #Convert taxids to names 
    for taxa in lineage:
        rank = ncbi.get_rank([taxa])[taxa]
        if rank in named_lineage.keys():
            named_lineage[rank] =  ncbi.get_taxid_translator([taxa])[taxa]
    return named_lineage

def expand_predictions(predictions, ncbi, reference_assembly):
    """
    Return taxonomic lineage of all Kraken2 predictions 
    Args:
        predictions: file path of Kraken2 output
        ncbi: NCBITaxa()
    """
    kraken_predictions = pd.read_csv(predictions, delimiter  = "\t", header = None)
    kraken_predictions = kraken_predictions[2].values
    unique_predictions = np.unique(np.array(kraken_predictions))
    lineages = {}
    #Get lineage for each unique prediction
    for prediction in unique_predictions:
        lineages[prediction] = expand_lineage(prediction, ncbi, reference_assembly)
    #Expand lineage of all predictions
    kraken2_expanded = [lineages[prediction] for prediction in kraken_predictions]
    kraken2_expanded = pd.DataFrame(kraken2_expanded)
    return kraken2_expanded

def fasta2kmer(fasta, report_path, report_name, ksize = 13):
    """
    K-merize FASTA file 
    Args:
        fasta: file path to fasta
        report_path: path to store classifer output
        report_name: file name of output
        ksize: size of k-mers to be created from sequence data
    """
    converted = []
    for record in SeqIO.parse(f"{fasta}", "fasta"):
            kmer = build_kmers(str(record.seq), ksize)
            converted.append(kmer + "\n")
    f = open(f"{report_path}/{report_name}_kmer.txt", "w")
    for kmer in converted:
        f.write(kmer)
    f.close()

def fasta2bwa(fasta, report_path, report_name):
   """
   Format FASTA file to allow faster taxonomic classification using BWA and Kraken2
   Args:
        fasta: file path to fasta
        report_path: path to store classifer output
        report_name: file name of output
   """
   f = open(f"{report_path}/{report_name}_bwa.fa", "w")
   counter = 0
   for record in SeqIO.parse(f"{fasta}", "fasta"):
        f.write(">" + record.id + "^" +  str(counter) + "\n" + str(record.seq) + "\n")
        counter +=1
   f.close()

def evaluation(report_path, report_name, model_path, min_threshold = 0.5):
    """
    Given Kraken2's genus classifications, generate species-level predictions using machine learning classifiers
    Args:
        report_path: path to store classifer output
        report_name: file name of output
        model_path: path in which models are stored
        min_threshold: minimum softmax score needed to use ML prediction
    """
    #Create tuple of sequences and position for k-merized FASTA file
    evaluation_file = open(f"{report_path}/{report_name}_kmer.txt", "r")
    counter = 0
    seqs = []
    for sequence in evaluation_file:
            seqs.append((sequence[:-1], counter))
            counter += 1
    
    ranks = ['phylum', 'class', 'order', 'family', 'genus', "species"]
    model_preds = np.zeros(shape = (counter, 6))
    model_preds = pd.DataFrame(model_preds, columns = ranks)

    #Generate predictions for only those genera in which trained classifiers are present
    trained_genus = [model.split("_")[0] for model in os.listdir(model_path)]
    reference = pd.read_csv(f"{report_path}/{report_name}_lineage_kraken.csv")
    reference["genus"] = reference["genus"].apply(lambda x: x if x in trained_genus else "nan")
    model_preds["genus"] = reference["genus"]
    
    #Create dictionary with structure: {Prediction1 => [(Kraken2_Prediction, position in FASTA) ... (Kraken2_Prediction, position in FASTA)], Prediction2 => [(Kraken2_Prediction, position in FASTA) ... (Kraken2_Prediction, position in FASTA)] }
    genus_preds_dict = {}
    position = 0
    for genus in model_preds["genus"].values:
        pred = genus
        if not(pred in genus_preds_dict.keys()):
            genus_preds_dict[pred] = []
        genus_preds_dict[pred].append((pred, position))
        position +=1

    #Create list of precictons with structure: => [(prediction, score, position in FASTA)...(prediction, score, position in FASTA)]
    pred_tracker = []
    for genus, seq_list in genus_preds_dict.items():
        if str(genus) == "nan":
            for seq in seq_list:
                pred = "NA"
                score = 0.49
                pred_tracker.append((pred, score, seq[1]))
            continue
        else:
            model = ft.load_model(f"{model_path}/{genus}_model.bin")
            for seq in seq_list:
                model_prediction = model.predict(seqs[seq[1]][0])
                pred = model_prediction[0][0].split("label__")[1]
                score = model_prediction[1][0]
                pred_tracker.append((pred, score, seq[1]))
    #Resort predictions based on original position in FASTA
    pred_tracker = sorted(pred_tracker, key = lambda model_output: model_output[2])
    species_pred = []
    #Keep predictions which meet the minimum threshold
    for model_output in pred_tracker:
        if model_output[1] >= min_threshold:
            species_pred.append(model_output[0])
        else:
            species_pred.append("NA")
    model_preds["species"] = species_pred
    return model_preds

def evaluation_bwa(report_path, report_name, specialized_path):
    """
    Given Kraken2's genus classifications, generate species-level predictions using BWA
    Args:
        report_path: path to store classifer output
        report_name: file name of output
        model_path: path in which BWA indices are stored
    """
    #Create tuple of sequences and position for BWA FASTA file
    evaluation_file = open(f"{report_path}/{report_name}_bwa.fa", "r")
    counter = 0
    seqs = []
    for record in SeqIO.parse(evaluation_file, "fasta"):
        seqs.append((str(record.seq), counter))
        counter += 1

    ranks = ['phylum', 'class', 'order', 'family', 'genus', "species"]
    model_preds = np.zeros(shape = (counter, 6))
    model_preds = pd.DataFrame(model_preds, columns = ranks)

    #Generate predictions for only those genera in which BWA indices are present
    trained_genus = np.unique([indices.split(".")[0] for indices in os.listdir(specialized_path)])
    reference = pd.read_csv(f"{report_path}/{report_name}_lineage_kraken.csv")
    reference["genus"] = reference["genus"].apply(lambda x: x if x in trained_genus else "nan")
    model_preds["genus"] = reference["genus"]

    #Create dictionary with structure: {Prediction1 => [(Kraken2_Prediction, position in FASTA) ... (Kraken2_Prediction, position in FASTA)], Prediction2 => [(Kraken2_Prediction, position in FASTA) ... (Kraken2_Prediction, position in FASTA)] }
    genus_preds_dict = {}
    position = 0
    for genus in model_preds["genus"].values:
        pred = genus
        if not(pred in genus_preds_dict.keys()):
            genus_preds_dict[pred] = []
        genus_preds_dict[pred].append((pred, position))
        position +=1

    #Create list of precictons with structure: => [(prediction, score, position in FASTA)...(prediction, score, position in FASTA)]
    pred_tracker = []
    for genus, seq_list in genus_preds_dict.items():
        if str(genus) == "nan":
            for seq in seq_list:
                pred = "NA"
                score = 0.49
                pred_tracker.append((pred, score, seq[1]))
            continue
        else:
            temp_fasta = open(f"{report_path}/temp.fasta", "w")
            counters_seen = []
            for seq in seq_list:
                temp_fasta.write(">"+"_"+str(seq[1])+"\n"+seqs[seq[1]][0]+"\n")
                counters_seen.append(seq[1])
            temp_fasta.close()
            os.system(f"bwa mem -t 80 {specialized_path}/{genus}.fa {report_path}/temp.fasta > {report_path}/{report_name}_{genus}_aligned.sam")
            os.system(f"awk '($1 !~ /^@/ && $3 != \"*\") {{score = \"NA\"; for(i = 12; i <= NF; i++) {{if($i ~ /^AS:i:/) {{split($i, arr, \":\"); score = arr[3]; break;}}}} print $1 \"\t\" $3 \"\t\" score}}' {report_path}/{report_name}_{genus}_aligned.sam > {report_path}/{report_name}_{genus}_reads_to_reference.tsv")
            read_2_reference = pd.read_csv(f"{report_path}/{report_name}_{genus}_reads_to_reference.tsv", header = None, sep = "\t")
            read_2_reference["species"] = read_2_reference[1].apply(lambda x: x.split("|")[1])
            read_2_reference = read_2_reference.sort_values(2, ascending=False).drop_duplicates((0, "species")).sort_index()
            r2r_dict_holder = {}
            for index, row in read_2_reference.iterrows():
                if row[0] not in r2r_dict_holder.keys():
                          r2r_dict_holder[row[0]] = (row[1], row[2])
                else:
                          if row[2] >  r2r_dict_holder[row[0]][1]:
                                    r2r_dict_holder[row[0]] = (row[1], row[2])
                          elif row[2] == r2r_dict_holder[row[0]][1]:
                                    r2r_dict_holder[row[0]] = ("|NA", row[2])
                          else:
                                    continue
            new_read_2_reference = pd.DataFrame(r2r_dict_holder).T
            new_read_2_reference = new_read_2_reference.reset_index()
            new_read_2_reference  = new_read_2_reference.rename(columns ={"index":0, 0:1, 1:2})
            read_2_reference = new_read_2_reference
            positions_mapped = [int(z.split("_")[1])for z in read_2_reference[0].values]
            positions_missing = list(set(counters_seen) - set(positions_mapped))
            positions_dict = dict(zip(positions_mapped, read_2_reference[1].values))
            for b in positions_missing:
                positions_dict[b] = "|NA"
            for seq in seq_list:
                pred_tracker.append((positions_dict[seq[1]].split("|")[1], 1, seq[1]))

    #Resort predictions based on original position in FASTA
    pred_tracker = sorted(pred_tracker, key = lambda model_output: model_output[2])
    species_pred = []
    #Stores Predictions in List
    for model_output in pred_tracker:
            species_pred.append(model_output[0])
    model_preds["species"] = species_pred
    return model_preds

def ensemble(report_path, report_name, mode):
    """
    Ensemble Kraken2-informed specialized classifiers predictions and Kraken2 classifications
    Args:
        report_path: path to store classifer output
        report_name: file name of output
    """
    if mode == "Kraken2_ML":
        specialized_species = pd.read_csv(f"{report_path}/{report_name}_ml.csv")
    else:
        specialized_species = pd.read_csv(f"{report_path}/{report_name}_bwa.csv")
    reference = pd.read_csv(f"{report_path}/{report_name}_lineage_kraken.csv")
    ensemble_preds = []
    for i in range(0, len(reference)):
        #If Kraken2 has no genus predictions, ensemble leaves species unclassified
        if str(reference["genus"].iloc[i]) == "nan":
            ensemble_preds.append("nan")
        #IF specialized model has no species prediction, ensemble uses Kraken2 prediction
        elif str(specialized_species["species"].iloc[i]) == "nan":
            ensemble_preds.append(reference["species"].iloc[i])
        #By default, uses ML predictions
        else:
            ensemble_preds.append(str(specialized_species["species"].iloc[i]))
    ensemble_preds = [str(i).replace(" ","_") for i in ensemble_preds]
    reference["species"] = ensemble_preds
    return reference

