import numpy as np
import pandas as pd
import fasttext as ft

from tqdm import tqdm

def evaluate_taxa(taxa, predictions_dict, predictions_df, model_path):
    """
    Generate taxonomic predictions for the next rank for all reads predicted to be of 'taxa' at current rank
    Args:
        taxa: predicted taxa
        predictions_dict: dictionary with predictions as keys, and a list composed of tuples of (sequence, position)
        prediction_df: dataframe of predictions with scores
        model_path: path in which models are saved
    """
    model = ft.load_model(f"{model_path}/{taxa}_model.bin")
    for seq in predictions_dict[taxa]:
        prediction = model.predict(seq[0])
        pred_label = prediction[0][0].split("label__")[1]
        score = prediction[1][0]
        predictions_df[taxa].iloc[seq[1]] = (pred_label, score)
        predictions_dict[pred_label].append(seq)
        return predictions_df, predictions_dict

def lcpn_prediction(lineages, sequence_file, report_path, report_name, model_path):
    """ 
    Generate predictions for all reads across all taxonomic ranks using LCPN architecture
    Args:
        lineages: dataframe containing entire lineage of all trained species
        sequence_file: K-merized sequence file to analyze
        report_path: path in which outputs are saved
        report_name: output name
        model_path: path in which models are stored
    """
    ranks = ["species", "genus", "family", "class", "order", "phylum"]
    all_taxa = ["phylum"]

    #create list of all taxa from phylum to species
    for rank in ranks[::-1]:
        all_taxa.extend(list(np.unique(lineages[rank].values)))
        if rank == "genus":
            #count total number of taxa, exluding species
            num_of_taxa = len(all_taxa)
    all_taxa = [taxa.replace(" ","_") for taxa in all_taxa]

    predictions = {taxa: [] for taxa in all_taxa}
    #binary classifiers can have negative class corresponding to "different"
    predictions["different"] = []
     
    #Generate tuples of sequence and positional counter
    evaluation_file = open(sequence_file, "r")
    counter = 0
    seqs = []
    for sequence in evaluation_file:
            seqs.append((sequence[:-1], counter))
            counter += 1
    arr = np.zeros(shape = (counter, len(all_taxa)))
    predictions_mtx = pd.DataFrame(arr, columns = all_taxa)
    
    #All sequences need to be classified by the initial phylum classifier
    predictions["phylum"] = seqs
    #Generate predictions for each predicted taxa
    for taxa in tqdm(all_taxa[:num_of_taxa]):
            predictions_mtx, predictions = evaluate_taxa(taxa, predictions, predictions_mtx, model_path)

    return predictions_mtx

def standard_lcl_prediction(sequence_file, report_path, report_name, model_path):
    """
    Generate predictions for all reads using standard LCL architecture
    Args:
        sequence_file: K-merized sequence file to analyze
        report_path: path in which outputs are saved
        report_name: output name
        model_path: path in which models are stored    
    """
    ranks = ["species", "genus", "family", "class", "order", "phylum"]
    #Generate tuples of sequence and positional counter
    evaluation_file = open(sequence_file, "r")
    counter = 0
    seqs = []
    for sequence in evaluation_file:
            seqs.append((sequence[:-1], counter))
            counter += 1
    arr = np.zeros(shape = (counter, len(ranks)))
    predictions_mtx = pd.DataFrame(arr, columns = ranks)
    #Create predictions across all ranks
    for rank in ranks:
        model = ft.load_model(f"{model_path}/{rank}_model.bin")
        for seq in tqdm(seqs):
            prediction = model.predict(seq[0])
            pred_label = prediction[0][0].split("label__")[1]
            score = prediction[1][0]
            predictions_mtx[rank].iloc[seq[1]] = (pred_label, score)
    return predictions_mtx


def hierarchical_lcl_prediction(sequence_file, report_path, report_name, model_path):
    """
    Generate predictions for all possible outputs across all ranks using the LCL architecture with heirarchy-informed evaluation
    Args:
        sequence_file: K-merized sequence file to analyze
        report_path: path in which outputs are saved
        report_name: output name
        model_path: path in which models are stored    
    """
    ranks = ["species", "genus", "family", "class", "order", "phylum"]
    evaluation_file = open(sequence_file, "r")
    #Generate tuples of sequence and positional counter
    counter = 0
    seqs = []
    for sequence in evaluation_file:
            seqs.append((sequence[:-1], counter))
            counter += 1
    #Create predictions across all ranks
    for rank in ranks:
        model = ft.load_model(f"{model_path}/{rank}_model.bin")
        num_labels = len(model.labels)
        predictions = []
        for seq in tqdm(seqs):
            all_label_outputs = {}
            #Keep softmax score for all possible labels 
            for label in model.labels:
                all_label_outputs[label] = []
            prediction = model.predict(seq[0], k = num_labels)
            labels_and_scores = list(zip(prediction[0], prediction[1]))
            #Create dict of structure {label A: score, label B:score}
            for label_and_score in labels_and_scores:
                all_label_outputs[label_and_score[0]] = (label_and_score[1])
            predictions.append(all_label_outputs)

        #Create CSV file for all possible label at a specific rank
        predictions_mtx = pd.DataFrame(predictions)
        predictions_mtx.to_csv(f"{report_path}/{report_name}_{rank}_lcl.csv")


def lcl_lcpn_prediction(sequence_file, report_path, report_name, model_path):
    """
    Generate predictions for all reads across all taxonomic ranks using the LCL-LCPN architecture
        sequence_file: K-merized sequence file to analyze
        report_path: path in which outputs are saved
        report_name: output name
        model_path: path in which models are stored 
    """
    evaluation_file = open(sequence_file, "r")
    #Generate tuples of sequence and positional counter
    counter = 0
    seqs = []
    for sequence in evaluation_file:
            seqs.append((sequence[:-1], counter))
            counter += 1
    ranks = ['phylum', 'class', 'order', 'family', 'genus', "species"]
    arr = np.zeros(shape = (counter, len(ranks)))
    predictions_mtx = pd.DataFrame(arr, columns = ranks)
    for rank in ranks:
        #if rank is not genus or species, use the LCL framework for classification
        if rank in ['phylum', 'class', 'order', 'family']: 
            model = ft.load_model(f"{model_path}/{rank}_model.bin")
            #Predictions for each rank are temporarily stored in prediction_tracker
            prediction_tracker = []
            for seq in seqs:
                prediction = model.predict(seq[0])
                pred_label = prediction[0][0].split("label__")[1]
                score = prediction[1][0]
                prediction_tracker.append((pred_label, score))
            predictions_mtx[rank] = prediction_tracker
        #Using lCPN framework for genus and species classifications
        else:
            #Create dictionary of structure { prediction: [(prediction, position),...,(prediction, position)]} to keep track of position
            prior_tracker = {}
            position = 0
            for output in prediction_tracker:
                pred_label = output[0]
                if not(pred_label in prior_tracker.keys()):
                    prior_tracker[pred_label] = []
                prior_tracker[pred_label].append((pred_label, position))
                position +=1
            prediction_tracker = []
            for prior_pred, pred_and_position in prior_tracker.items():
                #Output of "different" can arise in instances where binary classifiers are used for family -> genus predictions
                if prior_pred == "different":
                    for single_pred_and_position in pred_and_position:
                        pred_label = "na"
                        score = 1.0
                        prediction_tracker.append((pred_label, score, single_pred_and_position[1]))
                    continue
                else:
                    model = ft.load_model(f"{model_path}/{prior_pred}_model.bin")
                for single_pred_and_position in pred_and_position:
                    prediction = model.predict(seqs[single_pred_and_position[1]][0])
                    pred_label = prediction[0][0].split("label__")[1]
                    score = prediction[1][0]
                    prediction_tracker.append((pred_label, score, single_pred_and_position[1]))
            #Re-sort list of predictions to ensure it matches the ordering of the original file
            prediction_tracker = sorted(prediction_tracker, key = lambda model_output: model_output[2])
            predictions_mtx[rank]= [(pred_and_score[0], pred_and_score[1]) for pred_and_score in prediction_tracker]
    return predictions_mtx

def lcpn_evaluation(output, report_path, report_name, threshold = 0.5):
    """
    Inference scheme of predictions of LCPN output
    Args:
        output_path: path in which RefSeq sequnces are downloaded and processed
        report_path: path in which outputs are saved
        report_name: name of output file
        threshold: minimum softmax threshold needed to accept prediction
    """
    #Load data
    lcpn_output = pd.read_csv(f"{report_path}/{output}").drop(columns='Unnamed: 0')
    predictions = ["unclassified" for i in range(len(lcpn_output))]
    to_ignore = [0, "0.0"]
    #Evaluate predictions bottom-up
    for i in range(len(predictions)):
        single_prediction = list(lcpn_output.loc[i].values)
        for rank_prediction in single_prediction:
            if rank_prediction in to_ignore:
                continue
            #Get softmax score
            num = float(rank_prediction.split(",")[1][1:-1])
            #Get prediction
            name = rank_prediction.split(",")[0][2:-1]
            if num >= threshold and name != "different":
                predictions[i] = name
    #Write predictions to file
    f = open(f"{report_path}/{report_name}_lcpn_predictions_{threshold}.txt", "w")
    for prediction in predictions:
        f.write(prediction + "\n")
    f.close()


def standard_lcl_evaluation(output, report_path, report_name, threshold = 0.5):
    """
    Inference schehem of LCL output
    Args:
        output_path: path in which RefSeq sequances are downloaded and processed
        report_path: path in which outputs are saved
        report_name: name of output file
        threshold: minimum softmax threshold needed to accept prediction
     """
    #Load data
    lcl_output = pd.read_csv(f"{report_path}/{output}_lcl.csv").drop(columns='Unnamed: 0')
    predictions = ["unclassified" for i in range(len(lcl_output))]
    ranks = ["species", "genus", "family", "order", "class", "phylum"]
    #Evaluate predictions bottom-up
    for rank in ranks:
        for i in range(len(predictions)):
            #get softmax score
            num = float(lcl_output[rank].iloc[i].split(",")[1][1:-1])
            if float(num) >= threshold and predictions[i] == "unclassified":
                predictions[i] = lcl_output[rank].iloc[i].split(",")[0][2:-1]
    #Write predictions to file
    f = open(f"{report_path}/{report_name}_lcl_predictions_{threshold}.txt", "w")
    for prediction in predictions:
        f.write(prediction + "\n")
    f.close()



def hierarchical_lcl_evaluation(lineage, output_path, species_record, report_path, report_name, threshold = 0.5):
    """
    Inference scheme of Hierarchy-informed LCL output
    Args:
        lineage: dataframe containing entire lineage of all trained species
        output_path: path in which RefSeq sequences are downloaded and processed
        species_record: text file of species 
        report_path: path in which outputs are saved
        report_name: name of output file
        threshold: minimum softmax threshold needed to accept prediction
    """
    ranks = ["phylum", "class", "order", "family", "genus", "species"]
    #List to store processed predictions across all ranks
    data_for_all_ranks = []

    for rank in ranks:
        current_rank = pd.read_csv(f"{report_path}/{report_name}_{rank}_lcl.csv", index_col = 0)
        if rank == "phylum":
            lineage[rank] = lineage[rank].apply(lambda x: "__label__" + str(x))
            #Creating list with structure [(label corresponding to highest softmax score, position in dataframe, softmax score)...]
            predictions =  list(zip(current_rank.idxmax(axis = 1).values, current_rank.index.values, current_rank.max(axis = 1).values))
            prior = rank
        else:
            data_for_specific_rank = []
            current = rank
            #Create map of taxa at prior ranks to a list of descent taxa at current rank
            prior2current = {}
            if rank != "species":
                lineage[rank] = lineage[rank].apply(lambda x: "__label__" + str(x))
                for taxa in lineage[prior].values:
                    lowered = list(np.unique(lineage[lineage[prior] == taxa][current].values))
                    prior2current[taxa] = lowered
            else:
                species_data = [species.replace(" ","_") for species in species_record]
                genus_data = []
                for species in species_data:
                    genus_data.append("__label__" + species.split("_")[0])
                species_data = ["__label__"  + i for i in species_data]
                genus2species = list(zip(genus_data, species_data))
                genus2species = pd.DataFrame(genus2species, columns = ["genus", "species"])
                for taxa in lineage[prior].values:
                    lowered = list(genus2species[genus2species[prior] == taxa][current].values)
                    prior2current[taxa] = lowered
            for prior_taxa,current_taxa in prior2current.items():
                preds_labels = np.array([i[0] for i in predictions])
                preds_idx = np.array([i[1] for i in predictions])
                #Subset for the position of predictions that were predicted to be == prior_taxa at the higher rank
                reduced_pred = preds_idx[preds_labels == prior_taxa]
                if len(reduced_pred) == 0:
                   continue
                #With these indexs, subset into the predictions corresponding to descent taxa at the current rank
                subset_df = current_rank[current_taxa].iloc[reduced_pred]
                #Compute a re-normalized softmax score for these subset of predictions 
                subset_df = subset_df.div(subset_df.sum(axis=1), axis=0)
                #[(label corresponding to highest softmax score, position in dataframe, softmax score)]
                new_preds = list(zip(current_rank[current_taxa].iloc[reduced_pred].idxmax(axis = 1), reduced_pred, subset_df.max(axis = 1).values))
                data_for_specific_rank.extend(new_preds)
        prior = rank
        if rank != "phylum":
            predictions = data_for_specific_rank
        data_for_all_ranks.append(predictions)
    sorted_data_for_all_ranks = []
    #Resort predictions at each rank to match the positions of the original sequence file
    for i in data_for_all_ranks:
        sorted_data_for_all_ranks.append(sorted(i, key=lambda x: x[1]))
    #Evaluate predictions bottom up
    pred_dict = {i:"nan" for i in range(0, len(predictions))}
    for i in sorted_data_for_all_ranks:
        for j in i:
            if j[-1] >= threshold:
               pred_dict[j[1]] = j[0]
    #Write predictions to file
    pred_df = pd.DataFrame(pred_dict, index = [0]).T
    pred_df[0] = pred_df[0].apply(lambda x: x.replace("__label__",""))
    f = open(f"{report_path}/hic_lcl_{report_name}.txt", "w")
    for d in pred_df[0].values:
        f.write(d + "\n")
    f.close()


def lcl_lcpn_evaluation(output, report_path, report_name, threshold = 0.5):
    """
    Inference scheme of the LCL-LCPN output
    Args:
        output_path: path in which downloaded and processed assembly data is stored
        report_path: path in which outputs are saved
        report_name: name of output file
        threshold: minimum softmax threshold needed to accept prediction
    """
    #Load Data
    lcl_lcpn_output = pd.read_csv(f"{report_path}/{output}_lcl_lcpn.csv").drop(columns='Unnamed: 0')
    predictions = ["unclassified" for i in range(len(lcl_lcpn_output ))]
    ranks = ["species", "genus", "family", "order", "class", "phylum"]
    #Evaluate predictions bottom-up
    for rank in ranks:
        for i in range(len(predictions)):
            #get softmax score
            num = float(lcl_lcpn_output[rank].iloc[i].split(",")[1][1:-1])
            if float(num) >= threshold and predictions[i] == "unclassified" and lcl_lcpn_output[rank].iloc[i].split(",")[0][2:-1] != "different":
                if rank != "phylum" and lcl_lcpn_output[rank].iloc[i].split(",")[0][2:-1] == "na":
                    continue
                predictions[i] = lcl_lcpn_output[rank].iloc[i].split(",")[0][2:-1]
    #Write predictions to file
    f = open(f"{report_path}/{report_name}_lcl_lcpn_predictions_{threshold}.txt", "w")
    for prediction in predictions:
        f.write(prediction + "\n")
    f.close()

