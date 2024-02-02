# HiTaxon: A hierarchical ensemble framework for taxonomic classification of short reads

## Overview

**HiTaxon** is an automated framework for creating custom short-read taxonomic classifiers for various environments. Given a list of genera, HiTaxon downloads and processes assemblies from RefSeq such that it creates a set of non-redundant sequences for species encompassed in the genera. HiTaxon then creates a custom database for a primary reference-dependent classifier and uses the classifier to generate taxonomic predictions on input FASTA files. The reference dependent genus-level outputs are used to determine which HiTaxon built specialized classifiers for specific genus to set of species pairs are used for species predictions.

Please refer to our paper for more information:

> Verma, B. and Parkinson, J. HiTaxon: A hierarchical ensemble framework for taxonomic classification of short reads. Bioinformatics Advances. 2024. DOI: https://doi.org/10.1093/bioadv/vbae016

Data used to evaluate classifiers in our study can be accessed at: 

https://doi.org/10.5281/zenodo.8335901


## Installation

Clone the repository
```bash
git clone https://github.com/ParkinsonLab/HiTaxon
```
Navigate to the repository
```bash
cd HiTaxon
```
Create conda environment (Can take ~30 Minutes)
```bash
conda env create --name HiTaxon --file environment.yml
```
Activate conda environment 
```bash
conda activate HiTaxon
```
Build HiTaxon
```bash
pip install .
```



## Basic Setup

In order to use HiTaxon, you must create a configuration file using the structure below and save it as config.file in the HiTaxon directory

```
ASSEMBLY_SUMMARY=default
GENUS_NAMES=/path/to/list_of_genera
KRAKEN_NAME=desired_database_name
KRAKEN_PATH=/path/to/store_database
MODEL_PATH=/path/to/store_ML_models
NUM_OF_THREADS=number_of_threads
OUTPUT_PATH=/path/to/download_and_process_RefSeq_data
REPORT_PATH=/path/to/store_taxonomic_predictions
BWA_PATH=/path/to/store_bwa_indices
```

Note 1: If ASSEMBLY_SUMMARY is left as default, HiTaxon will download and use the latest file from NCBI. If you want to use a specific version, provide the path to the specified assembly_summary text file

Note 2: For directory path declarations (i.e variables with PATH in name), if directory does not exist, HiTaxon will create the directory 

The text file corresponding to GENUS_NAMES needs to be structured as below:

```
Bacillus
Enterococcus
Escherichia
Lactobacillus
Listeria
Staphylococcus
Salmonella
```

## Use Cases Higlighted in Publication

1. [Building the best taxonomic classifier for a particular dataset when there are moderate time constraints](#use-case-1-building-the-best-taxonomic-classifier-for-a-particular-dataset-moderate-time-constraints)
2. [Building the best taxonomic classifier for a particular dataset when there are significant time constrains](#use-case-2-building-the-best-taxonomic-classifier-for-a-particular-dataset-significant-time-constrains)

### Use Case 1: Building the best taxonomic classifier for a particular dataset (moderate time constraints)

To create the best taxonomic classifier for a particular dataset when there are moderate time constraints, we recommend employing Kraken2-HiTaxon-Align, which is an hierarchical ensemble consisting of a Kraken2 classifier paired with a HiTaxon constructed database and specialized genus to set of species BWA indices. In our publication, we highlight that Kraken2-HiTaxon-Align is the best performing taxonomic classifier amongst all that were tested.

In order to build this ensemble, execute the following steps in order:
1. Create a directory to store sequence data from RefSeq.
2. Within this directory, create a file called `taxon.txt` which lists all genera of interest, using the same format highlighted earlier in the documentation.
3. In the HiTaxon directory, create a `config.file` which defines a set of important parameters for HiTaxon, using the same format highlighted earlier in the documentation.
4. Download sequences from RefSeq that pertain to the set of genera listed in `taxon.txt`:
    ```bash
    ./HiTaxon.sh --collect
    ```
5. Cluster similiar assemblies and coding sequences:
    ```bash
    ./HiTaxon.sh --process
    ```
6. Create a custom database composed of non-redundant sequences aquired from Step 5 for Kraken2:
    ```bash
    ./HiTaxon.sh --build
    ```
7. Create specialized BWA indices for each genus to set of species pair using non-redundant sequences aquired from Step 5 for BWA:
    ```bash
    ./HiTaxon.sh --align
    ```
    Note: If execution is terminated mid-construction of index for specific genus to set of species pair, simply delete the .ann file corresponding to it and rerun the command 

8. Use the hierarchical ensemble to generate taxonomic predictions for short-reads in `input.fasta`:
    ```bash
    ./HiTaxon.sh --evaluate -f path/to/input.fasta -o name_of_output_report -m Kraken2_BWA
    ```
   Will generate an output of `{name_of_output_report}_ensemble_bwa.csv `

### Use Case 2: Building the best taxonomic classifier for a particular dataset (significant time constrains)
To create the best taxonomic classifier for a particular dataset when there are significant time constraints, we recommend employing Kraken2-HiTaxon-DB, which is consists solely of Kraken2 with a HiTaxon curated database. This approach has a small reduction in recall relative to Kraken2-HiTaxon-Align but requires much less time to generate predictions.

In order to build this classifier, execute the following steps in order:
1. Create a directory to store sequence data from RefSeq.
2. Within this directory, create a file called `taxon.txt` which lists all genera of interest, using the same format highlighted earlier in the documentation.
3. In the HiTaxon directory, create a `config.file` which defines a set of important parameters for HiTaxon, using the same format highlighted earlier in the documentation.
4. Download sequences from RefSeq pertaining to the set of genera listed in `taxon.txt`:
    ```bash
    ./HiTaxon.sh --collect
    ```
5. Cluster similiar assemblies and coding sequences:
    ```bash
    ./HiTaxon.sh --process
    ```
6. Create a custom database composed of non-redundant sequences aquired from Step 5 for Kraken2:
    ```bash
    ./HiTaxon.sh --build
    ```

7. Use the hierarchical ensemble to generate taxonomic predictions for short-reads in `input.fasta`:
    ```bash
    ./HiTaxon.sh --evaluate -f path/to/input.fasta -o name_of_output_report -m Kraken2
    ```
    Will generate an output of `{name_of_output_report}_lineage_kraken.csv` 

## Potential Other Usages

1. [Ensembling existing Kraken2 databases with HiTaxon-curated specialized classifiers](#potential-other-use-1-ensembling-existing-kraken2-database-with-hitaxon-curated-specialized-classifiers)
2. [Creating Training Data for custom ML classifiers](#potential-other-use-2-creating-training-data-for-custom-ml-classifier)

### Potential Other Use 1: Ensembling existing Kraken2 database with HiTaxon-curated specialized classifiers
To create a hierarchical ensemble using an existing Kraken2 database with HiTaxon-curated specialized classifiers, execute the following steps in order:

1. Create a directory to store sequence data from RefSeq.
2. Within this directory, create a file called `taxon.txt` which lists all genera of interest, using the same format highlighted earlier in the documentation.
3. In the HiTaxon directory, create a `config.file`, using the same format highlighted earlier in the documentation. **Make sure config.file references the path in which the pre-existing Kraken2 database is stored, alongside the name of the database**
4. Download sequences from RefSeq pertaining to the set of genera listed in `taxon.txt`:
    ```bash
    ./HiTaxon.sh --collect
    ```
5. Cluster similiar assemblies and coding sequences:
    ```bash
    ./HiTaxon.sh --process
    ```
6. Option A: Create specialized BWA indices for each genus to set of species pair using non-redundant sequences aquired from Step 5 for BWA:
    ```bash
    ./HiTaxon.sh --align
    ```

   Option B. Train specialized ML for each genus to set of species pair using non-redundant sequences aquired from Step 5:
     ```bash
    ./HiTaxon.sh --train
    ```
   Note: You will be prompted by the command line on as to whether the HiTaxon is in 1) data_creation  or 2) model_creation mode, **select option 2** . 
  


7. Option A. Use the hierarchical ensemble of Kraken2 and BWA to generate taxonomic predictions for short-reads in `input.fasta`:
    ```bash
    ./HiTaxon.sh --evaluate -f path/to/input.fasta -o name_of_output_report -m Kraken2_BWA
    ```

    Option B. Use the hierarchical ensemble of Kraken2 and ML classifiers to generate taxonomic predictions for short-reads in  `input.fasta`:
    ```bash
    ./HiTaxon.sh --evaluate -f path/to/input.fasta -o name_of_output_report -m Kraken2_ML
    ```

### Potential Other Use 2: Creating Training Data for Custom ML Classifier
To train ML classifiers for HiTaxon, we have built a robust framework for training both multi-class classifiers (i.e genus encompasses multiple species) and binary classifiers (i.e genus encompasses a single species). However, while we used FastText classifiers, we do acknowledge that researchers might prefer to test different algorithms for species classification. Consequently, we allow for HiTaxon to create and process training data without forcing the user to also train ML models.
To use this feature, execute the following steps in order:

1. Create a directory to store sequence data from RefSeq.
2. Within this directory, create a file called `taxon.txt` which lists all genera of interest, using the same format highlighted earlier in the documentation.
3. In the HiTaxon directory, create a `config.file`, using the same format highlighted earlier in the documentation.
4. Download sequences from RefSeq pertaining to the set of genera listed in `taxon.txt`:
    ```bash
    ./HiTaxon.sh --collect
    ```
5. Cluster similiar assemblies and coding sequences:
    ```bash
    ./HiTaxon.sh --process
    ```
6. Create and preprocess the training data 
    ```bash
    ./HiTaxon.sh --train
    ```
You will be prompted whether HiTaxon is in 1) data_creation  or  2) model_creation mode, **select  option 1**. This will generate a .txt for all genus to set of species pairs, for both multi-class and binary problems, which can be easily reformatted by researchers for their ML algorithm of choice.

Note: If you want to change the K-mer value from the default of K = 13, a single parameter change in train.py, which can be found in scripts/train/, is all that is needed


## Links to software used for HiTaxon
- [art](https://anaconda.org/bioconda/art)
- [biopython](https://github.com/biopython/biopython)
- [bwa](https://github.com/lh3/bwa)
- [cd-hit](https://github.com/weizhongli/cdhit)
- [ete](https://github.com/etetoolkit/ete)
- [FastANI](https://github.com/ParBLiSS/FastANI)
- [fastp](https://github.com/OpenGene/fastp)
- [fastq-join](https://github.com/brwnj/fastq-join)
- [fasttext](https://github.com/facebookresearch/fastText)
- [kaiju](https://github.com/bioinformatics-centre/kaiju)
- [kraken2](https://github.com/DerrickWood/kraken2)
- [ncbi datasets](https://github.com/ncbi/datasets)
- [numpy](https://github.com/numpy/numpy)
- [pandas](https://github.com/pandas-dev/pandas)
- [polyester](https://github.com/alyssafrazee/polyester)
- [prodigal](https://github.com/hyattpd/Prodigal)
- [scikit-learn](https://github.com/scikit-learn/scikit-learn)
- [seqkit](https://github.com/shenwei356/seqkit)
- [seqtk](https://github.com/lh3/seqtk)
