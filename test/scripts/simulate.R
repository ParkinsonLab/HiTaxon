library(polyester)
library(Biostrings)
library(stringr)

#Simulate test reads with Polyester
simulator <- function(sample_num, test_path) {
  sample_csv <- paste(paste(test_path,"/sample_", as.character(sample_num), sep = ""),".csv", sep = "")
  sample <- read.csv(sample_csv)
  my_range <- 1:nrow(sample)
  for(i in my_range) {
    genome <- sample["X0"][i,1]
    name <- sample["X1"][i,1]
    file_name <- (str_replace(name, " ", "_"))
    transcripts <- sample["num_transcript"][i,1]
    fold_change = sample(c(1, 1, 2, 0.5, 4, 0.25, 8, 0.125), size=transcripts, replace=TRUE, prob=c(0.34136827, 0.34136827, 0.13592719, 0.13592719, 0.02140428, 0.02140428, 0.00130026, 0.00130026))
    fasta = paste(sample["X0"][i,1], ".fna", sep = "")
    fasta_path = paste(test_path, file_name, fasta, sep = "/")
    fold_change_matrix = matrix(c(rep(1,transcripts), fold_change), nrow = transcripts)
    current_sample = paste(test_path,"/simulated_reads/", as.character(sample_num), sep = "")
    simulate_experiment(fasta_path, reads_per_transcript=sample["baselines"][i,1],
                      num_reps=c(1,1), fold_changes=fold_change, outdir = paste(current_sample, file_name, sep = "/"), readlen = 150, seed = sample_num + 1)
}
}



args <- commandArgs(trailingOnly = TRUE)
sample_range <- 1:10
for(sample_num in sample_range) {
  simulator(sample_num, args[1])
}

