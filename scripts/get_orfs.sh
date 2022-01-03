# -------------------------------------------------------------
# Script purpose: Identify the ORF that overlaps best with the nhmmer
#                 hit on each extracted nucleotide sequence. 
#
# Inputs:
#   * Fasta file containing nucleotide sequences of nhmmer hits
#     ("raw_ureC_nucleotide_extracts.fasta")
#   * Summary of nhmmer hits ("nhmmer_hit_summary.txt")
#   * Fasta of ORFs inferred using Genetic Code 11 ("orfs.fasta")
#
# Outputs:
#   * Fasta file that contains the best ORF for each nhmmer
#     hit ("best_orfs.fasta").
#   * Summary of nhmmer and getorf statistics for each extracted
#     sequence ("best_orfs_summary.csv")
#
# Notes:  Best ORF identified as the shortest ORF that completely
#         encompasses the nhmmer hit. If there is no such ORF,
#         then the ORF that minimizes the sum of the absolute
#         differences between the ORF start and stop positions
#         and the nhmmer hit start and stop positions is chosen
#         as the best ORF.
# -------------------------------------------------------------

# Identify all ORFs
DIR_PATH=$(echo "../data/extracted_seqs")

sed -i '/^>/ s/:.*$//' $DIR_PATH/raw_ureC_nucleotide_extracts.fasta
getorf -sequence $DIR_PATH/raw_ureC_nucleotide_extracts.fasta \
  -table 11 -minsize 300 -find 1 -outseq $DIR_PATH/orfs.fasta
  
# Identify best ORF for each extracted sequence
Rscript filter_best_orfs.R
