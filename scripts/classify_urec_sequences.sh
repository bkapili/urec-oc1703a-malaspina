# -------------------------------------------------------------
# Script purpose: Identify putative UreC sequences from the
#                 possible ORFs.
#
# Inputs:
#   * Gold Standard ureC nucleotide and amino acid sequences
#   * .fasta of OC1703A and Malaspina assemblies
#
# Outputs:
#   * Fasta of putative UreC sequences ("urec_aa_sequences.fasta")
#   * Updated summary of nhmmer and getorf statistics for each extracted
#     sequence including which ORFs were identified as UreC
#     ("best_orfs_summary.csv")
# -------------------------------------------------------------

# Make BLAST db of Gold Standard UreC sequences
makeblastdb -in ../data/ureC_aa_gold_standard.fasta \
  -dbtype prot \
  -parse_seqids \
  -out ../data/blastdb/dbUreC

# Perform all-against-all BLAST between ORFs and BLAST db
mkdir ../data/urec_filter
cd ../data/urec_filter

blastp -db ../blastdb/dbUreC -query ../extracted_seqs/best_orfs.fasta -outfmt 6 \
  -max_target_seqs 999999 -out urec_blast_output.txt -num_threads 4
  
# Run R script that identifies putative UreC sequences
cd ../../scripts
Rscript filter_evalue.R
