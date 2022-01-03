# -------------------------------------------------------------
# Script purpose: Identify E-value threshold for classification
#                 of UreC sequences.
#
# Inputs:
#   * Gold Standard ureC nucleotide and amino acid sequences
#   * .fasta of OC1703A and Malaspina assemblies
#
# Outputs:
#   * Summary of nhmmer hits ("nhmmer_hit_summary.txt")
# -------------------------------------------------------------

### Load required packages
# List required packages
cranPackages <- c("BiocManager", "dplyr", "tidyr", "ggplot2", "Biostrings")
biocPackages <- c("Biostrings")

# Install missing CRAN packages
installedCRANPackages <- cranPackages %in% rownames(installed.packages())
if (any(installedCRANPackages == FALSE)) {
  install.packages(cranPackages[!installedCRANPackages],
                   repos='http://cran.us.r-project.org')
}

# Install missing Bioconductor packages
installedBioPackages <- biocPackages %in% rownames(installed.packages())
if (any(installedBioPackages == FALSE)) {
  BiocManager::install(biocPackages[!installedBioPackages])
}

# Load packages
lapply(c(cranPackages, biocPackages), library, character.only = TRUE)

### Identify reasonable E-value threshold/cutoff
# Set minimum length requirement
gs <- readAAStringSet(filepath = "../data/ureC_aa_gold_standard.fasta")
minLength <- min(width(gs))/2

# Load BLAST results
dfBlast <- read.table(file = "../data/urec_filter/urec_blast_output.txt")
colnames(dfBlast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                       "qstart", "qend", "sstart", "send", "evalue", "bitscore")
dfBlast <- dfBlast %>%
  mutate(log_e = log(dfBlast$evalue, base = 10))

# Plot log10-transformed E-values for BLAST hits
# with alignment length greater than half the shortest
# UreC in the Gold Standard
dfPlot <- dfBlast %>%
  filter(length > minLength)

ggplot(data = dfPlot, aes(x = log_e)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = -70, linetype = "dashed") +
  theme_bw()

### Filter sequences using threshold
# Keep all BLAST hits with log10(E-value) < -70
keep <- dfBlast %>%
  filter(log_e < -70) %>%
  pull(qseqid) %>%
  unique

bestOrf <- readAAStringSet(filepath = "../data/extracted_seqs/best_orfs.fasta")
urec <- bestOrf[names(bestOrf) %in% keep]
urec <- urec[order(width(urec))]

### Export results
# Write putative UreC sequences to fasta
writeXStringSet(urec, filepath = "../data/urec_filter/urec_aa_sequences.fasta")

# Update log to track which ORFs were classified as UreC
dfBestOrfs <- read.csv(file = "../data/extracted_seqs/best_orfs_summary.csv", row.names = 1)
dfBestOrfs <- dfBestOrfs %>%
  mutate(urec = seq_names %in% names(urec))

write.csv(dfBestOrfs, file = "../data/urec_filter/urec_database_details.csv", quote = FALSE)