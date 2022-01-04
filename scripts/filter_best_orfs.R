# -------------------------------------------------------------
# Script purpose: Identify the ORF that overlaps best with each
#                 nhmmer hit.
#
# Inputs:
#   * Summary of nhmmer hits ("nhmmer_hit_summary.txt")
#   * Fasta of ORFs ("orfs.fasta")
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

library(Biostrings)
library(dplyr)

# Load nhmmer summary
dfSummary <- read.table(file = "../data/nhmmer/nhmmer_hit_summary.txt")
colnames(dfSummary) <- c("name", "nhmmer_ali_from", "nhmmer_ali_to", "seq_len", "strand",
                         "e_value", "extend_start", "extend_stop")
dfSummary$extend_start <- dfSummary$extend_start - 1

# Load ORFs
orfs <- readAAStringSet(filepath = "../data/extracted_seqs/orfs.fasta")

### Create tibble of ORF data
# Extract sequence names
seqNames <- names(orfs) %>%
  strsplit(" ") %>%
  sapply("[[", 1) %>%
  sub("[_][^_]+$", "", .)

# Extract ORF start indices
orf_start <- names(orfs) %>%
  strsplit(" ") %>%
  sapply("[[", 2) %>%
  gsub("\\[", "", .) %>%
  as.integer

# Extract ORF stop indices
orf_stop <- names(orfs) %>%
  strsplit(" ") %>%
  sapply("[[", 4) %>%
  gsub("\\]", "", .) %>%
  as.integer

# Match sequences to corresponding rows in nhmmer
# output table to determine num of bp to add to account
# for bases added during extraction
matchToSummary <- match(seqNames, dfSummary$name)

# Assemble into data frame
dfOrfs <- data.frame(orf_names = names(orfs),
                     seq_names = seqNames,
                     orf_start = orf_start + dfSummary$extend_start[matchToSummary],
                     orf_stop = orf_stop + dfSummary$extend_start[matchToSummary])

# Swap start and stop for ORFs on negative strand
indSwap <- grep("REVERSE", dfOrfs$orf_names)
dfOrfs[indSwap, c("orf_start", "orf_stop")] <- dfOrfs[indSwap, c("orf_stop", "orf_start")]

# Add nhmmer summary to ORF data frame
matchToSummary <- match(dfOrfs$seq_names, dfSummary$name)
dfSummaryJoin <- dfSummary[matchToSummary,] %>%
  select(nhmmer_ali_from, nhmmer_ali_to, seq_len, strand)

dfOrfs <- cbind(dfOrfs, dfSummaryJoin)

# Remove ORFs found on opposite strand of nhmmer hit
dfOrfs <- dfOrfs %>%
  filter(case_when(strand == "-" ~ grepl("REVERSE", orf_names),
                   strand == "+" ~ !grepl("REVERSE", orf_names)))

### Determine best ORF hits
# Determine which ORFs fully contain nhmmer hit
matchToSummary <- match(dfOrfs$seq_names, dfSummary$name)
completeCoverage <- dfOrfs$orf_start <= dfSummary$nhmmer_ali_from[matchToSummary] &
  dfOrfs$orf_stop >= dfSummary$nhmmer_ali_to[matchToSummary]

delta <- abs(dfOrfs$orf_start - dfSummary$nhmmer_ali_from[matchToSummary]) +
  abs(dfOrfs$orf_stop - dfSummary$nhmmer_ali_to[matchToSummary])

dfOrfs <- dfOrfs %>%
  mutate(orf_length = orf_stop-orf_start,
         complete_coverage = completeCoverage,
         delta = delta)

# Identify shortest fully-containing ORF for each sequence
dfCompleteCov <- dfOrfs %>%
  group_by(seq_names) %>%
  filter(complete_coverage == TRUE) %>%
  slice(which.min(orf_length))

# For sequences with no ORF fully containing the nhmmer hit,
# identify ORF that minimizes absolute differences between
# ORF and nhmmer hit start and stop positions
dfMinDelta <- dfOrfs %>%
  group_by(seq_names) %>%
  filter(!any(complete_coverage == TRUE)) %>%
  slice(which.min(delta))

# Combine best ORFs
dfBestOrfs <- rbind(dfCompleteCov, dfMinDelta)
dfOrfs[dfOrfs$seq_names %in% dfBestOrfs$seq_names == FALSE,]

# Write to fasta
bestOrfs <- orfs[match(dfBestOrfs$orf_names, names(orfs))]
matchInd <- match(names(bestOrfs), dfBestOrfs$orf_names)
names(bestOrfs) <- dfBestOrfs$seq_names[matchInd]
writeXStringSet(bestOrfs, filepath = "../data/extracted_seqs/best_orfs.fasta")

# Write summary of best ORFs to .csv
matchToSummary <- match(dfBestOrfs$seq_names, dfSummary$name)
dfSummaryJoin <- dfSummary[matchToSummary,] %>%
  select(nhmmer_ali_from, nhmmer_ali_to, seq_len, strand)

dfBestOrfs <- cbind(dfBestOrfs, dfSummaryJoin)
write.csv(dfBestOrfs, file = "../data/extracted_seqs/best_orfs_summary.csv", quote = FALSE)
