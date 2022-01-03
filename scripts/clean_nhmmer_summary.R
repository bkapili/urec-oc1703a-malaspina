# -------------------------------------------------------------
# Script purpose: Prepare table specifying sequences ranges to extract.
#
# Input:    Summary of nhmmer hits ("nhmmer_hit_summary.txt")
#
# Outputs:  Edited summary of nhmmer hits in which start and stop
#           HMM hit indices are swapped if on negative strand
#           and two additional columns indicating the sequence range
#           in the assembly to extract. Extraction indices are the HMM
#           hits extended by 2000 bp in both directions to ensure
#           the full gene is extracted ("nhmmer_hit_summary.txt").
# -------------------------------------------------------------

### Load required packages
# List required packages
cranPackages <- c("BiocManager", "dplyr", "tidyr")

# Install missing CRAN packages
installedCRANPackages <- cranPackages %in% rownames(installed.packages())
if (any(installedCRANPackages == FALSE)) {
  install.packages(cranPackages[!installedCRANPackages],
                   repos='http://cran.us.r-project.org')
}

# Load packages
lapply(c(cranPackages), library, character.only = TRUE)

### Edit nhmmer hit summary file
# Read reformated tblout
filePath <- "../data/nhmmer/nhmmer_hit_summary.txt"
df <- read.table(file = filePath)
colnames(df) <- c("name", "ali_from", "ali_to", "seq_len", "strand", "e_value")

# Swap start and stop for hits on negative strand
indSwap <- df$strand == "-"
df[indSwap, c("ali_from", "ali_to")] <- df[indSwap, c("ali_to", "ali_from")]

# Extend the region to extract on each side by 2000 nucleotides
# to ensure full gene is extracted
df <- df %>%
  mutate(extend_start = ali_from-2000,
         extend_stop = ali_to+2000)

# Set negative start indices to 1
df[which(df[, "extend_start"] < 1), "extend_start"] <- 1

# Set stop indices that are greater than seq length to equal seq length
setStop <- which(df[,"extend_stop"] > df[, "seq_len"])
df[setStop, "extend_stop"] <- df[setStop, "seq_len"]

# Write to table
write.table(df, file = filePath, sep = " ",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
