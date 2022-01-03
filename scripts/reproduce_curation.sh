# -------------------------------------------------------------
# Script purpose: Execute full pipeline for curating UreC dataset
#                 from OC1703A and Malaspina metagenome assemblies.
#
# Notes:          See individual scripts for more information
#                 about each step.
# -------------------------------------------------------------

# Run nhmmer
bash nhmmer.sh

# Extract sequences of nhmmer hits
bash extract_sequences.sh

# Identify ORFs on extracted sequences
bash get_orfs.sh

# Classify UreC sequences
bash classify_urec_sequences.sh
