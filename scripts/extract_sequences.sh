# -------------------------------------------------------------
# Script purpose: Extract nhmmer hit sequences. 
#
# Input:  Edited summary of nhmmer hits ("nhmmer_hit_summary.txt")
#
# Output: Fasta file containing nucleotide sequences of nhmmer hits
#         ("raw_ureC_nucleotide_extracts.fasta").
# -------------------------------------------------------------

# Make BLAST database
mkdir ../data/blastdb
makeblastdb -in ../data/OC1703_malaspina_assemblies.fasta \
  -dbtype nucl \
  -parse_seqids \
  -out ../data/blastdb/dbOM
  
# Extract sequences between extended region indices
mkdir ../data/extracted_seqs
cd ../data/extracted_seqs

while read hzp; do
  seqNAME=$(echo $hzp | cut -d ' ' -f1)
  START=$(echo $hzp | cut -d ' ' -f7)
  STOP=$(echo $hzp | cut -d ' ' -f8)

  blastdbcmd \
    -db ../blastdb/dbOM \
    -dbtype nucl \
    -entry $seqNAME \
    -range $START-$STOP >> raw_ureC_nucleotide_extracts.fasta
done < ../nhmmer/nhmmer_hit_summary.txt
