# -------------------------------------------------------------
# Script purpose: Identify ureC homologs in OC1703A and Malaspina assemblies
#                 using nhmmer.
#
# Inputs:
#   * Gold Standard ureC nucleotide and amino acid sequences
#   * Fasta file of OC1703A and Malaspina assemblies
#
# Outputs:
#   * Summary of nhmmer hits ("nhmmer_hit_summary.txt")
# -------------------------------------------------------------

# Make and change directory
mkdir ../data/nhmmer
cd ../data/nhmmer

# Create amino acid alignment of UreC gold standards
mafft --maxiterate 1000 --thread 8 --dash --localpair --originalseqonly \
  ../ureC_aa_gold_standard.fasta > ureC_aa_alignment.fasta
  
# Backalign ureC nucleoide sequences to AA alignment
pal2nal.pl ureC_aa_alignment.fasta ../ureC_nt_gold_standard.fasta \
  -output fasta -codontable 11 > ureC_nt_alignment.fasta
  
# Shorten sequence names so they can be formatted in BLAST db
sed -i '/^>/ s/\s.*$//' ../OC1703_malaspina_assemblies.fasta
sed -i '/^>/ s/OC1703_OC1703/OC1703/g' ../OC1703_malaspina_assemblies.fasta
maxLen=$(grep '>' ../OC1703_malaspina_assemblies.fasta | wc -L)

if [[ "$maxLen" -gt 50 ]]; then
  echo "Maximum sequence name length too long (>50 chars)."
  exit; fi
  
# Run nhmmer
nhmmer --dna --qformat afa -E 10.0 --incE 1e-2 --seed 767620 \
  -o nhmmer_output.txt --hmmout nhmmer_hmmout.hmm \
  --tblout nhmmer_tblout.txt --notextw --cpu 4 \
  ureC_nt_alignment.fasta ../OC1703_malaspina_assemblies.fasta
  
# Extract sequence name, alignment start index, alignment stop index,
# sequence length, strand, and E-value of nhmmer hits
sed 's/  */ /g' nhmmer_tblout.txt |\
  cut -d ' ' -f1,7,8,11,12,13 |\
  sed '/^#/d' > nhmmer_hit_summary.txt
  
# Run R script to edit nhmmer_hit_summary.txt in preparation
# for sequence extraction
cd ../../scripts
Rscript clean_nhmmer_summary.R
