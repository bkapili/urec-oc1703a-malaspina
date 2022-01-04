
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

## Overview

This repo contains the code to reproduce the UreC dataset curation from
OC1703A and Malaspina metagenome assemblies as performed in
Arandia-Gorostidi *et al.*, in prep. The sections below provide
instructions on how to:

  - Create a `conda` environment for data processing
  - Clone this repo
  - Download the metagenome assemblies
  - Execute the scripts

### conda environment setup and installations

The code below assumes you already have
[conda](https://docs.conda.io/en/latest/) installed. It will create a
new `conda` environment called `urec_oc1703a_malaspina`, activate it,
and install all the necessary packages (and the specific versions) used
during analysis. The environment with installations will require ~1.7 Gb
of disk space.

``` bash
# Set path for environment
CONDA_PATH=$(echo "SET_PATH_HERE")

# Create conda environment
conda create --prefix $CONDA_PATH/urec_oc1703a_malaspina
conda activate $CONDA_PATH/urec_oc1703a_malaspina

# Install packages
conda install -c bioconda blast=2.12.0  mafft=7.490 \
  hmmer=3.3.2 emboss=6.6.0 pal2nal=14.1
  
conda install -c conda-forge r-base=4.1.1
```

### Clone git repo

The code below will clone this repo into a subdirectory named
`urec-oc1703a-malaspina` in the directory you specify in `REPO_PATH`.

``` bash
# Set path for repo
REPO_PATH=$(echo "SET_PATH_HERE")

# Clone repo
mkdir $REPO_PATH/urec-oc1703a-malaspina
git clone https://github.com/bkapili/urec-oc1703a-malaspina.git $REPO_PATH/urec-oc1703a-malaspina
```

### Download metagenome assemblies

The code below will download the OC1703A and Malaspina metagenome
assemblies in the `data` directory. <span style="color: red;">**Update
to public repository link once data are
uploaded.**</span>

``` bash
cp $GROUP_SCRATCH/Nestor/OC1703_malaspina_assemblies.fasta $REPO_PATH/urec-oc1703a-malaspina/data/
```

### Run scripts

The code in this section will run the script `reproduce_curation.sh`,
which executes all the individual scripts in order. In order for it to
run properly, it should be executed directly from the scripts folder.
After it runs, a fasta file named `urec_aa_sequences.fasta` is saved in
a new subdirectory named `urec_filter` that contains the putative UreC
sequences. A csv file named `urec_database_details.csv` is also written
that contains detailed information about the UreC classification step
(*e.g.*, ORF start/stop positions, nhmmer hit start/stop positions, ORF
length, strand sense, *etc.*).

``` bash
# Change directory
cd $REPO_PATH/urec-oc1703a-malaspina/scripts

# Execute reproduce_curation.sh
bash reproduce_curation.sh
```

Additional information about each script is provided in the block of
comment code at the head of each script.
