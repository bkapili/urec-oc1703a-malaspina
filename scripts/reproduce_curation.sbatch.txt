#! /bin/bash
#
#SBATCH --job-name=reproduce_curation
#
#SBATCH --partition=serc
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

bash reproduce_curation.sh
