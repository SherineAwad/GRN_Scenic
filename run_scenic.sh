#!/bin/bash
#SBATCH --job-name scenic 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=50000
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=sherinem@umich.edu
#SBATCH --account=thahoang99 

conda activate scenic 

python GRN.py 

#pyscenic grn Rbpj_mCherry_filtered_scenic.loom test_tfs.txt -o test_adj.csv --num_workers 2
