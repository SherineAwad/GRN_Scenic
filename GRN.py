import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess

wdir ="/nfs/turbo/umms-thahoang/sherine/scenic/"
os.chdir( wdir )

f_mtx_dir =  wdir 

f_tfs = "allTFs_mm.txt"
f_loom_path_scenic = "Rbpj_mCherry_filtered_scenic.loom"
#f_tfs = "test_tfs.txt"
output_csv = "Rbpj_mCherry_adj.csv"
num_workers = 4

cmd = [
    "pyscenic",
    "grn",
    f_loom_path_scenic,
    f_tfs,
    "-o", output_csv,
    "--num_workers", str(num_workers)
]


result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
print(result.stdout)  # Prints the standard output of the command
print(result.stderr)  # Prints any error messages


adjacencies = pd.read_csv("Rbpj_mCherry_adj.csv", index_col=False, sep='\t')
adjacencies.head()


