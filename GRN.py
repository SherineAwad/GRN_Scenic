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


# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt = "Rbpj_mCherry_unfiltered.loom" # test dataset, n=500 cells

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "Rbpj_mCherry_unfiltered.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = "adata_prefilter.h5ad"

# path to pyscenic output
f_pyscenic_output = "pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = 'Rbpj_mCherry_scenic_integrated-output.loom'




f_mtx_dir =  wdir 

f_tfs = "allTFs_mm.txt"
cmd = f"pyscenic grn {f_loom_path_scenic} {f_tfs} -o Rbpj_mCherry_adj.csv --num_workers 4"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
print(result.stdout)  # Prints the standard output of the command
print(result.stderr)  # Prints any error messages


adjacencies = pd.read_csv("Rbpj_mCherry_adj.csv", index_col=False, sep='\t')
adjacencies.head()


