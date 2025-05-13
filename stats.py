import anndata as ad

adata = ad.read_loom("Rbpj_mCherry_filtered.loom")

with open("allTFs_mm.txt") as f:
    tfs = set(line.strip() for line in f)

genes = set(adata.var_names)
print(f"TFs in matrix: {len(tfs & genes)} / {len(tfs)}")

print(adata.shape)


type(adata.X)


import numpy as np
adata.X = np.array(adata.X)



