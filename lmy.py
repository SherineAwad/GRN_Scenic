import scanpy as sc
import numpy as np
import loompy
import scipy.sparse

# Load AnnData
adata = sc.read("anndata.h5ad")

# Use adata.raw if it exists and is populated, else fallback to adata.X
X = adata.raw.X if adata.raw is not None else adata.X
var_names = adata.raw.var_names if adata.raw is not None else adata.var_names

# Convert to dense if sparse
if scipy.sparse.issparse(X):
    X = X.toarray()

# Gene filter: keep genes expressed in >=3 cells and with total counts >= 10
gene_filter = (X > 0).sum(axis=0) >= 3
gene_filter &= X.sum(axis=0) >= 10

X_filt = X[:, gene_filter]
genes_filt = np.array(var_names)[gene_filter]

# Transpose for loom (genes x cells)
M = X_filt.T

# Cell metrics
nGene = (M > 0).sum(axis=0).astype(int)
nUMI = M.sum(axis=0).astype(int)

# Replace NaNs, infs (if any)
nGene = np.nan_to_num(nGene, nan=0, posinf=0, neginf=0)
nUMI = np.nan_to_num(nUMI, nan=0, posinf=0, neginf=0)

# Create loom
row_attrs = {"Gene": genes_filt}
col_attrs = {
    "CellID": np.array(adata.obs_names),
    "nGene": nGene,
    "nUMI": nUMI
}

loompy.create("Rbpj_mCherry_filtered_scenic.loom", M, row_attrs, col_attrs)
print("✅ Filtered loom created successfully.")

