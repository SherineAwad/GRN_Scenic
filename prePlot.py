import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt



adata = sc.read_h5ad("adata_prefilter.h5ad")

adata.var_names_make_unique()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=True)

x = adata.obs['n_genes']
x_lowerbound = 1500
x_upperbound = 2000
nbins=100



sns.histplot(x, ax=ax1, stat='density', kde=True, bins=nbins)
sns.histplot(x, ax=ax2, stat='density', kde=True, bins=nbins)
sns.histplot(x, ax=ax3, stat='density', kde=True, bins=nbins)



ax2.set_xlim(0,x_lowerbound)
ax3.set_xlim(x_upperbound, adata.obs['n_genes'].max() )

for ax in (ax1,ax2,ax3):
  ax.set_xlabel('')

ax1.title.set_text('n_genes')
ax2.title.set_text('n_genes, lower bound')
ax3.title.set_text('n_genes, upper bound')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')
fig.text(0.5, 0.0, 'Genes expressed per cell', ha='center', va='center', size='x-large')

fig.tight_layout()


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=True)

x = adata.obs['percent_mito']
x_lowerbound = [0.0, 0.07 ]
x_upperbound = [ 0.10, 0.3 ]
nbins=100



sns.histplot(x, ax=ax1, stat='density', kde=True, bins=nbins)
sns.histplot(x, ax=ax2, stat="density", bins=int(nbins/(x_lowerbound[1]-x_lowerbound[0])), kde=True)
sns.histplot(x, ax=ax3, stat="density", bins=int(nbins/(x_upperbound[1]-x_upperbound[0])), kde=True)



ax2.set_xlim(x_lowerbound[0], x_lowerbound[1])
ax3.set_xlim(x_upperbound[0], x_upperbound[1] )
for ax in (ax1,ax2,ax3): 
  ax.set_xlabel('')

ax1.title.set_text('percent_mito')
ax2.title.set_text('percent_mito, lower bound')
ax3.title.set_text('percent_mito, upper bound')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')
fig.text(0.5, 0.0, 'Mitochondrial read fraction per cell', ha='center', va='center', size='x-large')

fig.tight_layout()


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)


sns.histplot(adata.obs['n_genes'], ax=ax1, stat="density", bins=100, kde=True)
sns.histplot(adata.obs['n_counts'], ax=ax2, stat="density", bins=100, kde=True)
sns.histplot(adata.obs['percent_mito'], ax=ax3, stat="density", bins=100, kde=True)

ax1.title.set_text('Number of genes expressed per cell')
ax2.title.set_text('Counts per cell')
ax3.title.set_text('Mitochondrial read fraction per cell')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

fig.tight_layout()

fig.savefig('filtering_panel_prefilter.pdf', dpi=600, bbox_inches='tight')

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
    jitter=0.4, multi_panel=True,save="violin_preQC.png") 
sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_mito', save="scatter_preQC.png")


