#This script describes the workflow of single-cell subclustering (PT cells) used in Figure 2.
##Initial QC may have been performed before this script. (See Methods)

import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3          
sc.logging.print_header()
from matplotlib import rcParams
import matplotlib.font_manager
rcParams['font.sans-serif']=['Arial']
sc.settings.set_figure_params(dpi=100, facecolor='white',fontsize=12)
from matplotlib import rcParams
rcParams['axes.grid'] = False

adata = sc.read("read_anndata.h5ad")

sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata,n_top_genes=3000)
sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adata)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata,n_pcs=15, n_neighbors=40,metric='cosine')

sc.tl.umap(adata,min_dist=0.2,spread=2,maxiter=1000)

sc.pl.umap(adata, color=['check_gene_of_interest'])

sc.tl.leiden(adata)
sc.pl.umap(adata,color='leiden')
