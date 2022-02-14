#This script describes the workflow of single-cell subclustering used in Figure 6, including LoH and DCT/CNT/PC cells.
##Initial QC may have been performed before this script. (See Methods)

##for LoH cells
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

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.calculate_qc_metrics(adata,  percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],jitter=0.4, multi_panel=True)

sc.pp.normalize_total(adata, target_sum=1e3)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

sc.pp.regress_out(adata, ['total_counts'])

sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=40, n_pcs=40,metric='cosine')

sc.tl.umap(adata,min_dist=0.001,maxiter=2000)

sc.pl.umap(adata, color=['check_gene_of_interest'])

sc.tl.leiden(adata)
sc.pl.umap(adata,color='leiden')


##############


#for DCT/CNT/PC cells

adata = sc.read("read_anndata.h5ad")

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.calculate_qc_metrics(adata,  percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],jitter=0.4, multi_panel=True)

sc.pp.normalize_total(adata, target_sum=1e3)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

sc.pp.regress_out(adata, ['total_counts'])

sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=40, n_pcs=20,metric='cosine')

sc.tl.umap(adata,min_dist=0.1,maxiter=2000)

sc.pl.umap(adata, color=['check_gene_of_interest'])

sc.tl.leiden(adata)
sc.pl.umap(adata,color='leiden')
