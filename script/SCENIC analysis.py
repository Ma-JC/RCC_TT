 #!/usr/bin/env python
# coding: utf-8
library(Seurat)
library(harmony)

data = readRDS("./data/Mesenchymal_clean.rds")

data$All_third_cluster = factor(data$All_third_cluster,
                                labels = c("CD36+ Pericyte 1","CD36+ Pericyte 2","S100A4+ Pericyte 1","S100A4+ Pericyte 2","IFI+ Pericyte","Cycling Mesenchymal cell","vSMC/Pericyte",
                                                                    "RRAD+ vSMC","ITGA8+ vSMC","FOSB+ vSMC",
                                                                    "ADH1B+ Fibroblast","FBLN5+ Fibroblast","MFAP5+ Fibroblast","CYSLTR2+ Fibroblast","FAP+ Fibroblast",
                                                                    "Mesangial cell"),
                                levels = c("CD36+ Pericyte 1","CD36+ Pericyte 2","S100A4+ vSMC 1","S100A4+ vSMC 2","IFI+ Pericyte","Cycling Mesenchymal cell","vSMC/Pericyte",
                                                                    "RRAD+ vSMC","ITGA8+ vSMC","FOSB+ vSMC",
                                                                    "ADH1B+ Fibroblast","FBLN5+ Fibroblast","MFAP5+ Fibroblast","CYSLTR2+ Fibroblast","FAP+ Fibroblast",
                                                                    "Mesangial cell")
                                )

data$All_third_cluster = as.character(data$All_third_cluster)

if( !dir.exists("SCENIC")){
    dir.create("SCENIC")
}

if( !dir.exists("SCENIC/INPUT")){
    dir.create("SCENIC/INPUT")
}

if( !dir.exists("SCENIC/OUTPUT")){
    dir.create("SCENIC/OUTPUT")
}



Count_Data = as.data.frame(data@assays$RNA@counts)
write.table(Count_Data,file = "./SCENIC/INPUT/Mesenchymal_Count_SCENIC.txt",sep = "\t",quote = F,row.names = T,col.names = T)



Meta_Data = data@meta.data
write.table(Meta_Data,file = "./SCENIC/INPUT/Mesenchymal_Meta_SCENIC_meta.txt",sep = "\t",quote = F,row.names = T,col.names = T)


# ### Run SCEINC 

# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt



sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600) 
# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 32 



exprmat = pd.read_csv("./SCENIC/INPUT/Mesenchymal_Count_SCENIC.txt",sep = "\t",index_col=0) 
exprmat = exprmat.transpose() 
adata = sc.AnnData(exprmat)

cla = pd.read_csv("./SCENIC/INPUT/Mesenchymal_Meta_SCENIC_meta.txt",sep = "\t") 

if not os.path.exists("./SCENIC/OUTPUT/SCENIC_INPUT/"):
    os.makedirs("./SCENIC/OUTPUT/SCENIC_INPUT/")
    
if not os.path.exists("./SCENIC/OUTPUT/RESULT/"):
    os.makedirs("./SCENIC/OUTPUT/RESULT/")



# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt = "./SCENIC/OUTPUT/SCENIC_INPUT/SCENIC_Mesenchymal_unfiltered.loom" # test dataset, n=500 cells

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "./SCENIC/OUTPUT/SCENIC_INPUT/SCENIC_Mesenchymal_filtered_scenic.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = "./SCENIC/OUTPUT/SCENIC_INPUT/SCENIC_Mesenchymal_anndata.h5ad"
f_anndata_vis_path = "./SCENIC/OUTPUT/SCENIC_INPUT/SCENIC_Mesenchymal_anndata_vis.h5ad"
# path to pyscenic output
f_pyscenic_output = "./SCENIC/OUTPUT/RESULT/SCENIC_Mesenchymal_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = './SCENIC/OUTPUT/RESULT/SCENIC_Mesenchymal_scenic_integrated-output.loom'

row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) , 
    "Clusters": np.array(cla["All_third_cluster"]), 
    "Patient": np.array(cla["Patient"]), 
    "thrombus": np.array(cla["newmeta"]),
    "Sample_Tag": np.array(cla["Sample_Tag"]), 
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() , 
}

lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs ) 

# read unfiltered data from a loom file
adata = sc.read_loom( f_loom_path_unfilt ) 

nCountsPerGene = np.sum(adata.X, axis=0) 
nCellsPerGene = np.sum(adata.X>0, axis=0)

# Show info
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

nCells=adata.X.shape[0]

# pySCENIC thresholds
minCountsPerGene=3*.01*nCells 
print("minCountsPerGene: ", minCountsPerGene)

minSamples=.01*nCells 
print("minSamples: ", minSamples)

sc.pp.filter_cells(adata, min_genes=0) 

mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 

adata.obs['n_counts'] = adata.X.sum(axis=1).A1


# In[17]:


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=True) 
x = adata.obs['n_genes'] 
x_lowerbound = 2000
x_upperbound =4000
nbins=100

sns.distplot(x, ax=ax1, norm_hist=True, bins=nbins)
sns.distplot(x, ax=ax2, norm_hist=True, bins=nbins)
sns.distplot(x, ax=ax3, norm_hist=True, bins=nbins)

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


# In[18]:


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=True)

x = adata.obs['percent_mito']
x_lowerbound = [0.0, 0.05]
x_upperbound = [ 0.15, 0.2]
nbins=100

sns.distplot(x, ax=ax1, norm_hist=True, bins=nbins)
sns.distplot(x, ax=ax2, norm_hist=True, bins=int(nbins/(x_lowerbound[1]-x_lowerbound[0])) )
sns.distplot(x, ax=ax3, norm_hist=True, bins=int(nbins/(x_upperbound[1]-x_upperbound[0])) )

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


# In[19]:


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), dpi=150, sharey=False)

sns.distplot( adata.obs['n_genes'], ax=ax1, norm_hist=True, bins=100)
sns.distplot( adata.obs['n_counts'], ax=ax2, norm_hist=True, bins=100)
sns.distplot( adata.obs['percent_mito'], ax=ax3, norm_hist=True, bins=100)

ax1.title.set_text('Number of genes expressed per cell')
ax2.title.set_text('Counts per cell')
ax3.title.set_text('Mitochondrial read fraction per cell')

fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

fig.tight_layout()

fig.savefig('filtering_panel_prefilter.pdf', dpi=600, bbox_inches='tight')


sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_mito')

sc.pp.filter_genes(adata, min_cells=3 )


adata.write( f_anndata_path )

row_attrs = {
    "Gene": np.array(adata.var_names) ,
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "Clusters": np.array(adata.obs["Clusters"]),
    "Patient": np.array(adata.obs["Patient"]),
    "thrombus": np.array(adata.obs["thrombus"]),
    "Sample_Tag": np.array(adata.obs["Sample_Tag"]),
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() , 
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


adata = sc.read_h5ad( f_anndata_path )

# Total-count normalize (library-size correct) to 10,000 reads/cell
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4) # 1e4*gene/total_genes

# log transform the data.
sc.pp.log1p(adata) 

adata.raw = adata

# identify highly variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=4, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

# keep only highly variable genes:
adata = adata[:, adata.var['highly_variable']] 

# scale each gene to unit variance, clip values exceeding SD 10.
sc.pp.scale(adata, max_value=10)

# update the anndata file:
# adata.write( f_anndata_path )

# adata = sc.read_h5ad( f_anndata_path )
# principal component analysis
sc.tl.pca(adata, svd_solver='arpack',n_comps = 30)
sc.pl.pca_variance_ratio(adata, log=True,n_pcs = 30)
# adata.write( f_anndata_path )

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30) 

sc.tl.umap(adata)

sc.pl.umap(adata, color=['Patient'] ) 

sc.external.pp.bbknn(adata, batch_key='Patient',n_pcs = 30,n_trees = 100,trim = 20)

sc.tl.umap(adata)

sc.pl.umap(adata, color=['Patient'] ) 



sc.pl.umap(adata, color=['Clusters'] ) 
sc.pl.umap(adata, color=['Patient'] ) 


# tSNE
tsne = TSNE( n_jobs=32 )
adata.obsm['X_tsne'] = tsne.fit_transform( adata.X )


adata.write( f_anndata_vis_path )


# ## Step 1

f_tfs = "/home/ncpsb/new_mnt/Other/Reference_data/SCENIC/allTFs_hg38.txt"

get_ipython().system(' /home/ncpsb/anaconda3/envs/SCENIC/bin/arboreto_with_multiprocessing.py      {f_loom_path_scenic}      {f_tfs}      --method grnboost2      --output ./SCENIC/OUTPUT/RESULT/adj.tsv      --num_workers 24      --seed 777')



adjacencies = pd.read_csv("./SCENIC/OUTPUT/RESULT/adj.tsv", index_col=False, sep='\t')
adjacencies.loc[adjacencies.loc[:,"target"] == "POSTN"]


# ## Step 2

import glob
# ranking databases
f_db_glob = "/home/ncpsb/new_mnt/Other/Reference_data/SCENIC/hg38*.feather" 
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "/home/ncpsb/new_mnt/Other/Reference_data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

get_ipython().system('/home/ncpsb/anaconda3/envs/SCENIC/bin/pyscenic ctx ./SCENIC/OUTPUT/RESULT/adj.tsv      {f_db_names}      --annotations_fname {f_motif_path}      --expression_mtx_fname {f_loom_path_scenic}      --output ./SCENIC/OUTPUT/RESULT/reg.csv      --mask_dropouts      --num_workers 24')


# ## Step 3

get_ipython().system('/home/ncpsb/anaconda3/envs/SCENIC/bin/pyscenic aucell      {f_loom_path_scenic}      ./SCENIC/OUTPUT/RESULT/reg.csv      --output {f_pyscenic_output}      --num_workers 24')

import json
import zlib
import base64


lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


pd.DataFrame( lf.ra.Regulons, index=lf.ra.Gene).to_csv("/home/ncpsb/mnt/Analysis_for_each_Patient/04. Discovery in Mesenchymal Cell/SCENIC/OUTPUT/RESULT/Regulons.csv")

pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID).to_csv("/home/ncpsb/mnt/Analysis_for_each_Patient/04. Discovery in Mesenchymal Cell/SCENIC/OUTPUT/RESULT/RegulonsAUC.csv")
