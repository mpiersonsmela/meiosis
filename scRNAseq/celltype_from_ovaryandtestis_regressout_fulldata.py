#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import scanpy as sc
import os
import anndata

data_path = '/Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/scRNAseq_2024/fulldata'
atlas_path = 'atlas_testis_and_ovary_regressed.h5ad'

sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=600, figsize=(5,4), format='png')
sc.settings.figdir = 'Analysis_Output/FULLDATA_2024-05-11_ovary_and_testis_annotation_regressed'

#Filter parameters should be the same between both datasets

def filter_preprocessing(adata, gene_col = 'n_genes_by_counts', counts_col = 'total_counts'):
    print("Intial shape:" + str(adata.shape))
    sc.pp.filter_cells(adata, min_genes=1000)
    sc.pp.filter_genes(adata, min_cells=5)
    print("Filtered by min genes and min cells:" + str(adata.shape))
    
    #Check for too many mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter the data
    adata = adata[adata.obs[gene_col] < 14000,:]
    adata = adata[adata.obs[counts_col] < 100000,:]
    adata = adata[adata.obs.pct_counts_mt < 10,:]
    
    if "is_doublet" in adata.obs.columns.values:
        print("before doublet removal, shape:" + str(adata.shape))
        adata = adata[adata.obs.is_doublet == False,:]
        print("removed doublets from annotated data")
        
    print("Final shape:" + str(adata.shape))
    
    return adata

############### Read the atlas data ###############

#Use cached results from "celltype_from_ovaryandtestis_regressout.ipynb"


atlas = sc.read_h5ad(atlas_path)


# In[2]:


sc.pl.umap(atlas, color='celltypes_ordered', save = "_scanpy_ingest_atlas_testisandovary.png")


# In[5]:


############### Read my data ###############

gene_data = pd.read_csv(os.path.join(data_path,'all_genes.csv'))
cell_meta = pd.read_csv(os.path.join(data_path,'cell_metadata.csv.gz'))

# The DGE_filtered folder contains the expression matrix, genes, and files 
adata = sc.read_mtx(os.path.join(data_path,'count_matrix.mtx.gz'))

# add cell meta data to anndata object
adata.obs = cell_meta
adata.obs.set_index('bc_wells', inplace=True)
adata.obs.index.name = None
adata.obs_names_make_unique()

#add sample info
#adata.obs = adata.obs.merge(diff_key, on = 'sample', how = 'left')
adata.obs.index = adata.obs.index.map(str)

# find genes with nan values and filter
gene_data = gene_data[gene_data.gene_name.notnull()]
notNa = gene_data.index
notNa = notNa.to_list()

# remove genes with nan values and assign gene names
adata = adata[:,notNa]
adata.var = gene_data
adata.var.set_index('gene_name', inplace=True)
adata.var.index.name = None
adata.var_names_make_unique()

#save adata to load faster next time
#try:
#    adata.write(os.path.join(data_path,"adata_unfiltered.h5ad"))
#except Exception as e:
#    print(e)

#Filter the data
adata = filter_preprocessing(adata, gene_col = 'n_genes_by_counts', counts_col = 'total_counts')

############### Combine datasets ###############

#Subset the datasets to the intersection of the genes [note, could also do high-variance genes here]

var_names = atlas.var_names.intersection(adata.var_names)
atlas = atlas[:, var_names]
adata = adata[:, var_names]
atlas.shape


# In[6]:


#CPM normalize adata
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata #Save raw CPMs for plotting

#CPM normalize atlas (Already done)
#sc.pp.normalize_total(atlas, target_sum=1e6)
#sc.pp.log1p(atlas)
#atlas.raw = atlas #Save raw CPMs for plotting

#Do PCA in preparation for ingest
#sc.pp.pca(atlas)

#Check the output
#sc.pp.neighbors(atlas, n_pcs =50)
#sc.tl.umap(atlas)
#sc.pl.umap(atlas, color='celltype', save = "_scanpy_ingest_atlas.png")

#Do ingest!
sc.tl.ingest(adata, atlas, obs='celltype')

############### Save cell types ###############

#Also load other sample metadata

sample_metadata = pd.read_csv("sample_metadata.csv")

#Add sample metadata
adata.obs = adata.obs.merge(sample_metadata, on = "sample", how = "left")
adata.obs.index = adata.obs.index.astype(str)


# In[7]:


adata.obs.to_csv(os.path.join(sc.settings.figdir,"obs_with_celltypes_testis_and_ovary.csv"))


# In[8]:


gene_list = pd.read_csv("gene_list2.csv") #genes for plotting


# In[9]:


sample_metadata["sample"].tolist()


# In[10]:


############### Output ###############

#Plot the output
sc.pl.umap(adata, color=['celltype'], save='_celltypes_atlas_embedding.png')
sc.pl.umap(adata, color=['sample'], legend_fontsize=8, save='_samples_atlas_embedding.png')
sc.pl.umap(adata, color=['day'], legend_fontsize=8, save='_day_atlas_embedding.png')
sc.pl.umap(adata, color=['cell_line'], legend_fontsize=8, save='_cell_line_atlas_embedding.png')

#Plot by sample
sample_list=sample_metadata["sample"].tolist()

#Subset by sample
adata_by_sample = [None] * len(sample_list)
for i, sample in enumerate(sample_list):
    #Subset adata
    adata_by_sample[i] = adata[adata.obs['sample'] == sample,:]

adata_by_sample[0].shape

for i, sample in enumerate(sample_list):
    sc.pl.umap(adata_by_sample[i], color=['celltype'], title = sample, save='_scanpy_ingest_embedding_'+sample+'.png')
    
#Calculate percentages of cell types by sample
celltypes_by_sample = [None] * len(sample_list)

for i, sample in enumerate(sample_list):
    celltypes_by_sample[i] = adata_by_sample[i].obs.celltype.value_counts(normalize=True)
    #print(sample)
    #print(celltypes_by_sample[i])
    
celltypes = pd.concat(celltypes_by_sample, axis = 1)
celltypes.columns = sample_list
celltypes.to_csv(os.path.join(sc.settings.figdir,"celltypes_by_samples.csv"))

############### Gene plots ###############

gene_list = gene_list["Genes"].tolist() #['TEX15','SPATA22','SYCP3', 'SYCP1','REC8','mGreenLantern','DAZL', 'STRA8', 'HORMAD1', 'HOXB5', 'BCL2', 'ZGLP1', 'DDX4', 'KIT', 'POU5F1', 'NANOG', 'CD38', 'MLH1', 'SYCE1', 'SYCE2', 'SYCE3', 'YTHDC2','MAEL', 'MEIOSIN', 'MSH5', 'PIWIL1', 'SMC1B', 'SPO11', 'PRDM9', 'PRDM1', 'NANOS3', 'ZBTB16', 'pCXLE_WPRE_polyA']

for i in gene_list:
    try:
        sc.pl.umap(adata, color=i, color_map='viridis', legend_fontsize=8, use_raw = True, save='_library_'+i+'.png')
        sc.pl.umap(atlas, color=i, color_map='viridis', legend_fontsize=8, use_raw = True, save='_atlas_'+i+'.png')
    except Exception as e:
        print(e)
        continue


# In[11]:


#Calculate percentages of cell types by sample
celltypes_by_sample = [None] * len(sample_list)

for i, sample in enumerate(sample_list):
    celltypes_by_sample[i] = adata_by_sample[i].obs.celltype.value_counts()#normalize=True)
    #print(sample)
    #print(celltypes_by_sample[i])
    
celltypes = pd.concat(celltypes_by_sample, axis = 1)
celltypes.columns = sample_list
celltypes = celltypes.fillna(0).astype(int)
celltypes.to_csv(os.path.join(sc.settings.figdir,"celltypes_by_samples_unnormed.csv"))


# In[12]:


### Next Steps: ###
############### Rank genes for each cell type ###############


# In[ ]:




