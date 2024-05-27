#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
import seaborn as sns
#import igraph

data_path = '/Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/scRNAseq_2024/fulldata'
adata = sc.read_h5ad(os.path.join(data_path,"adata_final_2024-05-12.h5ad"))
sample_metadata = pd.read_csv("sample_metadata.csv")
gene_list = pd.read_csv("gene_list2.csv") #genes for plotting

sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=600, figsize=(5,4), format='png')
sc.settings.figdir = '2024-05-20_Analysis_Output_Fulldata'
adata


# In[2]:


obs_with_celltypes = pd.read_csv("/Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/scRNAseq_2024/code/Analysis_Output/FULLDATA_2024-05-11_ovary_and_testis_annotation_regressed/obs_with_celltypes_testis_and_ovary.csv", index_col=0)
obs_with_celltypes


# In[3]:


obs_with_celltypes['celltype']


# In[4]:


#add the cell type info
adata.obs['celltype'] = obs_with_celltypes.loc[:,'celltype'].values


# In[5]:


adata.obs


# In[6]:


#celltype_order_for_plot = ["Pluripotent Cells / PGCs", "Gonadal Germ Cells", "Mitotic Germ Cells", "Pre-Meiotic (STRA8+) Oogonia", "Meiotic Oogonia", "Pre-Oocytes", "Oocytes", "Pre-Spermatogonia", "Differentiating Spermatogonia", "Spermatogonial Stem Cells", "Early Meiotic Spermatocytes", "Late Meiotic Spermatocytes", "Round Spermatids", "Elongating Spermatids", "Sperm"]
celltype_order_for_plot = ["Pluripotent Cells / PGCs", "Gonadal Germ Cells", "Mitotic Germ Cells", "Pre-Meiotic (STRA8+) Oogonia", "Meiotic Oogonia", "Pre-Oocytes", "Oocytes", "Pre-Spermatogonia", "Spermatogonial Stem Cells", "Differentiating Spermatogonia", "Early Meiotic Spermatocytes", "Late Meiotic Spermatocytes", "Round Spermatids", "Elongating Spermatids", "Sperm"]

#Rename the cell types:
rename_dict = {"Primordial Germ Cells":"Pluripotent Cells / PGCs", "Germ Cells":"Gonadal Germ Cells", "STRA8+ Oogonia":"Pre-Meiotic (STRA8+) Oogonia"}
adata.obs['celltype'] = adata.obs['celltype'].replace(rename_dict)

#Order the cell type colum according to our order
adata.obs["celltype"] = pd.Categorical(
    values=adata.obs["celltype"], categories=celltype_order_for_plot, ordered=True
)

adata.obs


# In[7]:


### Make plots ###

#Plot by sample groups
#sc.pl.umap(adata, color=['sample'], legend_fontsize=8, save='_by_sample.png')
#sc.pl.umap(adata, color=['day'], legend_fontsize=8, save='_by_day.png') #note: don't use this, use the earlier one with viridis
sc.pl.umap(adata, color=['cell_line'], legend_fontsize=8, save='_by_cells.png')
sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save='_by_leiden.png')#, palette = sns.color_palette('colorblind'))


# In[8]:


sc.pl.umap(adata, color=['celltype'], legend_fontsize=8, save='_by_celltype.png')


# In[9]:


sc.pl.umap(adata, color=['day'], legend_fontsize=8, save='_by_day_rainbow2.png', palette = sns.color_palette("husl", 16)) #note: don't use this, use the earlier one with viridis



# In[16]:


#test custom colormap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Function to scale colormap values by square root
def scale_colormap(colormap, factor):
    """
    Scales the colormap values
    """
    # Extract colormap colors
    colors = plt.cm.get_cmap(colormap)(np.linspace(0, 1, 256))

    # Apply scaling
    scaled_colors = [colors[0],colors[int(255 - (256*(factor)))],
                     colors[int(255 - (256*(factor)**2))],colors[int(255 - (256*(factor)**3))],
                     colors[int(255 - (256*(factor)**4))],colors[int(255 - (256*(factor)**5))],
                     colors[int(255 - (256*(factor)**6))],colors[255]]

    # Create new colormap
    new_cmap = LinearSegmentedColormap.from_list(f'{colormap}_transformed', colors=scaled_colors)

    return new_cmap

viridis_transformed = scale_colormap('viridis', 2/3)


# In[11]:


gene_list = pd.read_csv("gene_list2.csv") #genes for plotting
s = gene_list['Category']
# Find the first and last index of each unique value
groups = s.groupby(s).groups
#result = [(min(indices), max(indices)) for indices in groups.values()]
positions = sorted([(min(indices), max(indices)) for indices in groups.values()])

#Generate labels
labels = s.drop_duplicates(keep='first').sort_index().tolist()

#Change "day" to categorical and reverse the order
adata.obs["day"]=pd.Categorical(adata.obs["day"], categories=[15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0], ordered = True)

#Violin plot: by day
sc.pl.stacked_violin(adata, gene_list["Genes"], groupby='day', bw = 0.5, figsize = (12,6), #num_categories=16, 
                     save = 'genes_by_day_violin_2024-05-13.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap="viridis")

#Violin plot: by cluster, with dendrogram
sc.pl.stacked_violin(adata, gene_list["Genes"], groupby='leiden', dendrogram=True,bw = 0.5, figsize = (12,6), #num_categories=16, 
                     save = 'genes_by_cluster_violin_dendrogram_2024-05-13.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap="viridis")


# In[18]:


plt.cm.get_cmap('viridis')(np.linspace(0, 1, 256))


# In[22]:


viridis_transformed(np.linspace(0, 1, 256))


# In[15]:


viridis_transformed = scale_colormap('viridis',2/3)


# In[20]:


sc.pl.stacked_violin(adata, gene_list["Genes"], groupby='day', bw = 0.5, figsize = (12,6), #num_categories=16, 
                     save = 'genes_by_day_violin_2024-05-13_viridis23_bestcolormap.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap=viridis_transformed)


# In[24]:


sc.pl.stacked_violin(adata, gene_list["Genes"], groupby='day', bw = 0.5, figsize = (6,12), #num_categories=16, 
                     save = 'genes_by_day_violin_2024-05-13_viridis23_bestcolormap_vertical.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap=viridis_transformed, swap_axes=True)


# In[21]:


sc.pl.dotplot(adata, gene_list["Genes"], groupby='day', figsize = (15,6), #num_categories=16, 
                     save = 'genes_by_day_violin_2024-05-20_dotplot.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap=viridis_transformed)


# In[25]:


#Change "day" to categorical and change the order back to normal
adata.obs["day"]=pd.Categorical(adata.obs["day"], categories=reversed([15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]), ordered = True)


sc.pl.dotplot(adata, gene_list["Genes"], groupby='day', figsize = (6,18), #num_categories=16, 
                     save = 'genes_by_day_violin_2024-05-20_dotplot_vertical.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap=viridis_transformed, swap_axes=True)


# In[ ]:




