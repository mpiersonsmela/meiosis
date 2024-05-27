import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
import seaborn as sns
import anndata

data_path = 'analysis/combined/all-sample/DGE_unfiltered'
sample_metadata = pd.read_csv("metadata/sample_metadata.csv")
gene_list = pd.read_csv("metadata/gene_list2.csv") #genes for plotting

sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=600, figsize=(5,4), format='png')
sc.settings.figdir = '2024-05-12_Analysis_Output_Fulldata'

#Filter parameters should be the same between both datasets

def filter_preprocessing(adata, gene_col = 'n_genes_by_counts', counts_col = 'total_counts'):

    #Check for too many mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pl.violin(adata, ['n_genes_by_counts'], save='_n_genes_prefilter', jitter=0.4)
    sc.pl.violin(adata, ['total_counts'], save='_total_counts_prefilter', jitter=0.4)
    sc.pl.violin(adata, ['pct_counts_mt'], save='_mito_pct_prefilter', jitter=0.4)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_gene_vs_transcript_counts_prefilter')

    print("Intial shape:" + str(adata.shape))
    sc.pp.filter_cells(adata, min_genes=1000)
    sc.pp.filter_genes(adata, min_cells=5)
    print("Filtered by min genes and min cells:" + str(adata.shape))
    
    # Filter the data
    adata = adata[adata.obs[gene_col] < 14000,:]
    adata = adata[adata.obs[counts_col] < 100000,:]
    adata = adata[adata.obs.pct_counts_mt < 10,:]
    
    if "is_doublet" in adata.obs.columns.values:
        print("before doublet removal, shape:" + str(adata.shape))
        adata = adata[adata.obs.is_doublet == False,:]
        print("removed doublets from annotated data")
        
    print("Final shape:" + str(adata.shape))

    sc.pl.violin(adata, ['n_genes_by_counts'], save='_n_genes_postfilter', jitter=0.4)
    sc.pl.violin(adata, ['total_counts'], save='_total_counts_postfilter', jitter=0.4)
    sc.pl.violin(adata, ['pct_counts_mt'], save='_mito_pct_postfilter', jitter=0.4)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_gene_vs_transcript_counts_postfilter')

    
    return adata

############### Read my data ###############

gene_data = pd.read_csv(os.path.join(data_path,'all_genes.csv'))
cell_meta = pd.read_csv(os.path.join(data_path,'cell_metadata.csv'))

# The DGE_filtered folder contains the expression matrix, genes, and files 
adata = sc.read_mtx(os.path.join(data_path,'count_matrix.mtx'))

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

#Filter the data
adata = filter_preprocessing(adata, gene_col = 'n_genes_by_counts', counts_col = 'total_counts')

#CPM normalize adata
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata #Save raw CPMs for plotting

#Regress out transcript counts
sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)

#PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='') # scanpy generates the filename automatically

#Clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)

#save adata to load faster next time
try:
    adata.write(os.path.join(data_path,"adata_pre_leiden.h5ad"))
except Exception as e:
    print(e)

#Perform Leiden clustering
sc.tl.leiden(adata, resolution = 0.67)

#Add sample metadata
adata.obs = adata.obs.merge(sample_metadata, on = "sample", how = "left")
adata.obs.index = adata.obs.index.astype(str)

#Plot by sample groups
sc.pl.umap(adata, color=['sample'], legend_fontsize=8, save='_by_sample.png')
sc.pl.umap(adata, color=['day'], legend_fontsize=8, save='_by_day.png')
sc.pl.umap(adata, color=['cell_line'], legend_fontsize=8, save='_by_cells.png')
sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save='_by_leiden.png')#, palette = sns.color_palette('colorblind'))

#Rank markers by cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# The head function returns the top n genes per cluster
top_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(15)
print(top_markers)
top_markers.to_csv(os.path.join(sc.settings.figdir,"2024-05-12_clusters_top_markers_regress.csv"))

#Save the gene ranking for each cluster
names_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])

scores_df = pd.DataFrame(adata.uns['rank_genes_groups']['scores'])

pvals_adj_df = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])

fc_df = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])

for i in range(len(fc_df.columns)):
    j = str(i)
    cluster_df = pd.concat((names_df[j],scores_df[j],fc_df[j],pvals_adj_df[j]),axis = 1, keys = ('Gene name', 'score','logfc','p_adj'))
    cluster_df.to_csv(os.path.join(sc.settings.figdir,"2024-05-12_cluster_"+j+"_markers.csv"))

adata.obs.index = adata.obs.index.astype(str) #squash bugs

#Plot gene expression UMAPs
for row in gene_list.iterrows():
    print(row)
    try:
        sc.pl.umap(adata, color=row[1]["Genes"], color_map='viridis', legend_fontsize=8, use_raw = True, save='_'+row[1]["Category"]+'_'+row[1]["Genes"]+'.svg')
    except Exception as e:
        print(e)
        continue

### Generate groups for plotting in the violin plot

s = gene_list['Category']
# Find the first and last index of each unique value
groups = s.groupby(s).groups
#result = [(min(indices), max(indices)) for indices in groups.values()]
positions = sorted([(min(indices), max(indices)) for indices in groups.values()])

#Generate labels
labels = s.drop_duplicates(keep='first').sort_index().tolist()

#Change "day" to categorical
adata.obs["day"]=adata.obs["day"].astype("category")

#Violin plot: by day
sc.pl.stacked_violin(adata, gene_list["Genes"], groupby='day', bw = 0.5, figsize = (12,6), #num_categories=16, 
                     save = 'genes_by_day_violin.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap="viridis")

#Violin plot: by cluster, with dendrogram
sc.pl.stacked_violin(adata, gene_list["Genes"], groupby='leiden', dendrogram=True,bw = 0.5, figsize = (12,6), #num_categories=16, 
                     save = 'genes_by_cluster_violin_dendrogram.svg', var_group_positions = positions,
                     var_group_labels = labels, cmap="viridis")

#Force-directed graph
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='leiden', save = 'FDG_by_cluster.svg')
sc.pl.draw_graph(adata, color='day', save = 'FDG_by_day.svg')

#Dendrogram
sc.tl.dendrogram(adata, groupby="leiden")
sc.pl.dendrogram(adata, groupby="leiden", save = 'dendrogram_by_cluster.svg')

#save adata to load faster next time
try:
    adata.write(os.path.join(data_path,"adata_final_2024-05-12.h5ad"))
except Exception as e:
    print(e)
