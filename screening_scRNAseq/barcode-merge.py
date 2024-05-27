import numpy as np
import pandas as pd
import scanpy as sc
import os
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

data_path = '/n/groups/church/mpsmela/scRNAseq_2023/analysis/combined-v3/all-sample/DGE_unfiltered'
bc_data_path = "/n/groups/church/mpsmela/scRNAseq_2023/analysis/bc_capture/capture_combined/all-sample/DGE_unfiltered"

sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=600, figsize=(5,4), format='png')
sc.settings.figdir = 'Analysis_Output/2023-10-06_bcmerge'

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

############### Read my data ###############

gene_data = pd.read_csv(os.path.join(data_path,'all_genes.csv.gz'))
cell_meta = pd.read_csv(os.path.join(data_path,'cell_metadata.csv.gz'))

# The DGE_filtered folder contains the expression matrix, genes, and files 
adata = sc.read_mtx(os.path.join(data_path,'count_matrix.mtx.gz'))

# add cell meta data to anndata object
adata.obs = cell_meta
adata.obs.set_index('bc_wells', inplace=True)
adata.obs.index.name = None
adata.obs_names_make_unique()

#add sample info
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

###Load the .h5ad file saved from the atlas analysis

adata = sc.read_h5ad(os.path.join(data_path,"adata_unfiltered.h5ad"))


############### Read barcode data ###############
# reading in gene and cell data
bc_gene_data = pd.read_csv(os.path.join(bc_data_path,'all_genes.csv'))

bc_cell_meta = pd.read_csv(os.path.join(bc_data_path,'cell_metadata.csv')) 

bc_adata = sc.read_mtx(os.path.join(bc_data_path,'count_matrix.mtx.gz')) #.gz on O2

# add cell meta data to barcode anndata object
bc_adata.obs = bc_cell_meta
bc_adata.obs.set_index('bc_wells', inplace=True)
bc_adata.obs.index.name = None
bc_adata.obs_names_make_unique()

#add sample info
bc_adata.obs.index = bc_adata.obs.index.map(str)

# find genes with nan values and filter
bc_gene_data = bc_gene_data[bc_gene_data.gene_name.notnull()]
notNa = bc_gene_data.index
notNa = notNa.to_list()

# remove genes with nan values and assign gene names
bc_adata = bc_adata[:,notNa]
bc_adata.var = bc_gene_data
bc_adata.var.set_index('gene_name', inplace=True)
bc_adata.var.index.name = None
bc_adata.var_names_make_unique()

############### Isolate barcode reads ###############
TF_list = ["AMP477_HDAC6", "AMP519_RAD51", "AMP529_SOD1", "AMP531_SRSF1", "AMP536_UBB", "BW88_RARG", "CK530_BOLL_1", "CK530_BOLL_2", "CK530_BOLL_3", "CK530_DAZL_1", "CK530_DAZL_2", "CK530_DAZL_3", "CK530_DDX4_1", "CK530_DDX4_2", "CK530_DDX4_3", "CK541_PRDM1_1", "CK541_PRDM1_2", "CK541_PRDM1_3", "CK542_SOX17_1", "CK542_SOX17_2", "CK542_SOX17_3", "CK543_ZGLP1_1", "CK543_ZGLP1_2", "CK543_ZGLP1_3", "CK555_NFYC_1", "CK555_NFYC_2", "CK555_NFYC_3", "CK578_NOBOX_1", "CK578_NOBOX_2", "CK578_NOBOX_3", "CK580_FIGLA_1", "CK580_FIGLA_2", "CK580_FIGLA_3", "CPO-1", "CPO-10", "CPO-11", "CPO-12", "CPO-13", "CPO-14", "CPO-15", "CPO-16", "CPO-17", "CPO-18", "CPO-19", "CPO-2", "CPO-20", "CPO-21", "CPO-22", "CPO-23", "CPO-24", "CPO-25", "CPO-26", "CPO-27", "CPO-28", "CPO-29", "CPO-3", "CPO-30", "CPO-31", "CPO-32", "CPO-33", "CPO-34", "CPO-4", "CPO-5", "CPO-6", "CPO-9", "DNMT1_CRISPRi", "TF101_EZH2", "TF102_KCTD19", "TF103_BRCA1", "TF104_DMRTB1", "TF105_HMGB2", "TF106_RAN", "TF107_HMGB1", "TF108_HESX1", "TF109_RNF138", "TF110_HOXB5", "TF111_HOXA5", "TF112_HOXA10", "TF113_HMGB3", "TF114_DMC1", "TF115_TERF1", "TF116_ZCWPW1", "TF117_RBM46", "TF118_RFX5", "TF119_CTCFL", "TF120_MYBL1", "TF121_DMRTC2", "TF122_ANHX", "TF123_RB1", "TF124_MEIOSIN", "TF125_SAP30", "TF126_HDAC2", "TF127_SMC3", "TF128_KDM5B", "TF129_SYMPK", "TF130_MSX1", "TF131_YTHDC2", "TF132_NELFE", "TF133_ZNF541", "TF134_NFYB", "TF135_MEIOC", "TF136_ESX1", "TF137_STRA8", "TF138_PRDM9", "TF144_DMRT1_1", "TF144_DMRT1_2", "TF145_TET1", "TF147_TET3", "TF148_PAX5", "TF149_PAX6a", "TF150_PAX6b", "TF151_DPPA3", "TF153_SMAD1-active", "TF154_SMAD9-active", "TF25_ELK1_1", "TF25_ELK1_2", "TF34_TAF4B_1", "TF34_TAF4B_2"]

#Filter out everything that isn't a TF barcode from the barcode capture data
bc_TFs_adata = bc_adata[:,TF_list]

#Remove any cells found in the TF barcode data that aren't in the main data
missing_cells = bc_TFs_adata.obs.index.difference(adata.obs.index)
good_cells = bc_TFs_adata.obs.index.difference(missing_cells)
bc_TFs_adata_filtered = bc_TFs_adata[good_cells,:]

#Convert to dataframe and save
tfcounts_df = pd.DataFrame.sparse.from_spmatrix(bc_TFs_adata_filtered.X, index=bc_TFs_adata_filtered.obs.index, columns = bc_TFs_adata_filtered.var.index)
tfcounts_df.to_csv(os.path.join(sc.settings.figdir,"captured_barcodes.csv"))

#Filter out everything that isn't a TF barcode from the main data
bc_main_adata = adata[:,TF_list]

#Convert to dataframe and save
tfcounts_main_df = pd.DataFrame.sparse.from_spmatrix(bc_main_adata.X, index=bc_main_adata.obs.index, columns = bc_main_adata.var.index)
tfcounts_main_df.to_csv(os.path.join(sc.settings.figdir,"main_barcodes.csv"))

############### Combine barcode reads (capture + main) ###############

#Add the two counts together
merged_barcodes_df = tfcounts_main_df.sparse.to_dense().add(tfcounts_df.sparse.to_dense(), fill_value = 0)

### Read from .csv backup

merged_barcodes_df = pd.read_csv(os.path.join(sc.settings.figdir,"merged_barcodes.csv"), index_col = 0)

#filter to cells in main data and convert to anndata
filtered_barcodes_df = merged_barcodes_df.loc[adata.obs.index] 
new_barcode_adata = anndata.AnnData(filtered_barcodes_df)

#Drop the TFs from the main adata
not_TFs = pd.Index.difference(adata.var.index,TF_list)
dropped_TFs = adata[:,not_TFs]

#Concatenate the merged barcode adata with the dropped data
merged_adata = anndata.concat([dropped_TFs, new_barcode_adata], axis = 1, join = "outer")

#Overwrite the old adata
adata = merged_adata

#Add TF info as observations
adata.obs[TF_list] = merged_barcodes_df
adata.obs = adata.obs.copy()

############### Save adata to load easily in the future ###############

try:
    adata.write(os.path.join(data_path,"adata_with_barcodes.h5ad"))
except Exception as e: #sometimes writing adata can be weird
    print(e)
    with open("adata_with_barcodes.pickle", 'wb') as f:
        pickle.dump(adata,f)

############### Filtering and basic analysis ###############

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
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save='_leiden_regress.png', palette = sns.color_palette('colorblind'))

#add back sample data (we lost it somehow along the way?)
adata.obs['sample'] = cell_meta['sample']
sc.pl.umap(adata, color=['sample'], legend_fontsize=8, save='_sample.png')
sample_list=["s01_ovaroid","s02_PGP1_topTF","s03_F3_topTF","s04_PGP1_topTF_pool_bulk","s05_PGP1_topTF_pool_DDX4","s06_PGP1_topTF_pool_REC8","s07_F3_topTF_pool_bulk","s08_F3_topTF_pool_DDX4","s09_F3_topTF_pool_SYCP3","s10_PGCLC_mCherry","s11_PGCLC_DNMTi","s12_PGP1_pool_bulk","s13_PGP1_pool_DDX4","s14_PGP1_pool_REC8","s15_F3_pool_bulk","s16_F3_pool_DDX4","s17_F3_pool_SYCP3","s18_noTF","s19_PGCLC_AM580","s20_PGCLC_TFs","s21_PGCLC_TFs_AM580"]

############### Save adata to load easily in the future ###############

try:
    adata.write(os.path.join(data_path,"adata_with_barcodes_regressed.h5ad"))
except Exception as e: #sometimes writing adata can be weird
    print(e)
    with open("adata_with_barcodes_regressed.pickle", 'wb') as f:
        pickle.dump(adata,f)

############### Filtering and basic analysis ###############

#Subset by sample
adata_by_sample = [None] * len(sample_list)

for i, sample in enumerate(sample_list):
    #Subset adata
    adata_by_sample[i] = adata[adata.obs['sample'] == sample,:]

#Plot by sample
for i, sample in enumerate(sample_list):
    sc.pl.umap(adata_by_sample[i], color=['leiden'], title = sample, save='_leiden_'+sample+'.png')
    #plt.title(sample)
    
    
#Gene plots
gene_list = ['TEX15','SPATA22','SYCP3', 'SYCP1','REC8','mGreenLantern','DAZL', 'STRA8', 'HORMAD1', 'HOXB5', 'BCL2', 'ZGLP1', 'DDX4', 'KIT', 'POU5F1', 'NANOG', 'CD38', 'MLH1', 'SYCE1', 'SYCE2', 'SYCE3', 'YTHDC2','MAEL', 'MEIOSIN', 'MSH5', 'PIWIL1', 'SMC1B', 'SPO11', 'PRDM9', 'PRDM1', 'NANOS3', 'ZBTB16', 'pCXLE_WPRE_polyA']

for i in gene_list:
    try:
        sc.pl.umap(adata, color=i, color_map='viridis', legend_fontsize=8, use_raw = True, save='_'+i+'.png')
        
    except Exception as e:
        print(e)
        continue
        
############### Look at marker genes in more detail ###############
        
gene_list2 = gene_list #['DAZL', 'DDX4', 'POU5F1','SOX17','PRDM1','NANOS3','MAEL','DPPA3', 'SPO11', 'SYCP3', 'REC8']

cutoff = 0

for gene in gene_list2:
    try:
        adata.obs[gene + '+'] = (adata.raw[:,gene].X.todense() > cutoff)
        adata.obs[gene + '+'] = adata.obs[gene + '+'].astype('category')
    except Exception as e:
        print(e)
        print(gene + ' not found')
        
sc.pl.umap(adata, color='SPO11+', color_map='viridis', legend_fontsize=8, use_raw = True, save='_SPO11+.png')

############### Barcode analysis ###############
#Need to filter out samples that don't have TFs, in order to avoid bias

samples_with_barcodes = ["s04_PGP1_topTF_pool_bulk","s05_PGP1_topTF_pool_DDX4","s06_PGP1_topTF_pool_REC8","s07_F3_topTF_pool_bulk","s08_F3_topTF_pool_DDX4","s09_F3_topTF_pool_SYCP3","s12_PGP1_pool_bulk","s13_PGP1_pool_DDX4","s14_PGP1_pool_REC8","s15_F3_pool_bulk","s16_F3_pool_DDX4","s17_F3_pool_SYCP3"]
adata_with_barcodes = adata[adata.obs["sample"].isin(samples_with_barcodes),:]
adata_with_barcodes.obs["sample"]


for gene in gene_list2:
    adata[adata.obs[gene+'+']==True, :].obs.to_csv(gene+"+_obs.csv")
    
    sc.tl.rank_genes_groups(adata_with_barcodes, gene + '+', method='t-test', key_added = 'rank_genes_groups_'+gene)

    # The head function returns the top n genes per cluster
    top_markers = pd.DataFrame(adata_with_barcodes.uns['rank_genes_groups_'+ gene]['names']).head(50)
    #print(top_markers)
    top_markers.to_csv(os.path.join(sc.settings.figdir,gene+"_markers.csv"))


    #Save the gene ranking for each cluster
    names_df = pd.DataFrame(adata_with_barcodes.uns['rank_genes_groups_'+gene]['names'])

    scores_df = pd.DataFrame(adata_with_barcodes.uns['rank_genes_groups_'+gene]['scores'])

    pvals_adj_df = pd.DataFrame(adata_with_barcodes.uns['rank_genes_groups_'+gene]['pvals_adj'])

    fc_df = pd.DataFrame(adata_with_barcodes.uns['rank_genes_groups_'+gene]['logfoldchanges'])

    for j in ['False', 'True']:
        cluster_df = pd.concat((names_df[j],scores_df[j],fc_df[j],pvals_adj_df[j]),axis = 1, keys = ('Gene name', 'score','logfc','p_adj'))
        cluster_df.to_csv(os.path.join(sc.settings.figdir,gene+"_ranking_"+j+".csv"))
