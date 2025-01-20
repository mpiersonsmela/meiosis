import os
import pandas as pd
import requests
from requests.structures import CaseInsensitiveDict
from dotenv import load_dotenv
from pathlib import Path


load_dotenv()

DATA_PATH = Path(os.getenv('OUTPUT_PATH')) / 'garcia_ATAC/volcano_plots'
ANNOTATION_PATH = Path(os.getenv('OUTPUT_PATH')) / 'garcia_ATAC/annotation'

#Read the input data
celltype = "GC"
fc_data = pd.read_csv(DATA_PATH / f"{celltype}_DESeq2.csv", index_col=0)
annotations = pd.read_csv(ANNOTATION_PATH / "strand_aware_annotation.final_annotation.csv", index_col=0)
annotations['name'] = annotations.index + ":" + annotations.peak_start.astype(str) + "-" + annotations.peak_end.astype(str)
annotations_unique = annotations.drop_duplicates(subset='name', keep='first')
annotations_unique.index = annotations_unique.name

fc_data = fc_data.merge(annotations_unique[["gene_name", "annotation_type"]], how="left", left_index=True, right_index=True)
fc_data["region"] = fc_data.index

fc_col = "log2FoldChange" 
pval_col = "padj"
fc_data = fc_data.dropna(subset = [fc_col, pval_col])

#Set up the output directory for GO enrichment
plot_dir = "GO_enrich/"
os.makedirs(plot_dir, exist_ok=True)

##### RUN GO ENRICHMENT #####

#Thresholds for enrichment
fc_thresh = 3
pval_thresh = 0.01
go_thresh = 0.05

#PANTHERDB (version 17) parameters. Note, the current version is 10.5281/zenodo.6399963
url_base = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="
#9606 = human
#GO:0008150 = biological process
#alternatively, GO:0003674 = molecular function
url_tail = "&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR"
#url_tail = "&organism=9606&annotDataSet=GO%3A0003674&enrichmentTestType=FISHER&correction=FDR"

#Request headers
headers = CaseInsensitiveDict()
headers["connection"] = "keep-alive"
headers["keep-alive"] = "timeout=3600, max=100"

fc_col = "log2FoldChange" 
pval_col = "padj"
fc_data = fc_data.dropna(subset = [fc_col, pval_col])

degs_up = fc_data[(fc_data[pval_col] < pval_thresh) & (fc_data[fc_col] > fc_thresh)][['gene_name', fc_col, pval_col]]

degs_up.to_csv(plot_dir + 'up_' + str(fc_thresh) + '.csv')

degs_down = fc_data[(fc_data[pval_col] < pval_thresh) & (fc_data[fc_col] < -fc_thresh)][['gene_name', fc_col, pval_col]]
degs_down = degs_down.dropna(subset=["gene_name"])
degs_down.to_csv(plot_dir + 'down_' + str(fc_thresh) + '.csv')

#Get enrichment data for upregulated genes
gene_list = degs_up['gene_name'].unique().tolist()
if len(gene_list) > 1:
    full_url = url_base + ','.join(gene_list) + url_tail

    print("Requesting enrichment data for " + str(len(gene_list)) + " upregulated targets")
    #print(full_url)

    try:
        r = requests.post(full_url, headers = headers)
    except Exception as e:
        print("Failed: " + full_url)
        print(e)

    if r.status_code == 200:
        response_df = pd.DataFrame(r.json()['results']['result'])
        response_df[['id','label']] = pd.json_normalize(response_df['term'])
        response_df.drop(columns = 'term', inplace = True)
        signif = response_df[response_df['fdr'] < go_thresh]
        signif.to_csv(plot_dir + '_' + str(fc_thresh) + "-fold_upregulated_GOs_biological_process.csv")
        #signif.to_csv(gene + '_' + str(fc_thresh) + "-fold_upregulated_GOs_molecular_function.csv")
    else: print("Failed: " + full_url)
    
#Get enrichment data for downregulated genes
gene_list = degs_down['gene_name'].unique().tolist()
# Identify and print any float values in gene_list
float_genes = [gene for gene in gene_list if isinstance(gene, float)]
if float_genes:
    print("Found float values in gene_list:", float_genes)


if len(gene_list) > 1:
    full_url = url_base + ','.join(gene_list) + url_tail

    print("Requesting enrichment data for " + str(len(gene_list)) + " downregulated targets")

    try:
        r = requests.post(full_url)
    except:
        print("Failed: " + full_url)

    if r.status_code == 200:
        response_df = pd.DataFrame(r.json()['results']['result'])
        response_df[['id','label']] = pd.json_normalize(response_df['term'])
        response_df.drop(columns = 'term', inplace = True)
        signif = response_df[response_df['fdr'] < go_thresh]
        signif.to_csv(plot_dir + '_' + str(fc_thresh) + "-fold_downregulated_GOs_biological_process.csv")
        #signif.to_csv(gene + '_' + str(fc_thresh) + "-fold_downregulated_GOs_molecular_function.csv")
    else: print("Request failed")
