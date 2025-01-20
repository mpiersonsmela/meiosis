import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from bioinfokit import visuz
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

##### MAKE VOLCANO PLOTS #####

fc_thresh = 1
pval_thresh = 10**-17

degs0 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]


fc_thresh = 4
pval_thresh = 10**-10

degs1 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]


fc_thresh = 6
pval_thresh = 10**-5

degs2 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]


fc_thresh = 8
pval_thresh = 0.05

degs3 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]

print(len(degs0), len(degs1), len(degs2), len(degs3))

degs = pd.concat([degs0, degs1, degs2, degs3])
degs = degs.drop_duplicates(subset="region", keep="first")

region_gene_map = dict(zip(degs["region"], degs["gene_name"]))

visuz.GeneExpression.volcano(df=degs, geneid = "region", lfc=fc_col, pv=pval_col,
                                    dotsize = 3, valpha = 0.5,
                                    color = ("blue", "gray", "red"), xlm = (-9.0, 5, 1), ylm = (0, 55,5), 
                                    sign_line = True,
                                    genenames = region_gene_map, gstyle = 2, 
                                    figtype = "svg",
                                    figname = str(DATA_PATH / f'volcano_{celltype}_labels_only'),
                                    title = f'{celltype}_our_cells_vs_atlas')