library(DESeq2)

# Read  pseudobulk matrix
counts <- read.csv("pseudobulk_export.csv", row.names=1)  # peak x celltype_sample
# Parse out cell_type and sample from row names
metadata <- data.frame(
  row.names = colnames(counts),
  celltype = sub("X.*", "", colnames(counts)),
  sample = sub(".*X", "", colnames(counts))
)

metadata$dataset = ifelse(grepl("24047.05", metadata$sample), "meiotic", "garcia_ATAC")

celltypes <- unique(metadata$celltype)

# Loop over each cell type
for (ct in celltypes) {
  if (ct == "oocyte"){
    next
  }
  # Identify which pseudobulk columns correspond to this cell type
  keep_cols <- rownames(metadata)[metadata$celltype == ct]
  
  # If there aren't at least two columns, skip
  if (length(keep_cols) < 2) {
    message("Skipping cell type ", ct, ": not enough pseudobulk columns.")
    next
  }
  
  # Subset the count matrix and the sample info
  counts_sub <- counts[, keep_cols, drop=FALSE]
  metadata_sub <- metadata[keep_cols, , drop=FALSE]

  design_formula <- ~ dataset
  
  dds_sub <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = metadata_sub,
    design    = design_formula
  )
  
  # Run the DESeq2 pipeline
  dds_sub <- DESeq(dds_sub)
  
  # Extract results comparing "meiotic" vs "garcia_ATAC"
  res <- results(dds_sub, contrast = c("dataset", "meiotic", "garcia_ATAC"))
  
  # Convert results to a data.frame for easier handling
  res_df <- as.data.frame(res)
  
  write.csv(res_df, paste0(ct, "_DESeq2.csv"), sep = ",")
  with(
    res_df,
    plot(
      log2FoldChange,
      -log10(pvalue),
      pch = 20,
      main = paste("Volcano Plot for Cell Type:", ct),
      xlab = "log2FoldChange",
      ylab = "-log10(pvalue)",
      xlim = c(-5, 5)
    )
  )
  
  # Add a horizontal line at FDR=0.05 => -log10(0.05) ~ 1.3
  abline(h = -log10(0.05), col = "red", lty = "dashed")
  
  # Calculate the % of peaks with p < 0.05
  fdr_threshold <- 0.05
  num_sig <- sum(res_df$pvalue < fdr_threshold, na.rm = TRUE)
  total <- sum(!is.na(res_df$pvalue))
  pct <- round((num_sig / total) * 100, 2)
  
  # Add the percentage to the top-left corner of the plot
  text_x <- min(res_df$log2FoldChange, na.rm = TRUE) * 0.8  # a bit left of center
  text_y <- max(-log10(res_df$pvalue), na.rm = TRUE) * 0.9     # near top
  text(
    x = text_x,
    y = text_y,
    labels = paste0(pct, "% of peaks\np-val<0.05"),
    col = "blue",
    pos = 4
  )
}