## Load Necessary Libraries
library(edgeR)
library(Seurat)

## Read Seurat Object
# Load the Seurat object
level4.seurat_3 <- readRDS("/PATH/TO/OUTPUT/level4.seurat_3.rds") #level4.seurat_3.rds from Figure-5B.R

## Aggregate Expression Data
# Create a new metadata column that combines cell type and sample name
level4.seurat_3$celltype_x_sample <- paste(Idents(level4.seurat_3), level4.seurat_3$Sample_Name, sep = "_")

# Aggregate expression data by the new metadata column
aggregated_data <- Seurat::AggregateExpression(level4.seurat_3, group.by = "celltype_x_sample", return.seurat = TRUE)

# Extract aggregated counts matrix
counts_matrix <- aggregated_data@assays$RNA$counts

## Differential Expression Analysis
# Create a factor representing the groups (Fibrotic and Alveolar)
group <- factor(gsub("-.*", "", colnames(counts_matrix)))
levels(group)

# Filter for Fibrotic and Alveolar groups
filtered_counts <- counts_matrix[, group %in% c("Fibrotic", "Alveolar")]
filtered_group <- factor(group[group %in% c("Fibrotic", "Alveolar")])

# Create DGEList object
dge <- DGEList(counts = filtered_counts, group = filtered_group)

# Normalize data
dge <- calcNormFactors(dge)

# Create design matrix
design <- model.matrix(~0 + filtered_group)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit model and perform differential expression analysis
fit <- glmQLFit(dge, design)
qlf.Alveolar_vs_Fibrotic <- glmQLFTest(fit, contrast = c(-1, 1))

# Adjust p-values for multiple testing using FDR
qlf.Alveolar_vs_Fibrotic$table$fdr <- p.adjust(qlf.Alveolar_vs_Fibrotic$table$PValue, "fdr")
Alveolar_vs_Fibrotic <- as.data.frame(qlf.Alveolar_vs_Fibrotic)

# Subset significant results
Sig_Alveolar_vs_Fibrotic <- subset.data.frame(Alveolar_vs_Fibrotic, fdr < 0.05)

## Prepare Gene List for GSEA
# Extract logFC values and name them by rownames
original_gene_list_Alveolar_vs_Fibrotic <- Sig_Alveolar_vs_Fibrotic$logFC
names(original_gene_list_Alveolar_vs_Fibrotic) <- rownames(Sig_Alveolar_vs_Fibrotic)

# Remove NA values and sort the gene list in decreasing order
gene_list_Alveolar_vs_Fibrotic <- na.omit(original_gene_list_Alveolar_vs_Fibrotic)
gene_list_Alveolar_vs_Fibrotic <- sort(gene_list_Alveolar_vs_Fibrotic, decreasing = TRUE)

# Convert the gene list to a data frame
gene_list_Alveolar_vs_Fibrotic_df <- as.data.frame(gene_list_Alveolar_vs_Fibrotic)

# Save the gene list to a file for GSEA
write.table(gene_list_Alveolar_vs_Fibrotic_df, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE, "/PATH/TO/OUTPUT/Figure-5A_GSEA.rnk")
