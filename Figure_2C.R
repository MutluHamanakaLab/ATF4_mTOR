## Load Necessary Libraries
library(tidyverse)
library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)

## Filter Significant Genes
# Subset the dataframe to include only significant genes (FDR < 0.05)
Sig_ctrl_tgfb_vs_atf4ko_tgfb <- subset.data.frame(ctrl_tgfb_vs_atf4ko_tgfb, fdr < 0.05)  # ctrl_tgfb_vs_atf4ko_tgfb from Figure_2B.R

## Merge Data
# Merge the filtered significant genes with the Table of Counts
DEGmatrix_ctrl_tgfb_vs_atf4ko_tgfb <- left_join(TableOfCounts_atf4, Sig_ctrl_tgfb_vs_atf4ko_tgfb, by = "Symbol")

# Drop rows with NA values
DEGmatrix_ctrl_tgfb_vs_atf4ko_tgfb <- DEGmatrix_ctrl_tgfb_vs_atf4ko_tgfb %>% drop_na()

## Prepare Gene List
# Extract logFC values and name them by the Symbol column
original_gene_list_ctrl_tgfb_vs_atf4ko_tgfb <- Sig_ctrl_tgfb_vs_atf4ko_tgfb$logFC
names(original_gene_list_ctrl_tgfb_vs_atf4ko_tgfb) <- DEGmatrix_ctrl_tgfb_vs_atf4ko_tgfb$Symbol

# Remove NA values and sort the gene list in decreasing order
gene_list_ctrl_tgfb_vs_atf4ko_tgfb <- na.omit(original_gene_list_ctrl_tgfb_vs_atf4ko_tgfb)
gene_list_ctrl_tgfb_vs_atf4ko_tgfb <- sort(gene_list_ctrl_tgfb_vs_atf4ko_tgfb, decreasing = TRUE)

# Convert the gene list to a data frame
gene_list_ctrl_tgfb_vs_atf4ko_tgfb_df <- as.data.frame(gene_list_ctrl_tgfb_vs_atf4ko_tgfb)

# Save the gene list to a file
write.table(gene_list_ctrl_tgfb_vs_atf4ko_tgfb_df, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE, "/PATH/TO/OUTPUT/GSEA_siControl_TGFb_vs_siATF4_TGFb.rnk")

## Pathway Analysis
# Convert gene symbols to Entrez IDs
react_gene_list_ctrl_tgfb_vs_atf4ko_tgfb <- gene_list_ctrl_tgfb_vs_atf4ko_tgfb
gene_id <- mapIds(org.Hs.eg.db, names(react_gene_list_ctrl_tgfb_vs_atf4ko_tgfb), 'ENTREZID', 'SYMBOL')
names(react_gene_list_ctrl_tgfb_vs_atf4ko_tgfb) <- gene_id

# Perform gene set enrichment analysis (GSEA) using Reactome pathways
ctrl_tgfb_vs_atf4ko_tgfb_React <- gsePathway(
  react_gene_list_ctrl_tgfb_vs_atf4ko_tgfb,
  organism = "human",
  pAdjustMethod = "none",
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Set the analysis results to be human-readable
ctrl_tgfb_vs_atf4ko_tgfb_React <- setReadable(ctrl_tgfb_vs_atf4ko_tgfb_React, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## Plot Results
# Create a dot plot for the GSEA results
gse_RE_ctrl_tgfb_vs_atf4ko_tgfb <- dotplot(
  ctrl_tgfb_vs_atf4ko_tgfb_React,
  showCategory = 13,
  split = ".sign",
  font = 17
) + facet_grid(. ~ .sign) + ggtitle("siControl_TGFb vs siATF4_TGFb")

# Save the plot to a file
ggsave("/PATH/TO/OUTPUT/Figure_2C.png", plot = gse_RE_ctrl_tgfb_vs_atf4ko_tgfb, width = 10, height = 5, limitsize = FALSE)
