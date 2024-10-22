# Load Necessary Libraries
library(dplyr)
library(tidyr)
library(clusterProfiler)  

# Filter significant genes from wt_vs_tgfb dataset
Sig_wt_vs_tgfb <- subset.data.frame(wt_vs_tgfb, abs(logFC) >= 0.5 & fdr < 0.05) # wt_vs_tgfb From Figure_1B.R

# Merge significant genes with count table based on gene symbol
DEGmatrix_wt_vs_tgfb <- left_join(TableOfCounts_mtor, Sig_wt_vs_tgfb, by = "Symbol")

# Remove any rows with missing values
DEGmatrix_wt_vs_tgfb <- DEGmatrix_wt_vs_tgfb %>% drop_na()

# Create original gene list with log fold change values
original_gene_list_wt_vs_tgfb <- Sig_wt_vs_tgfb$logFC
names(original_gene_list_wt_vs_tgfb) <- DEGmatrix_wt_vs_tgfb$Symbol

# Remove NA values and sort the gene list in decreasing order
gene_list_wt_vs_tgfb <- na.omit(original_gene_list_wt_vs_tgfb)
gene_list_wt_vs_tgfb <- sort(gene_list_wt_vs_tgfb, decreasing = TRUE)
gene_list_wt_vs_tgfb_df <- as.data.frame(gene_list_wt_vs_tgfb)

# Save the ranked gene list for GSEA for Figure 1C and 1D
write.table(gene_list_wt_vs_tgfb_df, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE, "/PATH/TO/OUTPUT/GSEA_wt_vs_tgfb_logFC0.5.rnk")

# Load organism-specific library
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# Load MSigDB gene sets
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

# Select Hallmark gene sets
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
head(m_t2g)

# Perform Gene Set Enrichment Analysis (GSEA)
HM_Sig_wt_vs_tgfb <- GSEA(gene_list_wt_vs_tgfb, TERM2GENE = m_t2g)

# Visualize GSEA results with a dot plot
gse_HM_Sig_wt_vs_tgfb <- dotplot(HM_Sig_wt_vs_tgfb, showCategory = 15, split = ".sign", font = 10) + 
  facet_grid(. ~ .sign) + 
  ggtitle("UT vs TGFb")

# Save the plot as a PNG file
ggsave("/PATH/TO/OUTPUT/Figure-1C.png", plot = gse_HM_Sig_wt_vs_tgfb, width = 8, height = 8, limitsize = FALSE)
