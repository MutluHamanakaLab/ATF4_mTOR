# Load Necessary Libraries
library(dplyr)
library(tidyr)
library(clusterProfiler)

# Filter significant genes from wt_vs_tgfb dataset
Sig_tgfb_vs_tgfb_torin <- subset.data.frame(tgfb_vs_tgfb_torin, fdr<0.05)

# Merge significant genes with count table based on gene symbol
DEGmatrix_tgfb_vs_tgfb_torin<- left_join(TableOfCounts_mtor, Sig_tgfb_vs_tgfb_torin, by = "Symbol")

# Remove any rows with missing values
DEGmatrix_tgfb_vs_tgfb_torin<- DEGmatrix_tgfb_vs_tgfb_torin%>%drop_na()

# Create original gene list with log fold change values
original_gene_list_tgfb_vs_tgfb_torin<- Sig_tgfb_vs_tgfb_torin$logFC
names(original_gene_list_tgfb_vs_tgfb_torin)<-DEGmatrix_tgfb_vs_tgfb_torin$Symbol

# Remove NA values and sort the gene list in decreasing order
gene_list_tgfb_vs_tgfb_torin<- na.omit(original_gene_list_tgfb_vs_tgfb_torin)
gene_list_tgfb_vs_tgfb_torin<-sort(gene_list_tgfb_vs_tgfb_torin, decreasing = TRUE)
gene_list_tgfb_vs_tgfb_torin_df<- as.data.frame(gene_list_tgfb_vs_tgfb_torin)

# Save the ranked gene list for GSEA for Figure 3C and 3D
write.table(gene_list_tgfb_vs_tgfb_torin_df, sep = "\t", quote = FALSE, row.names =TRUE, col.names = FALSE, "/PATH/TO/OUTPUT/GSEA_tgfb_vs_tgfb_torin.rnk")

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
HM_Sig_tgfb_vs_tgfb_torin <- GSEA(gene_list_tgfb_vs_tgfb_torin, TERM2GENE = m_t2g)

# Visualize GSEA results with a dot plot
gse_HM_Sig_tgfb_vs_tgfb_torin<- dotplot(HM_Sig_tgfb_vs_tgfb_torin, showCategory=15, split=".sign", font = 20) + facet_grid(.~.sign)+ggtitle("TGFb vs TGFb+Torin1")

# Save the plot as a PNG file
ggsave("/PATH/TO/OUTPUT/Figure-3C.png", plot = gse_HM_Sig_tgfb_vs_tgfb_torin, width = 12, height = 9, limitsize = FALSE)
