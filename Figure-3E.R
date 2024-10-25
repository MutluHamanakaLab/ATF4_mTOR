## Load Necessary Libraries
library(dorothea)
library(decoupleR)
library(tidyverse)
library(edgeR)
library(ggplot2)

## Load Dorothea Regulons
data("dorothea_hs", package = "dorothea")
regulons <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B", "C"))

## Get PROGENy network for human with top 500 interactions
net <- get_progeny(organism = 'human', top = 500)

## Extract top 5000 significant genes based on differential expression analysis
ttop_progeny_tgfb_torin <- topTags(qlf.tgfb_vs_tgfb_torin, n=5000)  # qlf.tgfb_vs_tgfb_torin from Figure_3A.R

## Calculate t-values for GSEA input
ttop_logFC <- ttop_progeny_tgfb_torin$table$logFC
ttop_p_values <- ttop_progeny_tgfb_torin$table$fdr
ttop_t_value <- -log10(ttop_p_values) * ttop_logFC
ttop_progeny_tgfb_torin <- cbind(ttop_progeny_tgfb_torin$table, t_value = ttop_t_value)

## Convert to data frame and rename columns
ttop_progeny_tgfb_torin <- as.data.frame(ttop_progeny_tgfb_torin)
ttop_progeny_tgfb_torin <- ttop_progeny_tgfb_torin %>% dplyr::rename(ID = Symbol)
rownames(ttop_progeny_tgfb_torin) <- NULL

## Prepare matrix for Dorothea analysis
ttop_progeny_tgfb_torin_df <- as.data.frame(ttop_progeny_tgfb_torin[, c(1, 8)]) %>% column_to_rownames('ID')
ttop_progeny_tgfb_torin_matrix <- as.matrix(ttop_progeny_tgfb_torin_df)

ttop_progeny_tgfb_torin_matrix_G <- ttop_progeny_tgfb_torin_matrix
ttop_progeny_tgfb_torin_matrix_G$Gene <- rownames(ttop_progeny_tgfb_torin_matrix)

## Identify clusters from previous heatmap
mTOR_row_cluster <- data.frame(cluster = cutree(mTOR_heatmap$tree_row, k = 2))  # mTOR_heatmap from Figure_3A.R
Sig.zscore.mat_mTOR_export <- as.data.frame(Sig.zscore.mat_mTOR)
Sig.zscore.mat_mTOR_export$Cluster <- mTOR_row_cluster$cluster

## Separate genes by clusters
mTOR_cluster1 <- Sig.zscore.mat_mTOR_export[Sig.zscore.mat_mTOR_export$Cluster == 1, ]
mTOR_cluster1$Cluster <- NULL
mTOR_cluster2 <- Sig.zscore.mat_mTOR_export[Sig.zscore.mat_mTOR_export$Cluster == 2, ]
mTOR_cluster2$Cluster <- NULL

cluster1_genes <- rownames(mTOR_cluster1)  # Check
cluster2_genes <- rownames(mTOR_cluster2)  # Check

## Define clusters for further analysis
cluster_down <- cluster2_genes
cluster_up <- cluster1_genes  # Use cluster_up in Figure_3E.R later

## Generate matrix for VIPER analysis of downregulated cluster
ttop_progeny_mTOR_down_matrix <- as.matrix(ttop_progeny_tgfb_torin_matrix[ttop_progeny_tgfb_torin_matrix_G$Gene %in% cluster_down, "t_value", drop = FALSE])

## Run VIPER algorithm to infer transcription factor activities
tf_activities_stat_down <- dorothea::run_viper(
  ttop_progeny_mTOR_down_matrix, regulons,
  options = list(minsize = 5, eset.filter = FALSE, cores = 1, verbose = FALSE, nes = TRUE)
)

## Convert VIPER results to a DataFrame and rename columns
tf_activities_stat_down_df <- tf_activities_stat_down %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t_value")

## Get top 5 transcription factors with positive and negative NES
top5_positive_NES <- tf_activities_stat_down_df %>%
  filter(NES > 0) %>%
  top_n(5, wt = NES) %>%
  arrange(desc(NES)) %>%
  mutate(GeneID = factor(GeneID))

top5_negative_NES <- tf_activities_stat_down_df %>%
  filter(NES < 0) %>%
  top_n(5, wt = abs(NES)) %>%
  arrange(NES) %>%
  mutate(GeneID = factor(GeneID))

## Combine the results into one data frame for plotting
tf_activities_stat_down_top5 <- bind_rows(top5_positive_NES, top5_negative_NES) %>%
  arrange(desc(NES))

## Create bar plot for top 5 transcription factors
NES_plot_down_top5 <- ggplot(tf_activities_stat_down_top5, aes(x = reorder(GeneID, NES), y = NES)) +
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face= "bold"),
    axis.text.y = element_text(size = 14, face= "bold"),
    panel.grid.minor = element_blank()
  ) +
  xlab("Transcription Factors") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("TGFb vs TGFb+Torin1")

## Set the file path and name to save the plot
NES_plot_cluster_down_png <- "/PATH/TO/OUTPUT/Figure-3E.png"

## Save the plot using ggsave
ggsave(NES_plot_cluster_down_png, plot = NES_plot_down_top5, width = 8, height = 6, limitsize = FALSE)
