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
ttop_progeny_tgfb <- topTags(qlf.wt_vs_tgfb, n = 5000) # qlf.wt_vs_tgfb from Figure_1B

## Calculate t-values for GSEA input
ttop_logFC <- ttop_progeny_tgfb$table$logFC
ttop_p_values <- ttop_progeny_tgfb$table$fdr # can use PValue as well
ttop_t_value <- -log10(ttop_p_values) * ttop_logFC
ttop_progeny_tgfb$table <- cbind(ttop_progeny_tgfb$table, t_value = ttop_t_value)

## Convert to data frame and rename columns
ttop_progeny_tgfb <- as.data.frame(ttop_progeny_tgfb)
ttop_progeny_tgfb <- ttop_progeny_tgfb %>% dplyr::rename(ID = Symbol)
rownames(ttop_progeny_tgfb) <- NULL

## Prepare matrix for Dorothea analysis
ttop_progeny_tgfb_df <- as.data.frame(ttop_progeny_tgfb[, c(1, 8)]) %>% column_to_rownames('ID')
ttop_progeny_ut_tgfb_matrix <- as.matrix(ttop_progeny_tgfb_df)

## Subset significant genes
Sig_wt_vs_tgfb <- subset.data.frame(wt_vs_tgfb, abs(logFC) >= 0.5 & fdr < 0.05)
sig_G <- Sig_wt_vs_tgfb$Symbol
ttop_progeny_ut_tgfb_matrix_G <- ttop_progeny_ut_tgfb_matrix
ttop_progeny_ut_tgfb_matrix_G$Gene <- rownames(ttop_progeny_ut_tgfb_matrix)

## Filter matrix for significant genes
ttop_progeny_mtor_DEG_TGFb_matrix <- as.matrix(
  ttop_progeny_ut_tgfb_matrix[ttop_progeny_ut_tgfb_matrix_G$Gene %in% sig_G, "t_value", drop = FALSE]
)

## Run VIPER algorithm to infer transcription factor activities
tf_activities_stat_tgfb <- dorothea::run_viper(
  ttop_progeny_mtor_DEG_TGFb_matrix, regulons,
  options = list(minsize = 5, eset.filter = FALSE, cores = 1, verbose = FALSE, nes = TRUE)
)

## Convert VIPER results to a DataFrame and rename columns
tf_activities_stat_tgfb_df <- tf_activities_stat_tgfb %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t_value")

## Get top 10 transcription factors with positive and negative NES
top10_positive_NES <- tf_activities_stat_tgfb_df %>%
  filter(NES > 0) %>%
  top_n(10, wt = NES) %>%
  arrange(desc(NES)) %>%
  mutate(GeneID = factor(GeneID))

top10_negative_NES <- tf_activities_stat_tgfb_df %>%
  filter(NES < 0) %>%
  top_n(10, wt = abs(NES)) %>%
  arrange(NES) %>%
  mutate(GeneID = factor(GeneID))

## Combine the results into one data frame for plotting
tf_activities_stat_tgfb_top10 <- bind_rows(top10_positive_NES, top10_negative_NES) %>%
  arrange(desc(NES))

## Create bar plot for top 10 transcription factors
NES_plot_tgfb_top10 <- ggplot(tf_activities_stat_tgfb_top10, aes(x = reorder(GeneID, NES), y = NES)) +
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
  ggtitle("UT vs TGFb")

## Set the file path and name to save the plot
NES_plot_tgfb_png <- "/PATH/TO/OUTPUT/Figure-1E.png"

## Save the plot using ggsave
ggsave(NES_plot_tgfb_png, plot = NES_plot_tgfb_top10, width = 8, height = 6, limitsize = FALSE)
