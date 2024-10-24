## Load Necessary Libraries
library(dorothea)
library(decoupleR)
library(tidyverse)
library(edgeR)
library(ggplot2)

# --- If not done already
## Load Dorothea Regulons
data("dorothea_hs", package = "dorothea")
regulons <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B", "C"))

## Get PROGENy network for human with top 500 interactions
net <- get_progeny(organism = 'human', top = 500)

# ----
## Generate matrix for VIPER analysis of upregulated cluster
# ttop_progeny_tgfb_torin_matrix, ttop_progeny_tgfb_torin_matrix_G, and cluster_up are from Figure_3E.R
ttop_progeny_mTOR_up_matrix <- as.matrix(ttop_progeny_tgfb_torin_matrix[ttop_progeny_tgfb_torin_matrix_G$Gene %in% cluster_up, "t_value", drop = FALSE])

## Run VIPER algorithm to infer transcription factor activities
tf_activities_stat_up <- dorothea::run_viper(ttop_progeny_mTOR_up_matrix, regulons,
    options =  list(minsize = 5, eset.filter = FALSE, 
    cores = 1, verbose = FALSE, nes = TRUE))

## Convert VIPER results to a DataFrame and rename columns
tf_activities_stat_up_df <- tf_activities_stat_up %>%
    as.data.frame() %>%
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "t_value")

## Get top 5 transcription factors with positive and negative NES
top5_positive_NES <- tf_activities_stat_up_df %>%
    filter(NES > 0) %>%
    top_n(5, wt = NES) %>%
    arrange(desc(NES)) %>%
    mutate(GeneID = factor(GeneID))

top5_negative_NES <- tf_activities_stat_up_df %>%
    filter(NES < 0) %>%
    top_n(5, wt = abs(NES)) %>%
    arrange(NES) %>%
    mutate(GeneID = factor(GeneID))

## Combine the results into one data frame for plotting
tf_activities_stat_up_top5 <- bind_rows(top5_positive_NES, top5_negative_NES) %>%
    arrange(desc(NES))

## Create bar plot for top 5 transcription factors
NES_plot_up_top5 <- ggplot(tf_activities_stat_up_top5,aes(x = reorder(GeneID, NES), y = NES)) + 
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
        mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = 
            element_text(angle = 45, hjust = 1, size =12, face= "bold"),
        axis.text.y = element_text(size =14, face= "bold"),panel.grid.minor = element_blank()) +
    xlab("Transcription Factors")+ylab("Normalized Enrichment Score (NES)")+
    ggtitle("TGFb vs TGFb+Torin1")

# Set the file path and name where you want to save the plot
NES_plot_cluster_up_png<- "PATH/TO/OUTPUT/Figure-3F.png"
# Save the plot using ggsave
ggsave(NES_plot_cluster_up_png , plot = NES_plot_up_top5 , width = 8, height = 6, limitsize = FALSE)
