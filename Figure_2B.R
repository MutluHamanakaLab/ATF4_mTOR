 ## Load Necessary Libraries
library(edgeR)
library(EnhancedVolcano)
library(ggplot2)

## Prepare Data for Visualization
ctrl_tgfb_vs_atf4ko_tgfb # ctrl_tgfb_vs_atf4ko_tgfb from Figure_2A.R

# Genes to label on the volcano plot
gene_label <- c("PSAT1","PHGDH","ASNS","MTHFD2","SLC6A9","GARS1","SLC7A5","EIF4EBP1","DDIT4","WARS1")

## Create Volcano Plot
ctrl_tgfb_vs_atf4ko_tgfb_volcano  <- EnhancedVolcano(ctrl_tgfb_vs_atf4ko_tgfb,
    lab = rownames(ctrl_tgfb_vs_atf4ko_tgfb),
    x = 'logFC',
    y = 'fdr',
    selectLab = gene_label,
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 4.0,
    labSize = 12.0,
    col=c('black', 'blue', 'grey', 'red3'),
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    title = 'siControl_TGFb vs siATF4_TGFb', 
    axisLabSize = 20)

## Save Volcano Plot
ctrl_tgfb_vs_atf4ko_tgfb_img<- "/PATH/TO/OUTPUT/Figure-2B.png"
plot_width <- 15
plot_height <- 15

ggsave(ctrl_tgfb_vs_atf4ko_tgfb_img , plot = ctrl_tgfb_vs_atf4ko_tgfb_volcano , width = plot_width, height = plot_height, limitsize = FALSE)
