## Load Necessary Libraries
library(edgeR)
library(EnhancedVolcano)
library(ggplot2)

## Prepare Data for Visualization
# Convert results to data frame
tgfb_vs_tgfb_torin <- as.data.frame(qlf.tgfb_vs_tgfb_torin) # qlf.tgfb_vs_tgfb_torin from Figure_3A.R
row.names(tgfb_vs_tgfb_torin) <- tgfb_vs_tgfb_torin$Symbol

# Genes to label on the volcano plot
genes_2_label <- c("MTHFD2","PSAT1","ASNS","PHGDH","GARS1","SLC6A9","SLC1A4","AARS1","AKNA","RGS9","ACKR4","TKT","CXCL6","PTGS2","ABCA9","NR4A2","PTN")

## Create Volcano Plot
tgfb_vs_tgfb_torin_volcano <- EnhancedVolcano(tgfb_vs_tgfb_torin,
    lab = rownames(tgfb_vs_tgfb_torin),
    x = 'logFC',
    y = 'fdr',
    selectLab = genes_2_label,
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 4.0,
    labSize = 12.0,  # Increased font size
    col = c('black', 'blue', 'grey', 'red3'),
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    title = 'TGFb vs TGFb_Torin1',
    axisLabSize = 20)  

## Save Volcano Plot
tgfb_vs_tgfb_torin_img<- "/PATH/TO/OUTPUT/Figure-3B.png"
plot_width <- 15
plot_height <- 15

# Save the plot using ggsave
ggsave(tgfb_vs_tgfb_torin_img , plot = tgfb_vs_tgfb_torin_volcano , width = plot_width, height = plot_height, limitsize = FALSE)
