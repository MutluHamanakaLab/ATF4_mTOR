## Load Necessary Libraries
library(edgeR)
library(pheatmap)
library(ggplot2)

## Import Data
counts.in_all_mTOR <- TableOfCounts_mTOR_modified  # TableOfCounts_mTOR_modified from Table_of_Counts_mTOR.R

## Define Metadata
counts.metadata_all_mTOR <- data.frame(
  dataset= colnames(counts.in_all_mTOR),
  Treatment= c("UT", "UT", "UT", "UT",
               "TGFb", "TGFb", "TGFb", "TGFb",
               "Torin1", "Torin1", "Torin1", "Torin1",
               "TGFb+Torin1", "TGFb+Torin1", "TGFb+Torin1", "TGFb+Torin1"
               ),
  stringsAsFactors = FALSE
)

## Define Groups
group_all_mTOR <- counts.metadata_all_mTOR$Treatment

## Create DGEList Object
y_all_mTOR <- DGEList(counts=counts.in_all_mTOR,
                      genes=row.names.data.frame(counts.in_all_mTOR),
                      group=group_all_mTOR)

## Filter Low-Expressed Genes
keep_all_mTOR <- rowSums(cpm(y_all_mTOR)>1) >= 1
table(keep_all_mTOR)
y_all_mTOR <- y_all_mTOR[keep_all_mTOR, , keep.lib.sizes=FALSE]

## Normalize the Data
y_all_mTOR <- calcNormFactors(y_all_mTOR, method = "TMM")
counts_all_mTOR <- as.matrix(y_all_mTOR$counts)
logCPM_all_mTOR <- cpm(counts_all_mTOR, prior.count=1, log=TRUE)

## Calculate Z-Score
ZScore_all_mTOR <- t(scale(t(logCPM_all_mTOR)))
ZScore_all_mTOR <- as.data.frame(ZScore_all_mTOR)

## Define Annotation Colors
ann.colors_all_mTOR <- list(
  Treatment= c(`UT` = "red", `TGFb` = "blue", `Torin1`="yellow", `TGFb+Torin1`="green")
)

## Subset Significant Genes
sig.qlf.tgfb_vs_tgfb_torin <- subset.data.frame(qlf.tgfb_vs_tgfb_torin, abs(logFC) >= 0.5 & PValue < 0.05)
sig_G <- sig.qlf.tgfb_vs_tgfb_torin$Symbol

## Generate Z-Score Matrix for Significant Genes
Sig.zscore_mTOR <- ZScore_all_mTOR[sig_G, ]
Sig.zscore.mat_mTOR <- as.matrix(Sig.zscore_mTOR)
Sig.zscore.mat_mTOR <- Sig.zscore.mat_mTOR[complete.cases(Sig.zscore.mat_mTOR),]

## Prepare Heatmap Annotation
heat.annotation_mTOR <- data.frame(counts.metadata_all_mTOR[,2])
colnames(heat.annotation_mTOR) <- "Treatment"
row.names(heat.annotation_mTOR) <- counts.metadata_all_mTOR[,1]

## Set Column Names for Heatmap
d.colnames_mTOR <- c(counts.metadata_all_mTOR[,1])
colnames(Sig.zscore.mat_mTOR) <- d.colnames_mTOR

## Italicize Row Names
newnames_mTOR_mTOR <- lapply(
  rownames(Sig.zscore.mat_mTOR),
  function(x) bquote(italic(.(x)))
)

## Create Heatmap
mTOR_heatmap <- pheatmap(
  Sig.zscore.mat_mTOR,
  annotation_col = heat.annotation_mTOR,
  cluster_cols = FALSE,
  main="",
  annotation_colors = ann.colors_all_mTOR,
  show_colnames = FALSE,
  labels_row = as.expression(newnames_mTOR_mTOR),
  fontsize_row = 1,
  cutree_rows = 2  # cut off by row to 2 main clusters
)

## Save Heatmap Plot
ggsave("/PATH/TO/OUTPUT/HeatMap_Figure-3A.png", plot = mTOR_heatmap, width = 10, height = 15, limitsize = FALSE)
