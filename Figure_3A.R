## Load Necessary Libraries
library(edgeR)
library(pheatmap)
library(ggplot2)

## Import Data
counts.in_all_mtor <- TableOfCounts_mtor_modified  # TableOfCounts_mtor_modified from Table_of_Counts_mTOR.R

## Define Metadata
counts.metadata_all_mtor <- data.frame(
  dataset= colnames(counts.in_all_mtor),
  Treatment= c("UT", "UT", "UT", "UT",
               "TGFb", "TGFb", "TGFb", "TGFb",
               "Torin1", "Torin1", "Torin1", "Torin1",
               "TGFb+Torin1", "TGFb+Torin1", "TGFb+Torin1", "TGFb+Torin1"
               ),
  stringsAsFactors = FALSE
)

## Define Groups
group_all_mtor <- counts.metadata_all_mtor$Treatment

## Create DGEList Object
y_all_mtor <- DGEList(counts=counts.in_all_mtor,
                      genes=row.names.data.frame(counts.in_all_mtor),
                      group=group_all_mtor)

## Filter Low-Expressed Genes
keep_all_mtor <- rowSums(cpm(y_all_mtor)>1) >= 1
table(keep_all_mtor)
y_all_mtor <- y_all_mtor[keep_all_mtor, , keep.lib.sizes=FALSE]

## Normalize the Data
y_all_mtor <- calcNormFactors(y_all_mtor, method = "TMM")
counts_all_mtor <- as.matrix(y_all_mtor$counts)
logCPM_all_mtor <- cpm(counts_all_mtor, prior.count=1, log=TRUE)

## Calculate Z-Score
ZScore_all_mtor <- t(scale(t(logCPM_all_mtor)))
ZScore_all_mtor <- as.data.frame(ZScore_all_mtor)

## Define Annotation Colors
ann.colors_all_mtor <- list(
  Treatment= c(`UT` = "red", `TGFb` = "blue", `Torin1`="yellow", `TGFb+Torin1`="green")
)

## Subset Significant Genes
sig.qlf.tgfb_vs_tgfb_torin <- subset.data.frame(qlf.tgfb_vs_tgfb_torin, abs(logFC) >= 0.5 & PValue < 0.05)
sig_G <- sig.qlf.tgfb_vs_tgfb_torin$Symbol

## Generate Z-Score Matrix for Significant Genes
GSEA.zscore_mTOR <- ZScore_all_mtor[sig_G, ]
GSEA.zscore.mat_mTOR <- as.matrix(GSEA.zscore_mTOR)
GSEA.zscore.mat_mTOR <- GSEA.zscore.mat_mTOR[complete.cases(GSEA.zscore.mat_mTOR),]

## Prepare Heatmap Annotation
heat.annotation_mTOR <- data.frame(counts.metadata_all_mtor[,2])
colnames(heat.annotation_mTOR) <- "Treatment"
row.names(heat.annotation_mTOR) <- counts.metadata_all_mtor[,1]

## Set Column Names for Heatmap
d.colnames_mTOR <- c(counts.metadata_all_mtor[,1])
colnames(GSEA.zscore.mat_mTOR) <- d.colnames_mTOR

## Italicize Row Names
newnames_mTOR_mtor <- lapply(
  rownames(GSEA.zscore.mat_mTOR),
  function(x) bquote(italic(.(x)))
)

## Create Heatmap
overlapped_mtor_heatmap <- pheatmap(
  GSEA.zscore.mat_mTOR,
  annotation_col = heat.annotation_mTOR,
  cluster_cols = FALSE,
  main="",
  annotation_colors = ann.colors_all_mtor,
  show_colnames = FALSE,
  labels_row = as.expression(newnames_mTOR_mtor),
  fontsize_row = 1,
  cutree_rows = 2  # cut off by row to 2 main clusters
)

## Save Heatmap Plot
ggsave("/PATH/TO/OUTPUT/HeatMap_Figure-3A.png", plot = overlapped_mtor_heatmap, width = 10, height = 15, limitsize = FALSE)
