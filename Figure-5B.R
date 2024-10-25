## Load Necessary Libraries
library(Seurat)
library(ggplot2)

## Read Seurat Object
# Load the Seurat object
level3.seurat <- readRDS("PATH/TO/FILE/fibroblast-lean-seurat-20240308.RDS")

## Clustering Analysis
# Perform clustering at a resolution of 0.2
level3.seurat <- FindClusters(level3.seurat, resolution = 0.2)
# Cluster sizes: 0    1    2    3    4    5    6
#               1817 1113  375  299  260  252  195

## Subset and Preprocess Data
# Keep only cluster 0
level4.seurat <- subset(level3.seurat, idents = c(0))

# Normalize data and regress out mitochondrial content
level4.seurat <- SCTransform(level4.seurat, vars.to.regress = "mitoRatio")

# Perform PCA
level4.seurat <- RunPCA(level4.seurat)

## PCA Analysis
# Calculate percent of variation associated with each PC
pct <- level4.seurat[["pca"]]@stdev / sum(level4.seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC is less than 5%
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
co2

# Minimum of the two calculations
pcs <- min(co1, co2)
pcs  # 12

## Dimensionality Reduction and Clustering
# Run UMAP using the determined number of PCs
level4.seurat <- RunUMAP(level4.seurat, dims = 1:pcs)

# Find neighbors using the same PCs
level4.seurat <- FindNeighbors(level4.seurat, dims = 1:pcs)

# Perform clustering at a higher resolution of 3
level4.seurat <- FindClusters(level4.seurat, resolution = 3)
table(level4.seurat@meta.data$seurat_clusters)
# Cluster sizes: 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
#              148 125 123 120 115 100  98  94  92  81  77  74  73  68  66  62  51  48  47  41  35  33  18  14  14

## Cluster Annotation
# Define cluster annotations
annotation <- c(
  "Inflammatory", "Alveolar", "Inflammatory", "Inflammatory", "Inflammatory", "Alveolar",
  "Fibrotic", "Alveolar", "Inflammatory", "Fibrotic", "Alveolar", "Fibrotic",
  "Fibrotic", "Inflammatory", "Inflammatory", "Inflammatory", "Fibrotic", "Inflammatory",
  "Inflammatory", "Inflammatory", "Inflammatory", "Inflammatory", "Inflammatory", "Inflammatory",
  "Inflammatory"
)

# Assign names to the annotations
names(annotation) <- levels(level4.seurat)

# Rename cluster identities using the annotations
level4.seurat_3 <- RenameIdents(level4.seurat, annotation)

# Save the annotated Seurat object
saveRDS(level4.seurat_3, "/PATH/TO/OUTPUT/level4.seurat_3.rds")

## Plotting
# Plot with legends
scDimplot <- DimPlot(level4.seurat_3, label = TRUE, label.size = 6)

# Clean plot with no legend version
scDimplot <- DimPlot(level4.seurat_3) + NoLegend()

# Save the plot
ggsave("/PATH/TO/OUTPUT/Figure-5B.png", scDimplot, width = 10, height = 10, units = "cm")
