level3.seurat <- readRDS("./prj/fibroblast-lean-seurat-20240308.RDS")
level3.seurat <- FindClusters(level3.seurat, resolution = 0.2)

#   0    1    2    3    4    5    6 
#1817 1113  375  299  260  252  195 
# Keep only 0
level4.seurat <- subset(level3.seurat, idents = c(0))

level4.seurat <- SCTransform(level4.seurat, vars.to.regress = "mitoRatio")
level4.seurat <- RunPCA(level4.seurat)

# Determine percent of variation associated with each PC
pct <- level4.seurat[["pca"]]@stdev / sum(level4.seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs #12

level4.seurat <- RunUMAP(level4.seurat, dims = 1:pcs)
level4.seurat <- FindNeighbors(level4.seurat, dims = 1:pcs)

level4.seurat <- FindClusters(level4.seurat, resolution = 3)
table(level4.seurat@meta.data$seurat_clusters)
# Res = 3
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24 
#148 125 123 120 115 100  98  94  92  81  77  74  73  68  66  62  51  48  47  41  35  33  18  14  14 

annotation <- c("Inflammatory","Alveolar","Inflammatory","Inflammatory","Inflammatory","Alveolar",
                "Fibrotic","Alveolar","Inflammatory","Fibrotic","Alveolar","Fibrotic",
                "Fibrotic","Inflammatory","Inflammatory","Inflammatory","Fibrotic","Inflammatory",
                "Inflammatory","Inflammatory","Inflammatory","Inflammatory","Inflammatory","Inflammatory",
                "Inflammatory")
names(annotation) <- levels(level4.seurat)
level4.seurat_3 <- RenameIdents(level4.seurat, annotation)

saveRDS(level4.seurat_3,"./output/rds/level4.seurat_3.rds")

# plot with legend
scDimplot <- DimPlot(level4.seurat_3, label = T,label.size = 6)

# clean plot with no legend version
scDimplot <- DimPlot(level4.seurat_3)+ NoLegend()

ggsave("./output/2024-10-14/scDimplot_level4_3.png",scDimplot, width = 10, height = 10, units = "cm")
