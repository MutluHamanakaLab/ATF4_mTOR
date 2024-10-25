level3.seurat <- readRDS("./prj/fibroblast-lean-seurat-20240308.RDS")
level3.seurat <- FindClusters(level3.seurat, resolution = 0.2)




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
