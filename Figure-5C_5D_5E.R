

ATF4_target_paper <- c("ATF4","PHGDH","SHMT2","ALDH18A1","PYCR1","AARS","GARS","WARS","EIF4EBP1","DDIT4")
Glycolysis_paper <- c("PFKL","ENO1","LDHA","SLC16A3","PGM1")  
Oxphos_paper <- c("NDUFA3","NDUFA4","NDUFAB1","NDUFB2","NDUFB7","SDHB","UQCRQ","UQCR10","COX7B","COX17","ATP6V0C","ATP5MC1")

W <- length(ATF4_target_paper) *0.7

scDotplot <- DotPlot(level4.seurat_3, features = ATF4_target_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("./output/2024-10-09/scDotplot_level4_3_atf4_20241009.png",scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave("./output/2024-10-09/scDotplot_level4_3_atf4_LOW_20241009.png",scDotplot, height = H, width = W)

W <- length(Glycolysis_paper) *1.1

scDotplot <- DotPlot(level4.seurat_3, features =Glycolysis_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("./output/2024-10-07/scDotplot_level4_3_glycolysis_20241007.png",scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("./output/2024-10-07/scDotplot_level4_3_glycolysis_LOW_20241007.png",scDotplot, height = H, width = W)

W <- length(Oxphos_paper) *0.6

scDotplot <- DotPlot(level4.seurat_3, features =Oxphos_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("./output/2024-10-09/scDotplot_level4_3_OxPhos_20241009.png",scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("./output/2024-10-09/scDotplot_level4_3_OxPhos_LOW_20241009.png",scDotplot, height = H, width = W)

plot_width <- 3.5
plot_height <- 3.5

GeneOfInterest_paper <- c(ATF4_target_paper,Glycolysis_paper,Oxphos_paper)
available_features <- rownames(level4.seurat_3)

for (gene in GeneOfInterest_paper) {
  if (gene %in% available_features) {
    plot <- FeaturePlot(level4.seurat_3, features = gene, reduction = "umap", label = F, order = T)
    ggsave(filename = paste0("./output/2024-10-24/featureplots_level4/FeaturePlots_level4_", gene, ".png"), plot = plot, width = plot_width, height = plot_height)
  } else {
    message(paste("The gene", gene, "is missing."))
  }
}

