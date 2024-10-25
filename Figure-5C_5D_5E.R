level4.seurat_3 <- readRDS(file = "/PATH/TO/INPUT/level4.seurat_3.rds")

ATF4_target_paper <- c("ATF4","PHGDH","SHMT2","ALDH18A1","PYCR1","AARS","GARS","WARS","EIF4EBP1","DDIT4")
Glycolysis_paper <- c("PFKL","ENO1","LDHA","SLC16A3","PGM1")  
Oxphos_paper <- c("NDUFA3","NDUFA4","NDUFAB1","NDUFB2","NDUFB7","SDHB","UQCRQ","UQCR10","COX7B","COX17","ATP6V0C","ATP5MC1")

W <- length(ATF4_target_paper) *0.7

scDotplot <- DotPlot(level4.seurat_3, features = ATF4_target_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("/PATH/TO/OUTPUT/Figure-5C_DotPlot_HIGH.png",scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave("/PATH/TO/OUTPUT/Figure-5C_DotPlot_LOW.png",scDotplot, height = H, width = W)

W <- length(Glycolysis_paper) *1.1

scDotplot <- DotPlot(level4.seurat_3, features =Glycolysis_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("/PATH/TO/OUTPUT/Figure-5D_DotPlot_High.png",scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("/PATH/TO/OUTPUT/Figure-5D_DotPlot_LOW.png",scDotplot, height = H, width = W)

W <- length(Oxphos_paper) *0.6

scDotplot <- DotPlot(level4.seurat_3, features =Oxphos_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("/PATH/TO/OUTPUT/Figure-5E_DotPlot_High.png",scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("/PATH/TO/OUTPUT/Figure-5E_DotPlot_LOW.png",scDotplot, height = H, width = W)

plot_width <- 3.5
plot_height <- 3.5

GeneOfInterest_paper <- c(ATF4_target_paper,Glycolysis_paper,Oxphos_paper)
available_features <- rownames(level4.seurat_3)

for (gene in GeneOfInterest_paper) {
  if (gene %in% available_features) {
    plot <- FeaturePlot(level4.seurat_3, features = gene, reduction = "umap", label = F, order = T)
    ggsave(filename = paste0("/PATH/TO/OUTPUT/Figure-5_FeaturePlot_", gene, ".png"), plot = plot, width = plot_width, height = plot_height)
  } else {
    message(paste("The gene", gene, "is missing."))
  }
}

