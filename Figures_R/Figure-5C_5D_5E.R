## Load Necessary Libraries
library(Seurat)
library(ggplot2)

## Read Seurat Object
level4.seurat_3 <- readRDS(file = "/PATH/TO/FILE/level4.seurat_3.rds")

# Set the cluster identities to the desired order
Idents(level4.seurat_3) <- factor(Idents(level4.seurat_3), levels = c("Fibrotic", "Inflammatory", "Alveolar"))

## Define Gene Sets
# Genes related to ATF4 targets
ATF4_target_paper <- c("ATF4", "PHGDH", "SHMT2", "ALDH18A1", "PYCR1", "AARS", "GARS", "WARS", "EIF4EBP1", "DDIT4")

# Genes related to Glycolysis
Glycolysis_paper <- c("PFKL", "ENO1", "LDHA", "SLC16A3", "PGM1")

# Genes related to Oxidative Phosphorylation (Oxphos)
Oxphos_paper <- c("NDUFA3", "NDUFA4", "NDUFAB1", "NDUFB2", "NDUFB7", "SDHB", "UQCRQ", "UQCR10", "COX7B", "COX17", "ATP6V0C", "ATP5MC1")

## Dot Plots for ATF4 Target Genes
# Calculate width based on number of genes
W <- length(ATF4_target_paper) * 0.7

# Create DotPlot for ATF4 target genes and save high and low height versions
scDotplot <- DotPlot(level4.seurat_3, features = ATF4_target_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) * 1.4
ggsave("/PATH/TO/OUTPUT/Figure-5C_DotPlot_HIGH.png", scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave("/PATH/TO/OUTPUT/Figure-5C_DotPlot_LOW.png", scDotplot, height = H, width = W)

## Dot Plots for Glycolysis Genes
# Calculate width based on number of genes
W <- length(Glycolysis_paper) * 1.1

# Create DotPlot for Glycolysis genes and save high and low height versions
scDotplot <- DotPlot(level4.seurat_3, features = Glycolysis_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) * 1.4
ggsave("/PATH/TO/OUTPUT/Figure-5D_DotPlot_High.png", scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave("/PATH/TO/OUTPUT/Figure-5D_DotPlot_LOW.png", scDotplot, height = H, width = W)

## Dot Plots for Oxphos Genes
# Calculate width based on number of genes
W <- length(Oxphos_paper) * 0.6

# Create DotPlot for Oxphos genes and save high and low height versions
scDotplot <- DotPlot(level4.seurat_3, features = Oxphos_paper) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) * 1.4
ggsave("/PATH/TO/OUTPUT/Figure-5E_DotPlot_High.png", scDotplot, height = H, width = W)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave("/PATH/TO/OUTPUT/Figure-5E_DotPlot_LOW.png", scDotplot, height = H, width = W)

## Feature Plots for Genes of Interest
plot_width <- 3.5
plot_height <- 3.5

# Combine all genes of interest
GeneOfInterest_paper <- c(ATF4_target_paper, Glycolysis_paper, Oxphos_paper)
available_features <- rownames(level4.seurat_3)

# Loop through genes and create FeaturePlot for each gene
for (gene in GeneOfInterest_paper) {
  if (gene %in% available_features) {
    plot <- FeaturePlot(level4.seurat_3, features = gene, reduction = "umap", label = FALSE, order = TRUE)
    ggsave(filename = paste0("/PATH/TO/OUTPUT/Figure-5_CDE_FeaturePlot", gene, ".png"), plot = plot, width = plot_width, height = plot_height)
  } else {
    message(paste("The gene", gene, "is missing."))
  }
}
