## Load Necessary Libraries
library(tidyverse)
library(Seurat)
library(decoupleR)
library(ggplot2)

## Read Seurat Object
# Load the Seurat object
level4.seurat_3 <- readRDS("/PATH/TO/FILE/level4.seurat_3.rds")

# Set the cluster identities to the desired order
Idents(level4.seurat_3) <- factor(Idents(level4.seurat_3), levels = c("Fibrotic", "Inflammatory", "Alveolar"))

## Network Analysis
# Retrieve CollecTRI network data for human
net_tri <- get_collectri(organism = 'human', split_complexes = FALSE)

# Extract counts from the RNA assay
mat <- as.matrix(level4.seurat_3@assays$RNA@data)

# Run ULM
acts <- run_ulm(mat = mat, net = net_tri, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)

## Data Preparation
data <- level4.seurat_3

# Extract ULM results and store them in a new assay called tfsulm
data[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change the default assay to tfsulm
DefaultAssay(object = data) <- "tfsulm"

# Scale the data
data <- ScaleData(data)

# Transfer the scaled data to the data slot
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data

## Visualization of Transcription Factor Activities
# Define transcription factors of interest
TF_Interest <- c("SMAD3", "SMAD4", "ATF4", "HIF1A", "YY1")

# Create FeaturePlot for the transcription factors and save the plot
feature_plots <- (FeaturePlot(data, features = TF_Interest) & scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red'))
ggsave("/PATH/TO/OUTPUT/Figure-5F_FeaturePlot.png", feature_plots, width = 7, height = 10.5)

## Dot Plot for Transcription Factor Activities
# Calculate plot width based on the number of transcription factors
W <- length(TF_Interest) * 1.1

# Create DotPlot for the transcription factors and save high and low height versions
scDotplot <- (DotPlot(data, features = TF_Interest) & scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) * 1.4
ggsave("/PATH/TO/OUTPUT/Figure-5F_Dotplot_HIGH.png", scDotplot, width = W, height = H)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave("/PATH/TO/OUTPUT/Figure-5F_Dotplot_LOW.png", scDotplot, width = W, height = H)
