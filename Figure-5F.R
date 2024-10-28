level4.seurat_3 <- readRDS("/PATH/TO/FILE/level4.seurat_3.rds")

net_tri <- get_collectri(organism='human', split_complexes=FALSE)
# Extract the normalized log-transformed counts
mat <- as.matrix(level4.seurat_3@assays$RNA@data)

# Run ulm
acts <- run_ulm(mat=mat, net=net_tri, .source='source', .target='target',
                .mor='mor', minsize = 5)

data <- level4.seurat_3

# Extract ulm and store it in tfsulm in pbmc
data[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "tfsulm"

# Scale the data
data <- ScaleData(data)
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data

TF_Interest <- c("SMAD3", "SMAD4","ATF4", "HIF1A", "YY1")

feature_plots <- (FeaturePlot(data, features = TF_Interest) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red'))
ggsave("/PATH/TO/OUTPUT/Figure-5FscFeaturePlot_level4_3_TF_act.png", feature_plots, width = 7, height = 10.5)

W <- length(TF_Interest) * 1.1


scDotplot<- (DotPlot(data, features = TF_Interest) &
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  RotatedAxis()
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("./output/2024-10-09/scDotplot_level4_3_TF_act_20241009.png",scDotplot, width = W, height = H)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("./output/2024-10-09/scDotplot_level4_3_TF_act_LOW_20241009.png",scDotplot, width = W, height = H)
