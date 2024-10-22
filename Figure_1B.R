## Load Necessary Libraries
library(edgeR)
library(EnhancedVolcano)
library(ggplot2)

## Import Count Data
# Read the filtered counts data
TableOfCounts_mtor <- read.csv("/PATH/TO/FILE/TableOfCounts_mTOR_filtered.csv") # Table_of_Counts_mTOR.R

# Read sample metadata
samples_mtor <- read.table("/PATH/TO/FILE/TGFb_Torin_sample.txt", header = TRUE)

## Define Groups
# Create a factor for the experimental groups
group_mtor <- factor(c(samples_mtor$Group))

## Create DGEList Object
# Construct a DGEList object for differential expression analysis
DEGall_mtor <- DGEList(counts = TableOfCounts_mtor, group = group_mtor)

# Filter low-expressed genes
keep_mtor <- filterByExpr(DEGall_mtor)
DEGall_mtor <- DEGall_mtor[keep_mtor, , keep.lib.sizes = FALSE]

# Normalize the data
DEGall_mtor <- calcNormFactors(DEGall_mtor)

## Design Matrix and Dispersion Estimation
# Create the design matrix for the group comparison
design_mtor <- model.matrix(~0 + group_mtor)
design_mtor

# Estimate dispersions
DEGall_mtor <- estimateDisp(DEGall_mtor, design_mtor)

## Fit the Model and Test
# Fit a negative binomial generalized log-linear model
fit <- glmQLFit(DEGall_mtor, design_mtor)

# Plot Biological Coefficient of Variation (BCV)
plotBCV(DEGall_mtor)

# Perform a quasi-likelihood F-test for differential expression
qlf.wt_vs_tgfb <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0))

# Adjust p-values for multiple testing using the False Discovery Rate (FDR)
qlf.wt_vs_tgfb$table$fdr <- p.adjust(qlf.wt_vs_tgfb$table$PValue, "fdr")

## Prepare Data for Visualization
# Convert results to data frame
wt_vs_tgfb <- as.data.frame(qlf.wt_vs_tgfb)
row.names(wt_vs_tgfb) <- wt_vs_tgfb$Symbol

# Genes to label on the volcano plot
genes_2_label <- c("ADH1B", "NPTX1", "BDKRB2", "TWSG1", "AKR1C1", "CTSK", "FER1L6", 
                   "BCL2A1", "RASL11A", "SERPINE1", "MTHFD2", "ASNS", "PSAT1", 
                   "ACTA2", "NOX4", "ITGA11", "IGFBP3", "IL11")

## Create Volcano Plot
wt_vs_tgfb_volcano <- EnhancedVolcano(
  wt_vs_tgfb,
  lab = rownames(wt_vs_tgfb),
  x = 'logFC',
  y = 'fdr',
  selectLab = genes_2_label,
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 4.0,
  labSize = 12.0,  # Increased font size
  col = c('black', 'blue', 'grey', 'red3'),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  title = 'UT vs TGFb',
  axisLabSize = 20
)

## Save Volcano Plot
wt_vs_tgfb_img <- "/PATH/TO/OUTPUT/Figure-1B.png"
plot_width <- 15
plot_height <- 15

# Save the plot using ggsave
ggsave(wt_vs_tgfb_img, plot = wt_vs_tgfb_volcano, width = plot_width, height = plot_height, limitsize = FALSE)
