# Load Necessary Libraries
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(edgeR)
library(sva)
library(ggplot2)

# Initialize variables for human use
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- na.omit(tx2gene)

# Load read abundance data for ATF4 and mTOR from CSV files
read_abundance_atf4 <- read.csv("/PATH/TO/THE/FILE/ATF4KD_abundance_filtered.csv") # From Table_of_Counts_atf4.R
read_abundance_mtor <- read.csv("/PATH/TO/THE/FILE/mTOR_abundance_filtered.csv") # From Table_of_Counts_mTOR.R

# Remove Symbol and entrezgene_id columns to combine read_abundances
read_abundance_atf4$Symbol <- NULL
read_abundance_atf4$entrezgene_id <- NULL

# Merge read_abundance data for ATF4 and mTOR
read_abundance_mtor_atf4_combined <- merge(read_abundance_mtor, read_abundance_atf4, by="ensembl_gene_id", all=TRUE)

# Save combined read_abundance data to a CSV file
write_csv(read_abundance_mtor_atf4_combined, "./PATH/TO/OUTPUT/TableOfCounts_mTOR_ATF4_combined.csv")

# Keep only Symbol column for further analysis
TableOfCounts_combined <- read_abundance_mtor_atf4_combined[-c(1,3)]
TableOfCounts_combined_modified <- TableOfCounts_combined[-c(1)]
rownames(TableOfCounts_combined_modified) <- TableOfCounts_combined$Symbol

# Define groups and batch information for the combined dataset
counts.metadata_DEGall_combined <- data.frame(
  dataset = colnames(TableOfCounts_combined_modified),
  Treatment = c("WT","WT","WT","WT", "WT_TGFb","WT_TGFb","WT_TGFb", "WT_TGFb",
                "WT_Torin1","WT_Torin1","WT_Torin1","WT_Torin1",
                "WT_TGFb_Torin1","WT_TGFb_Torin1","WT_TGFb_Torin1", "WT_TGFb_Torin1",
                "siControl_UT","siControl_UT","siControl_UT",
                "siControl_TGFb","siControl_TGFb","siControl_TGFb",
                "siATF4_UT","siATF4_UT","siATF4_UT",
                "siATF4_TGFb","siATF4_TGFb","siATF4_TGFb"),
  batch = c("mTOR","mTOR","mTOR","mTOR",
            "mTOR","mTOR","mTOR","mTOR",
            "mTOR","mTOR","mTOR","mTOR",
            "mTOR","mTOR","mTOR","mTOR",
            "ATF4","ATF4","ATF4",
            "ATF4","ATF4","ATF4",
            "ATF4","ATF4","ATF4",
            "ATF4","ATF4","ATF4"),
  stringsAsFactors = FALSE
)

group_combined <- as.factor(counts.metadata_DEGall_combined$Treatment)
batch_combined <- as.factor(counts.metadata_DEGall_combined$batch)

# Create a new variable that combines batch and treatment group information
combined_group <- factor(paste(batch_combined, group_combined, sep="_"))

# Create a DGEList object for differential expression analysis
DEGall_combined <- DGEList(counts=TableOfCounts_combined, group=group_combined)

# Filter low-expressed genes
keep_combined <- filterByExpr(DEGall_combined)
DEGall_combined <- DEGall_combined[keep_combined,,keep.lib.sizes=FALSE]

# Normalize the data
DEGall_combined_Norm <- calcNormFactors(DEGall_combined, method = "TMM")
DEGall_combined_Norm_logcpm <- DEGall_combined_Norm
DEGall_combined_Norm_logcpm$counts <- edgeR::cpm(DEGall_combined_Norm$counts, log=TRUE)

# Remove batch effects using ComBat
DEGall_combined_Norm_batch <- DEGall_combined_Norm_logcpm
DEGall_combined_Norm_batch$counts <- ComBat(dat= DEGall_combined_Norm_logcpm, batch=batch_combined, mod=NULL, par.prior=TRUE, prior.plots=TRUE)

# Design matrix for group comparison
design_combined <- model.matrix(~0+combined_group)
DEGall_combined_Norm <- estimateDisp(DEGall_combined_Norm, design_combined)
fit_combined <- glmQLFit(DEGall_combined_Norm, design_combined)

# Batch effect visualization setup
conditions_quad <- c("WT_UT", "WT_TGFb", "WT_Torin1", "WT_TGFb_Torin1")
conditions_trip <- c("siControl_UT", "siControl_TGFb", "siATF4_UT", "siATF4_TGFb")
condition_colors_quad <- c("red", "blue", "green", "purple")
condition_colors_trip <- c("orange", "pink", "brown", "gray")

# Create vectors for pch with separate values for quadruplicates and triplicates
pch_quad <- rep(16, length(conditions_quad))
pch_trip <- rep(17, length(conditions_trip))
all_pch <- c(rep(pch_quad, each = 4), rep(pch_trip, each = 4))

# Plot MDS with batch effect removed
opar <- par(no.readonly = TRUE)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))

MDS_plot_Normed_Batch_Removed <- plotMDS(
  DEGall_combined_Norm_batch$counts,
  labels = NULL,
  col = c(rep(condition_colors_quad, each = 4), rep(condition_colors_trip, each = 3)),
  pch = all_pch,
  cex = 2,
  gene.selection = "common"
)

# Add legend to the plot
legend(par("usr")[2], par("usr")[4], legend = c(conditions_quad, conditions_trip), col = c(condition_colors_quad, condition_colors_trip), pch = c(16,16,16,16,17,17,17,17), cex = 1, bty = "n")
par(opar)
