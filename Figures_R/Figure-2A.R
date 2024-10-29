library(pheatmap)
library(edgeR)
library(ggplot2)

# Read sample metadata
samples_atf4 <- read.table("/PATH/TO/INPUT/ATF4KO_sample.txt", header = TRUE)

## Define Groups
# Create a factor for the experimental groups
group_atf4 <- factor(c(samples_atf4$Group))

## Create DGEList Object
# Construct a DGEList object for differential expression analysis
DEGall_atf4 <- DGEList(counts=TableOfCounts_atf4, group=group_atf4) # TableOfCounts_atf4 from Table_of_Counts_atf4.R

# Filter low-expressed genes
keep_atf4 <- filterByExpr(DEGall_atf4)
DEGall_atf4 <- DEGall_atf4[keep_atf4,,keep.lib.sizes=FALSE] 

# Normalize the data
DEGall_atf4 <- calcNormFactors(DEGall_atf4)

## Design Matrix and Dispersion Estimation
# Create the design matrix for the group comparison
design_atf4 <- model.matrix(~0+group_atf4)

# Estimate dispersions
DEGall_atf4 <- estimateDisp(DEGall_atf4,design_atf4)

## Fit the Model and Test
# Fit a negative binomial generalized log-linear model
fit_atf4 <- glmQLFit(DEGall_atf4,design_atf4)

# Plot Biological Coefficient of Variation (BCV)
plotBCV(DEGall_atf4)

# Perform a quasi-likelihood F-test for differential expression
qlf.ctrl_tgfb_vs_atf4ko_tgfb <- glmQLFTest(fit_atf4, contrast = c(0,1,0,-1))

# Adjust p-values for multiple testing using the False Discovery Rate (FDR)
qlf.ctrl_tgfb_vs_atf4ko_tgfb$table$fdr <- p.adjust(qlf.ctrl_tgfb_vs_atf4ko_tgfb$table$PValue,"fdr")

## Prepare Data
# Convert results to data frame
ctrl_tgfb_vs_atf4ko_tgfb<- as.data.frame(qlf.ctrl_tgfb_vs_atf4ko_tgfb)
row.names(ctrl_tgfb_vs_atf4ko_tgfb) <- ctrl_tgfb_vs_atf4ko_tgfb$Symbol

## Import Data
counts.in_all_atf4 <- TableOfCounts_atf4_modified # TableOfCounts_atf4_modified from Table_of_Counts_atf4.R

## Define Metadata
counts.metadata_all_atf4 <- data.frame(
  dataset= c(colnames(counts.in_all_atf4)),
  Treatment= c("siControl_UT","siControl_UT","siControl_UT",
               "siControl_TGFb","siControl_TGFb","siControl_TGFb",
               "siATF4_UT","siATF4_UT","siATF4_UT",
               "siATF4_TGFb","siATF4_TGFb","siATF4_TGFb"
               ),
  stringsAsFactors = FALSE
)

## Define Groups
group_all_atf4 <- counts.metadata_all_atf4$Treatment

## Create DGEList Object
y_all_atf4 <- DGEList(counts=counts.in_all_atf4,
                                  genes=row.names.data.frame(counts.in_all_atf4),
                                  group=group_all_atf4)

## Filter Low-Expressed Genes
keep_all_atf4 <- rowSums(edgeR::cpm(y_all_atf4)>1)>=1
table(keep_all_atf4) # False:10331 True:14653
y_all_atf4<- y_all_atf4[keep_all_atf4, , keep.lib.sizes=FALSE]

## Normalize the Data
y_all_atf4<- calcNormFactors(y_all_atf4, method = "TMM")
counts_all_atf4 <- as.matrix(y_all_atf4$counts)
logCPM_all_atf4 <- edgeR::cpm(counts_all_atf4, prior.count=1, log=TRUE)

## Calculate Z-Score
ZScore_all_atf4 <- t(scale(t(logCPM_all_atf4)))
ZScore_all_atf4 <- as.data.frame(ZScore_all_atf4)

## Subset Significant Genes
Sig_ctrl_tgfb_vs_atf4ko_tgfb <- subset.data.frame(ctrl_tgfb_vs_atf4ko_tgfb, abs(logFC) >= 0.5 & fdr<0.05)
ATF4_G <- Sig_ctrl_tgfb_vs_atf4ko_tgfb$Symbol

## Generate Z-Score Matrix for Significant Genes
Sig.zscore_atf4 <- ZScore_all_atf4[ATF4_G  ,]#change for each cluster
Sig.zscore.mat_atf4 <- as.matrix(Sig.zscore_atf4)
Sig.zscore.mat_atf4 <- Sig.zscore.mat_atf4[complete.cases(Sig.zscore.mat_atf4),]

## Prepare Heatmap Annotation
heat.annotation_atf4 <- data.frame(counts.metadata_all_atf4[,2])
colnames(heat.annotation_atf4) <- "Treatment"
row.names(heat.annotation_atf4) <- counts.metadata_all_atf4[,1]

## Set Column Names for Heatmap
d.colnames_atf4 <- c(counts.metadata_all_atf4[,1])
colnames(Sig.zscore.mat_atf4) <- d.colnames_atf4

## Italicize Gene Names
newnames_atf4 <- lapply(
  rownames(Sig.zscore.mat_atf4),
  function(x) bquote(italic(.(x))))

## Define Annotation Colors
ann.colors_all_atf4 <- list(Treatment= c(`siControl_UT` = "red",`siControl_TGFb`="pink",
                                         `siATF4_UT` = "blue", `siATF4_TGFb`="lightblue"))


## Create Heatmap
atf4_heatmap <- pheatmap(Sig.zscore.mat_atf4 ,
                         annotation_col = heat.annotation_atf4, 
                         cluster_cols = FALSE,
                         main="",
                         annotation_colors = ann.colors_all_atf4,
                         show_colnames = FALSE,
                         labels_row = as.expression(newnames_atf4),
                         fontsize_row = 12
                         
)

## Save Heatmap Plot
ggsave("/PATH/TO/OUTPUT/Figure-2A.png", plot = atf4_heatmap, width = 10, height = 15, limitsize = FALSE)
