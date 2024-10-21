TableOfCounts_mtor <- read.csv("/PATH/TO/FILE/TableOfCounts_mTOR_filtered.csv")
samples_mtor <- read.table("TGFb_Torin_sample.txt", header = TRUE)
#DefineGroups 
group_mtor <- factor(c(samples_mtor$Group))

DEGall_mtor <- DGEList(counts=TableOfCounts_mtor, group=group_mtor)
keep_mtor <- filterByExpr(DEGall_mtor)
DEGall_mtor <- DEGall_mtor[keep_mtor,,keep.lib.sizes=FALSE] 
DEGall_mtor <- calcNormFactors(DEGall_mtor)

#group design matrix
design_mtor <- model.matrix(~0+group_mtor)
design_mtor

DEGall_mtor <- estimateDisp(DEGall_mtor,design_mtor) # sometimes DEGall is called as NormData
fit <- glmQLFit(DEGall_mtor,design_mtor)
plotBCV(DEGall_mtor)

qlf.wt_vs_tgfb <- glmQLFTest(fit, contrast = c(-1,1,0,0))
qlf.wt_vs_tgfb$table$fdr <- p.adjust(qlf.wt_vs_tgfb$table$PValue,"fdr")
write.csv(qlf.wt_vs_tgfb, "/PATH/TO/OUTPUT/qlf.wt_vs_tgfb.csv")

wt_vs_tgfb<- as.data.frame(qlf.wt_vs_tgfb)
row.names(wt_vs_tgfb) <- wt_vs_tgfb$Symbol
genes_2_label <- c("ADH1B","NPTX1","BDKRB2","TWSG1","AKR1C1","CTSK","FER1L6","BCL2A1","RASL11A","SERPINE1","MTHFD2","ASNS","PSAT1","ACTA2","NOX4","ITGA11","IGFBP3","IL11")

wt_vs_tgfb_volcano <- EnhancedVolcano(wt_vs_tgfb,
    lab = rownames(wt_vs_tgfb),
    x = 'logFC',
    y = 'PValue',
    selectLab = genes_2_label,
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 4.0,
    labSize = 12.0,  # Increased font size
    col = c('black', 'blue', 'grey', 'red3'),
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    title = 'UT vs TGFb',
    axisLabSize = 20) 

wt_vs_tgfb_img <- "/PATH/TO/OUTPUT/Figure-1B.png"
plot_width <- 15
plot_height <- 15

ggsave(wt_vs_tgfb_img , plot = wt_vs_tgfb_volcano , width = plot_width, height = plot_height, limitsize = FALSE)
