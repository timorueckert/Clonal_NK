# Updated downstream analysis including data from mtASAP5

path <- "~/"

### NKG2C+ vs NKG2C- in CMVneg: Add two extra donors (CMVneg3+4) ###

# Integrate the two CMV- donors with the rest
integrated <- readRDS(paste0(path, "integrated_processed"))
CMVpos1 <- subset(integrated, subset = donor == "CMVpos1")
CMVpos2 <- subset(integrated, subset = donor == "CMVpos2")
CMVpos3 <- subset(integrated, subset = donor == "CMVpos3")
CMVpos4 <- subset(integrated, subset = donor == "CMVpos4")
CMVneg1 <- subset(integrated, subset = donor == "CMVneg1")
CMVneg2 <- subset(integrated, subset = donor == "CMVneg2")
CMVneg3 <- readRDS(paste0(path, "CMVneg3_processed"))
CMVneg4 <- readRDS(paste0(path, "CMVneg4_processed"))

Single_donors <- list(CMVpos1, CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVneg3, CMVneg4)
rm(CMVpos1, CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVneg3, CMVneg4)
rm(integrated)


# Process single donors
Single_donors <- lapply(Single_donors, function(seu) {
  DefaultAssay(seu) <- "MACS2"
  return(seu)
})


Single_donors <- lapply(Single_donors, function(seu) {
  seu <- FindTopFeatures(seu, min.cutoff = 10)
  return(seu)
})


Single_donors <- lapply(Single_donors, function(seu) {
  seu <- RunTFIDF(seu)
  return(seu)
})

Single_donors <- lapply(Single_donors, function(seu) {
  seu <- RunSVD(seu)
  return(seu)
})

Single_donors <- lapply(Single_donors, function(seu) {
  seu <- RunUMAP(seu, reduction = "lsi", dims = 2:30)
  return(seu)
})

lapply(Single_donors, function(seu) {
  DimPlot(seu, group.by = "experiment")
})

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = Single_donors,
  anchor.features = rownames(Single_donors[[2]]),
  reduction = "rlsi",
  dims = 2:30
)

# Merge and calculate LSI
All_merged <- merge(Single_donors[[1]], y = c(Single_donors[[2]],
                                                 Single_donors[[3]],
                                                 Single_donors[[4]],
                                                 Single_donors[[5]],
                                                 Single_donors[[6]],
                                                 Single_donors[[7]],
                                                 Single_donors[[8]]))
DefaultAssay(All_merged) <- "MACS2"
All_merged <- FindTopFeatures(All_merged, min.cutoff = 10)
All_merged <- RunTFIDF(All_merged)
All_merged <- RunSVD(All_merged)

# integrate LSI embeddings
All_integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = All_merged[["lsi"]],
  anchor.features = rownames(Single_donors[[2]]),
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

All_integrated <- RunUMAP(All_integrated,
                             reduction = "integrated_lsi", dims = 2:10)

rm(Single_donors)
# Normalize gene activity
All_integrated <- NormalizeData(
  object = All_integrated,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(All_integrated$nCount_RNA_atac)
)

# Subset CMVneg
CMVneg_integrated <- subset(All_integrated, subset = serostatus == "CMVneg")

# Dimred and clustering
CMVneg_integrated  <- RunUMAP(CMVneg_integrated,
                          reduction = "integrated_lsi", dims = 2:10, min.dist = 0.4)

CMVneg_integrated <- FindNeighbors(CMVneg_integrated, reduction = "integrated_lsi", dims = 2:10)
CMVneg_integrated <- FindClusters(CMVneg_integrated, resolution = 0.1)
DimPlot(CMVneg_integrated)
# Subcluster to resolve early dim
FeaturePlot(CMVneg_integrated, features = c("CD56", "CD16", "CD62L", "CD57"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE)&NoLegend()&NoAxes()
temp <- subset(CMVneg_integrated, idents = "1")
temp <- FindNeighbors(temp, reduction = "integrated_lsi", dims = 2:10)
temp  <- FindClusters(temp, resolution = 0.1)
DimPlot(temp)
Idents(CMVneg_integrated, cells = WhichCells(temp, idents = "0")) <- "EarlyCD56dim"
Idents(CMVneg_integrated, cells = WhichCells(temp, idents = "1")) <- "CD56bright"
CMVneg_integrated <- RenameIdents(CMVneg_integrated, "0" = "CD56dim")
levels(CMVneg_integrated) <- c("CD56bright", "EarlyCD56dim", "CD56dim")
DimPlot(CMVneg_integrated, cols = clusterpal)
FeaturePlot(CMVneg_integrated, features = "Anti-PE", cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

saveRDS(CMVneg_integrated, paste0(path, "CMVneg_processed_new"))

# Repeat also CMVpos as higher cell numbers in CMVneg enable less downsampling
# => higher power to detect differential gene scores
CMVpos <- readRDS(paste0(path, "CMVpos_processed"))

pdf(file = "Plots/new/Fig2_ATAC_Dimplot_CMVneg_updated.pdf", paper = "a4")
DimPlot(CMVneg_integrated, cols = c(clusterpal[1:3], adaptivepal))&NoLegend()
dev.off()

# Plot NKG2C from ASAP experiments
ASAP_CMVpos <- WhichCells(CMVpos, expression = experiment%in%c("mtASAP1", "mtASAP2"))
ASAP_CMVneg <- WhichCells(CMVneg_integrated, expression = experiment%in%c("mtASAP1", "mtASAP2", "mtASAP5"))

pdf(file = "Plots/new/Fig2_ATAC_NKG2C_CMVneg_updated.pdf", paper = "a4")
FeaturePlot(CMVneg_integrated, features = "Anti-PE",
            cells = ASAP_CMVneg, cols = colorscale, min.cutoff = "q1",
            max.cutoff = "q99",order = TRUE, pt.size = 0.4)&NoLegend()&NoAxes()&ggtitle(label = NULL)
dev.off()

#### Differential gene activity between NKG2C+ vs NKG2C-

mtASAP_CMVneg <- subset(CMVneg_integrated, cells = ASAP_CMVneg)
mtASAP_CMVpos <- subset(CMVpos, cells = ASAP_CMVpos)

RidgePlot(mtASAP_CMVneg, features = "Anti-PE")
RidgePlot(mtASAP_CMVpos, features = "Anti-PE")

mtASAP_CMVneg$population <- "between"
mtASAP_CMVpos$population <- "between"

mtASAP_CMVneg$population[mtASAP_CMVneg[["ADT"]]@data["Anti-PE",]>0.5] <- "NKG2Cpos"
mtASAP_CMVneg$population[mtASAP_CMVneg[["ADT"]]@data["Anti-PE",]<0.25] <- "NKG2Cneg"

mtASAP_CMVpos$population[mtASAP_CMVpos[["ADT"]]@data["Anti-PE",]>0.5] <- "NKG2Cpos"
mtASAP_CMVpos$population[mtASAP_CMVpos[["ADT"]]@data["Anti-PE",]<0.25] <- "NKG2Cneg"
table(mtASAP_CMVneg$population)
table(mtASAP_CMVpos$population)

# Downsample to equal numbers for DAP analysis
sample_cells_CMVneg <- c(sample(Cells(mtASAP_CMVneg[,mtASAP_CMVneg$population=="NKG2Cpos"]), size = 4300),
                         sample(Cells(mtASAP_CMVneg[,mtASAP_CMVneg$population=="NKG2Cneg"]), size = 4300))

sample_cells_CMVpos <- c(sample(Cells(mtASAP_CMVpos[,mtASAP_CMVpos$population=="NKG2Cpos"]), size = 4300),
                         sample(Cells(mtASAP_CMVpos[,mtASAP_CMVpos$population=="NKG2Cneg"]), size = 4300))
# export those for consistency
saveRDS(sample_cells_CMVpos, "CMVposdownsampled")
saveRDS(sample_cells_CMVneg, "CMVnegdownsampled")

# Read in
sample_cells_CMVpos <- readRDS(paste0(path, "Clonal_NK1/CMVposdownsampled"))
sample_cells_CMVneg <- readRDS(paste0(path, "Clonal_NK1/CMVnegdownsampled"))

mtASAP_CMVneg_down <- subset(mtASAP_CMVneg, cells = sample_cells_CMVneg)
mtASAP_CMVpos_down <- subset(mtASAP_CMVpos, cells = sample_cells_CMVpos)


DefaultAssay(mtASAP_CMVneg_down) <- "RNA_atac"
DefaultAssay(mtASAP_CMVpos_down) <- "RNA_atac"
Idents(mtASAP_CMVneg_down) <- mtASAP_CMVneg_down$population
Idents(mtASAP_CMVpos_down) <- mtASAP_CMVpos_down$population
activity2C_CMVpos <- FindMarkers(mtASAP_CMVpos_down, ident.1 = "NKG2Cpos", ident.2 = "NKG2Cneg", logfc.threshold = 0)
activity2C_CMVneg <- FindMarkers(mtASAP_CMVneg_down, ident.1 = "NKG2Cpos", ident.2 = "NKG2Cneg", logfc.threshold = 0)
activity2C_CMVpos$gene <- rownames(activity2C_CMVpos)
activity2C_CMVneg$gene <- rownames(activity2C_CMVneg)

# Number of significantly different genes
table(activity2C_CMVpos$p_val_adj<0.05 & abs(activity2C_CMVpos$avg_log2FC)>0.1)
table(activity2C_CMVneg$p_val_adj<0.05 & abs(activity2C_CMVneg$avg_log2FC)>0.1)

# Add Diffexpr column
activity2C_CMVpos$diffexp <- "NO"
activity2C_CMVpos$diffexp[activity2C_CMVpos$p_val_adj < 0.05 & activity2C_CMVpos$avg_log2FC > 0.1] <- "UP"
activity2C_CMVpos$diffexp[activity2C_CMVpos$p_val_adj < 0.05 & activity2C_CMVpos$avg_log2FC < -0.1] <- "DOWN"

activity2C_CMVneg$diffexp <- "NO"
activity2C_CMVneg$diffexp[activity2C_CMVneg$p_val_adj < 0.05 & activity2C_CMVneg$avg_log2FC > 0.1] <- "UP"
activity2C_CMVneg$diffexp[activity2C_CMVneg$p_val_adj < 0.05 & activity2C_CMVneg$avg_log2FC < -0.1] <- "DOWN"

# Layer so that significant genes are on top
layer2_CMVpos <- activity2C_CMVpos[activity2C_CMVpos$diffexp%in%c("UP", "DOWN"),]
layer1_CMVpos <- activity2C_CMVpos[activity2C_CMVpos$diffexp%in%c("NO"),]

layer2_CMVneg <- activity2C_CMVneg[activity2C_CMVneg$diffexp%in%c("UP", "DOWN"),]
layer1_CMVneg <- activity2C_CMVneg[activity2C_CMVneg$diffexp%in%c("NO"),]

# Plot
genestolabel <- c("ZBTB16", "ZBTB38", "SYK", "IL12RB2", "IFNG", "CADM1", "ZNF516", "ZEB2", "AIM2", "RPTOR", "JAKMIP1", "PARD6G", "CADM1")
ggplot()+
  geom_point(data = layer1_CMVpos, aes(x = avg_log2FC, y = -log10(p_val_adj)), color = "black")+
  geom_point(data = layer2_CMVpos, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene),
                  data = activity2C_CMVpos[activity2C_CMVpos$gene%in%c(genestolabel)&!activity2C_CMVpos$diffexp=="NO",],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 150))+
  scale_x_continuous(limits = c(-0.6, 0.6))+
  xlab(label = "NKG2Cpos vs NKG2Cneg (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

ggplot()+
  geom_point(data = layer1_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj)), color = "black")+
  geom_point(data = layer2_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene),
                  data = activity2C_CMVneg[activity2C_CMVneg$gene%in%c(genestolabel)&!activity2C_CMVneg$diffexp=="NO",],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 150))+
  scale_x_continuous(limits = c(-0.6, 0.6))+
  xlab(label = "NKG2Cpos vs NKG2Cneg (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()





# Plot representative gene activities associated to some DAPs
DefaultAssay(CMVneg_integrated) <- "RNA_atac"
cells_downsampled2 <- sample(Cells(CMVneg_integrated), size = 3500)

pdf(file = "Plots/new/Fig2_ATAC_representatives_CMVneg.pdf")
FeaturePlot(CMVneg_integrated, features = c("JAKMIP1", "ZBTB38", "ZNF516", "ZBTB16"),
            cols = blueyellow,  pt.size = 0.2, min.cutoff = "q1", max.cutoff = "q95", coord.fixed = T,order = TRUE,cells = cells_downsampled2)&NoLegend()&NoAxes()
FeaturePlot(CMVneg_integrated, features = c("JAKMIP1", "ZBTB38", "ZNF516", "ZBTB16"),
            cols = blueyellow,  pt.size = 0.2, min.cutoff = "q1", max.cutoff = "q95", coord.fixed = T,order = F,cells = cells_downsampled2)&NoLegend()&NoAxes()
dev.off()





