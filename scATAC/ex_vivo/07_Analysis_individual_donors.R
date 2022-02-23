##### Donor-specific processing as preparation for subcluster 
##### analysis and mitochondrial genotyping ######

path <- "~/"

# Import datasets
integrated <- readRDS(paste0(path, "integrated_processed"))

# Analyze individual donors
DefaultAssay(integrated) <- "MACS2"

# Keep only datasets with mitochondrial information
Human_NK_genotyping <- subset(integrated, experiment%in%c("CMVpos2", "CMVneg2", "mtASAP2"))
rm(integrated)


CMVpos2 <- subset(Human_NK_genotyping, donor == "CMVpos2")
CMVpos3 <- subset(Human_NK_genotyping, donor == "CMVpos3")
CMVpos4 <- subset(Human_NK_genotyping, donor == "CMVpos4")

Human_NK_genotyping <- list(CMVpos2, CMVpos3, CMVpos4)
rm(CMVpos2, CMVpos3, CMVpos4)

# Normalize
Human_NK_genotyping <- lapply(Human_NK_genotyping, function(seu) {
  seu <- RunTFIDF(seu)
  return(seu)
})

Human_NK_genotyping <- lapply(Human_NK_genotyping, function(seu) {
  seu <- FindTopFeatures(seu, min.cutoff = "5")
  return(seu)
})

Human_NK_genotyping <- lapply(Human_NK_genotyping, function(seu) {
  seu <- RunSVD(seu)
  return(seu)
})

# Check for correlation with depth to select dims for further Analysis
DepthCor(Human_NK_genotyping[[1]])
DepthCor(Human_NK_genotyping[[2]])
DepthCor(Human_NK_genotyping[[3]])


##### CMVpos4/3 were multiplexed in experiment mtASAP2 => clean-up using homoplasmic
# mutations to correct for imperfect hashtagging 

#### Mitochondrial Genotyping to improve demultiplexing #####

CMVpos3 <- Human_NK_genotyping[[2]]
CMVpos4 <- Human_NK_genotyping[[3]]

# Filter cells with good mitochondrial coverage
VlnPlot(CMVpos3, features = "mtDNA_depth")
CMVpos3 <- subset(CMVpos3, mtDNA_depth > 5)

VlnPlot(CMVpos4, features = "mtDNA_depth")
CMVpos4 <- subset(CMVpos4, mtDNA_depth > 5)

# Calling of variable sites and selection of variants with high confidence
mito.refallele <- readRDS(paste0(path, "mito.refallele"))

variable.sites_CMVpos3 <- IdentifyVariants(CMVpos3, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos3, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 2)


variable.sites_CMVpos4 <- IdentifyVariants(CMVpos4, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos4, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 2)


# Select ubiquituous mutations which can serve to examine/improve hashtag-based demux
ubiq_CMVpos3 <- subset(
  variable.sites_CMVpos3, subset = n_cells_over_20 >= 3800 &
    strand_correlation >= 0.65 &
    vmr < 0.01 & vmr > 0.001
)


ubiq_CMVpos4 <- subset(
  variable.sites_CMVpos4, subset = n_cells_over_20 >= 3060 &
    strand_correlation >= 0.65 &
    vmr < 0.01 & vmr > 0.001
)



# Are there any overlapping mutations? No!
intersect(ubiq_CMVpos3$variant, ubiq_CMVpos4$variant)

# Take only the X most abundant mutations
top_ubiq_CMVpos3 <- ubiq_CMVpos3 %>% top_n(5, wt = mean)
top_ubiq_CMVpos4 <- ubiq_CMVpos4 %>% top_n(5, wt = mean)

# Compute variant frequency per cell for ubiquitous features
CMVpos3 <- AlleleFreq(
  object = CMVpos3,
  variants = c(top_ubiq_CMVpos3$variant, top_ubiq_CMVpos4$variant),
  assay = "mito"
)

CMVpos4 <- AlleleFreq(
  object = CMVpos4,
  variants = c(top_ubiq_CMVpos3$variant, top_ubiq_CMVpos4$variant),
  assay = "mito"
)

CMVpos3[["alleles"]]
CMVpos4[["alleles"]]

DefaultAssay(CMVpos3) <- "alleles"
DefaultAssay(CMVpos4) <- "alleles"


# # # Find Clonotypes based on ubiq. mutations to identify contaminants and exclude them;
# # They are very few, so k has to be low and clustering resolution quite high
CMVpos3 <- FindClonotypes(CMVpos3, k = 2)
CMVpos3 <- FindClusters(CMVpos3, resolution = 1, group.singletons = F)
VlnPlot(CMVpos3, features = c(top_ubiq_CMVpos3$variant, top_ubiq_CMVpos4$variant), ncol = 4) # Cluster 8 contains contaminating cells
table(Idents(CMVpos3)) # 3 cells
CMVpos3 <- subset(CMVpos3, idents = "8", invert = TRUE)



CMVpos4 <- FindClonotypes(CMVpos4, k = 2)
CMVpos4 <- FindClusters(CMVpos4, resolution = 2, group.singletons = F)
VlnPlot(CMVpos4, features = c(top_ubiq_CMVpos3$variant, top_ubiq_CMVpos4$variant)) # Cluster 8 contains contaminating cells
table(Idents(CMVpos4)) # 27 cells
CMVpos4 <- subset(CMVpos4, idents = "9", invert = TRUE)


###### Pre-processing of CMVpos3  #####
# (CMVpos2/4 are analyzed together with later time point and loaded from there)

DefaultAssay(CMVpos3) <- "MACS2"
CMVpos3 <- RunUMAP(CMVpos3, reduction = 'lsi', dims = 2:30, n.neighbors = 20)
CMVpos3 <- FindNeighbors(CMVpos3, dims = 2:30, reduction = "lsi")
CMVpos3 <- FindClusters(CMVpos3, resolution = 0.4)


# Subcluster to identify early dim

DefaultAssay(CMVpos3) <- "ADT"
FeaturePlot(CMVpos3, features = c("CD56", "CD57", "anti-Biotin", "CD62L"), cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order =TRUE)&NoLegend()
DefaultAssay(CMVpos3) <- "MACS2"
temp <- subset(CMVpos3, idents = 0)
temp  <- FindNeighbors(temp,  dims = 2:30, reduction = "lsi")
temp <- FindClusters(temp, resolution = 0.4)
DimPlot(temp)
Idents(CMVpos3, cells= WhichCells(temp, idents = "3")) <- "EarlyCD56dim"
DimPlot(CMVpos3)


DefaultAssay(Human_NK_genotyping[[3]]) <- "MACS2"

# Annotate
CMVpos3 <- RenameIdents(CMVpos3, "4" = "CD56bright", "0" = "CD56dim", "1" = "Adaptive1", "2" = "Adaptive2", "3" = "Adaptive3")
levels(CMVpos3) <- c("CD56bright", "EarlyCD56dim", "CD56dim", "Adaptive1", "Adaptive2", "Adaptive3")
DimPlot(CMVpos3, label = TRUE)
CMVpos3$annotation <- Idents(CMVpos3)


# Export data processed per individual
saveRDS(CMVpos2, file = "CMVpos2_processed")
saveRDS(CMVpos3, file = "CMVpos3_processed")
saveRDS(CMVpos4, file = "CMVpos4_processed")

# Join early and late time points of CMVpos2/4
CMVpos2_late <- readRDS("CMVpos2_late_processed")
CMVpos4_late <- readRDS("CMVpos4_late_processed")


# Merge early and late time points
CMVpos2_full <- merge(CMVpos2_early, CMVpos2_late)
CMVpos4_full <- merge(CMVpos4_early, CMVpos4_late)

DefaultAssay(CMVpos2_full) <- "MACS2"
DefaultAssay(CMVpos4_full) <- "MACS2"

# LSI on merged dataset
CMVpos2_full <- RunTFIDF(CMVpos2_full)
CMVpos2_full <- FindTopFeatures(CMVpos2_full, min.cutoff = '5')
CMVpos2_full <- RunSVD(CMVpos2_full)

CMVpos4_full <- RunTFIDF(CMVpos4_full)
CMVpos4_full <- FindTopFeatures(CMVpos4_full, min.cutoff = '5')
CMVpos4_full <- RunSVD(CMVpos4_full, )

# UMAP
DepthCor(CMVpos2_full)
CMVpos2_full <- RunUMAP(object = CMVpos2_full, reduction = 'lsi', dims = 2:30)

DepthCor(CMVpos4_full)
CMVpos4_full <- RunUMAP(object = CMVpos4_full, reduction = 'lsi', dims = 2:30)

# Inspect UMAP
DimPlot(CMVpos2_full, group.by = "experiment", shuffle = T)
DimPlot(CMVpos4_full, group.by = "experiment", shuffle = T)


# There are some  small batch effects => correct with harmony,
library(harmony)
CMVpos2_full <- RunHarmony(
  object = CMVpos2_full ,
  group.by.vars = 'experiment',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE,
  epsilon.harmony = -Inf
)


# since effects are minimal for CMVpos4, run simple linear regression to avoid integration artefacts (nclust = 1) 
CMVpos4_full <- RunHarmony(
  object = CMVpos4_full ,
  group.by.vars = 'experiment',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE,
  epsilon.harmony = -Inf,
  nclust = 1
)

CMVpos2_full <- RunUMAP(object = CMVpos2_full, reduction = 'harmony', dims = 2:30)
CMVpos4_full <- RunUMAP(object = CMVpos4_full, reduction = 'harmony', dims = 2:30)

# Inspect UMAP
DimPlot(CMVpos2_full, group.by = "experiment", shuffle = T)
DimPlot(CMVpos4_full, group.by = "experiment", shuffle = T)

DimPlot(CMVpos2_full, split.by = "experiment", shuffle = T, group.by = "experiment")
DimPlot(CMVpos4_full, split.by = "experiment", shuffle = T, group.by = "experiment")


# CMVpos2
CMVpos2_full <- RunUMAP(object = CMVpos2_full, reduction = 'harmony', dims = 2:30, n.neighbors = 20)
CMVpos2_full <- FindNeighbors(CMVpos2_full, reduction = "harmony", dims = 2:30)


CMVpos2_full <- FindClusters(CMVpos2_full, verbose = FALSE, algorithm = 3, resolution = 0.3)
DimPlot(CMVpos2_full, label = TRUE)

FeaturePlot(CMVpos2_full, features = c("CD56", "anti-Biotin", "CD16", "CD62L", "CD127", "CD57"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE, cells = WhichCells(CMVpos2_full, expression = experiment == "mtASAP3"))

# Exlude CD127 ILCs (Cluster 5) and recalculate umap
CMVpos2_full <- subset(CMVpos2_full, idents = "5", invert = TRUE)
CMVpos2_full <- RunUMAP(object = CMVpos2_full, reduction = 'harmony', dims = 2:30, n.neighbors = 20)
DimPlot(CMVpos2_full, label = TRUE)

# Subcluster to resolve early dim
temp <- subset(CMVpos2_full, idents = "3")
temp <- FindNeighbors(temp, reduction = "harmony", dims = 2:30)
temp <- FindClusters(temp, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(temp)
Idents(CMVpos2_full, cells = WhichCells(temp, idents = c("1", "2"))) <- "EarlyCD56dim"
DimPlot(CMVpos2_full)



# CMVpos4
CMVpos4_full <- RunUMAP(object = CMVpos4_full, reduction = 'harmony', dims = 2:30, n.neighbors = 20)

CMVpos4_full <- FindNeighbors(CMVpos4_full, reduction = "harmony", dims = 2:30)
CMVpos4_full <- FindClusters(CMVpos4_full, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(CMVpos4_full, label = TRUE)
FeaturePlot(CMVpos4_full, features = c("CD56", "anti-Biotin", "CD16", "CD62L", "CD127", "CD57"), cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE)

# Subcluster the few early CD56dim
temp <- subset(CMVpos4_full, idents = "4")

temp <- FindNeighbors(temp, reduction = "harmony", dims = 2:30, k.param = 5)
temp <- FindClusters(temp, verbose = FALSE, algorithm = 3, resolution = 0.3)
DimPlot(temp)
Idents(CMVpos4_full, cells = WhichCells(temp, idents = 1)) <- "EarlyCD56dim"
DimPlot(CMVpos4_full, label = TRUE)

# Annotate clusters
CMVpos2_full <- RenameIdents(CMVpos2_full, "3" = "CD56bright", "0" = "CD56dim", "1" = "Adaptive1", "2" = "Adaptive2", "4" = "Adaptive3")
CMVpos4_full <- RenameIdents(CMVpos4_full, "4" = "CD56bright", "0" = "CD56dim", "1" = "Adaptive1", "2" = "Adaptive2", "3" = "Adaptive3")

levels(CMVpos2_full) <- c("CD56bright", "EarlyCD56dim", "CD56dim", "Adaptive1", "Adaptive2", "Adaptive3")
levels(CMVpos4_full) <- c("CD56bright", "EarlyCD56dim", "CD56dim", "5", "Adaptive1", "Adaptive2", "Adaptive3")

CMVpos2_full$annotation <- Idents(CMVpos2_full)
CMVpos4_full$annotation <- Idents(CMVpos4_full)

# Add RNA assay to merged objects
DefaultAssay(CMVpos2_full) <- "MACS2"
gene.activities_adaptive <- GeneActivity(CMVpos2_full)
CMVpos2_full[["RNA_atac"]] <- CreateAssayObject(counts = gene.activities_adaptive)
DefaultAssay(CMVpos2_full) <- "RNA_atac"
CMVpos2_full <- NormalizeData(
  object = CMVpos2_full,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(CMVpos2_full$nCount_RNA_atac)
)
# 
DefaultAssay(CMVpos4_full) <- "MACS2"
Annotation(CMVpos4_full) <- Annotation(integrated)
gene.activities_adaptive <- GeneActivity(CMVpos4_full)
CMVpos4_full[["RNA_atac"]] <- CreateAssayObject(counts = gene.activities_adaptive)
DefaultAssay(CMVpos4_full) <- "RNA_atac"
CMVpos4_full <- NormalizeData(
  object = CMVpos4_full,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(CMVpos4_full$nCount_RNA_atac)
)

# Export
saveRDS(CMVpos2_full, file = paste0(path, "CMVpos2_full"))
saveRDS(CMVpos4_full, file = paste0(path, "CMVpos4_full"))


#### CMVpos1 ####

# Process first experiment
CMVpos1_early <- subset(integrated, subset = donor == "CMVpos1")
CMVpos1_early <- RunTFIDF(CMVpos1_early)
CMVpos1_early <- FindTopFeatures(CMVpos1_early, min.cutoff = 20)
CMVpos1_early <- RunSVD(CMVpos1_early)
CMVpos1_early <- RunUMAP(CMVpos1_early, reduction = "lsi",
                         dims = 2:30)
CMVpos1_early <- FindNeighbors(CMVpos1_early, reduction = "lsi", dims = 2:30)
CMVpos1_early <- FindClusters(CMVpos1_early, resolution = 0.7)
DimPlot(CMVpos1_early, label = T)

# Import second experiment
CMVpos1_late <- readRDS(paste0(path, "CMVpos1_late"))

# Merge and integrate
CMVpos1_full <- merge(CMVpos1_early, CMVpos1_late)
CMVpos1_full <- RunTFIDF(CMVpos1_full)

# Use overlapping variable features reduce batch effects
CMVpos1_early <- FindTopFeatures(CMVpos1_early, min.cutoff = 30)
CMVpos1_late <- FindTopFeatures(CMVpos1_late, min.cutoff = 20)
var_features <- intersect(VariableFeatures(CMVpos1_early),VariableFeatures(CMVpos1_late))
VariableFeatures(CMVpos1_full) <- var_features


CMVpos1_full <- RunSVD(CMVpos1_full)
CMVpos1_full <- RunUMAP(CMVpos1_full, reduction = "lsi", dims = 2:30)
DimPlot(CMVpos1_full, group.by = "experiment")

# Integrate with harmony
library(harmony)
CMVpos1_full <- RunHarmony(
  object = CMVpos1_full,
  group.by.vars = 'experiment',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE,
  epsilon.harmony = -Inf,
  nclust = 1
)
# Simple linear regression is not sufficient
CMVpos1_full <- RunUMAP(CMVpos1_full, reduction = "harmony", dims = 2:30)
DimPlot(CMVpos1_full, group.by = "experiment")


# Integrate without setting nclust = 1
library(harmony)
CMVpos1_full <- RunHarmony(
  object = CMVpos1_full,
  group.by.vars = 'experiment',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE,
  epsilon.harmony = -Inf
)
CMVpos1_full <- RunUMAP(CMVpos1_full, reduction = "harmony", dims = 2:30)
DimPlot(CMVpos1_full, group.by = "experiment")

DefaultAssay(CMVpos1_full) <- "ADT"
FeaturePlot(CMVpos1_full, features = c("CD56", "CD62L", "CD57", "anti-Biotin",
                                       "Anti-PE", "KIR2DL1-S1-S3-S5", "KIR2DL2-L3",
                                       "KIR3DL1", "NKp30"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99",
            order = TRUE)&NoLegend()&NoAxes()

DefaultAssay(CMVpos1_full) <- "MACS2"
CMVpos1_full <- FindNeighbors(CMVpos1_full, reduction = "harmony",
                              dims = 2:30)
CMVpos1_full <- FindClusters(CMVpos1_full, resolution = 0.4)
DimPlot(CMVpos1_full)

# Resolve early dim
temp <- subset(CMVpos1_full, idents = "5")
temp <- FindNeighbors(temp, reduction = "harmony",
                      dims = 2:30)
temp <- FindClusters(temp, resolution = 0.2)
DimPlot(temp)
Idents(CMVpos1_full, cells = WhichCells(temp, idents = "1")) <- "CD56bright"
Idents(CMVpos1_full, cells = WhichCells(temp, idents = "0")) <- "EarlyCD56dim"

#Annotate
CMVpos1_full <- RenameIdents(CMVpos1_full, "0" = "CD56dim", "1" = "Adaptive1",
                             "2" = "Adaptive2", "3" = "Adaptive3", "4" = "Adaptive4",
                             "6" = "Adaptive5")
levels(CMVpos1_full) <- c("CD56bright", "EarlyCD56dim", "CD56dim", "Adaptive1",
                          "Adaptive2", "Adaptive3","Adaptive4", "Adaptive5")
DimPlot(CMVpos1_full, cols = c(clusterpal[1:3], adaptivepal), pt.size = 0.4, shuffle = T)&NoAxes()

CMVpos1_full$annotation <- Idents(CMVpos1_full)


# Export
saveRDS(CMVpos1_full, file = paste0(path, "CMVpos1_full"))
CMVpos1_full <- readRDS( paste0(path, "CMVpos1_full"))


#### CMVpos3: Integration of both timepoints ####
CMVpos3_early <- readRDS(file = paste0(path, "CMVpos3_processed"))
CMVpos3_late <- readRDS(file = paste0(path, "CMVpos3_late"))
# Check phenotypic similarity of two time points before integration
DefaultAssay(CMVpos3_late) <- "ADT"
DefaultAssay(CMVpos3_early) <- "ADT"
FeaturePlot(CMVpos3_late, cols = colorscale, min.cutoff = "q1",
            max.cutoff = "q99", order = TRUE,
            features = c("Anti-PE", "anti-Biotin", "NKp30",
                         "KIR2DL1-S1-S3-S5", "KIR2DL2-L3", "KIR3DL1",
                         "CD328", "KLRG1", "CD2"))
FeaturePlot(CMVpos3_early, cols = colorscale, min.cutoff = "q1",
            max.cutoff = "q99", order = TRUE,
            features = c("Anti-PE", "anti-Biotin", "NKp30",
                         "KIR2DL1-S1-S3-S5", "KIR2DL2-L3", "KIR3DL1",
                         "CD328", "KLRG1", "CD2"))

# Two time points seem very similar => merge and integrate with harmony
DefaultAssay(CMVpos3_late) <- "MACS2"
DefaultAssay(CMVpos3_early) <- "MACS2"
CMVpos3_full <- merge(CMVpos3_early, CMVpos3_late)
CMVpos3_full <- RunTFIDF(CMVpos3_full)
CMVpos3_full <- FindTopFeatures(CMVpos3_full, min.cutoff = 10)
CMVpos3_full <- RunSVD(CMVpos3_full)
CMVpos3_full <- RunUMAP(CMVpos3_full, reduction = "lsi", dims = 2:30)
DimPlot(CMVpos3_full, group.by = "experiment")

# Already without integration, overlap between time points is nicely visible
# => Integrate
library(harmony)
CMVpos3_full <- RunHarmony(
  object = CMVpos3_full ,
  group.by.vars = 'experiment',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE,
  epsilon.harmony = -Inf,
  nclust = 1
)
CMVpos3_full <- RunUMAP(CMVpos3_full, reduction = "harmony", dims = 2:30, n.neighbors = 20)
DimPlot(CMVpos3_full, group.by = "experiment")
# Clustering
CMVpos3_full <- FindNeighbors(CMVpos3_full, reduction = "harmony",
                              dims = 2:30)
CMVpos3_full <- FindClusters(CMVpos3_full, resolution = 0.4)
DimPlot(CMVpos3_full)

# Resolve early dim
FeaturePlot(CMVpos3_full, features = c("CD56", "CD62L", "CD57", "anti-Biotin"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99",
            order = TRUE)
temp <- subset(CMVpos3_full, idents = "4")
temp <- FindNeighbors(temp, reduction = "harmony",
                      dims = 2:30)
temp <- FindClusters(temp, resolution = 0.4)
DimPlot(temp)
Idents(CMVpos3_full, cells = WhichCells(temp, idents = "0")) <- "CD56bright"
Idents(CMVpos3_full, cells = WhichCells(temp, idents = "1")) <- "EarlyCD56dim"
rm(temp)
DimPlot(CMVpos3_full)

# Annotate clusters and re-order levels
CMVpos3_full <- RenameIdents(CMVpos3_full, "0" = "CD56dim", "1" = "Adaptive1", "2" = "Adaptive2", "3" = "Adaptive3")
levels(CMVpos3_full) <- c("CD56bright", "EarlyCD56dim", "CD56dim", "Adaptive1", "Adaptive2", "Adaptive3")
DimPlot(CMVpos3_full, cols = c(clusterpal[1:3], adaptivepal), split.by = "experiment")
CMVpos3_full$annotation <- Idents(CMVpos3_full)

FeaturePlot(CMVpos3_full, cols = colorscale, min.cutoff = "q1",
            max.cutoff = "q99", order = TRUE,
            features = c("Anti-PE", "anti-Biotin"), split.by = "experiment")

# Export
saveRDS(CMVpos3_full, file = paste0(path, "CMVpos3_full"))
