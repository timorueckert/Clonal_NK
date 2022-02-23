## Mitochondrial genotyping and clonotype clustering



CMVpos2_full <- readRDS(paste0(path, "CMVpos2_full"))
CMVpos4_full <- readRDS(paste0(path, "CMVpos4_full"))
CMVpos3_full <- readRDS(paste0(path, "CMVpos3_full"))
CMVpos1_full <- readRDS(paste0(path, "CMVpos1_full"))




CMVpos2_full <- subset(CMVpos2_full, mtDNA_depth >= 5)
CMVpos3_full <- subset(CMVpos3_full, mtDNA_depth >= 5)
CMVpos4_full <- subset(CMVpos4_full, mtDNA_depth >= 5)


# Color palettes from archr
large_pal <- sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)], size = 100)

##### Mitochondrial genotyping ####

# Stash cluster annotations
CMVpos2_full$annotation <- Idents(CMVpos2_full)
CMVpos3_full$annotation <- Idents(CMVpos3_full)
CMVpos4_full$annotation <- Idents(CMVpos4_full)
CMVpos1_full$annotation <- Idents(CMVpos1_full)


#### CMVpos2

# Calling of variable sites and selection of variants with high confidence
variable.sites_CMVpos2 <- IdentifyVariants(CMVpos2_full, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos2, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
# Variable features
high.conf_CMVpos2 <- subset(
  variable.sites_CMVpos2, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)


# Compute variant frequency per cell
CMVpos2_full <- AlleleFreq(
  object = CMVpos2_full,
  variants = c(high.conf_CMVpos2$variant),
  assay = "mito"
)

CMVpos2_full[["alleles"]]

DefaultAssay(CMVpos2_full) <- "alleles"

#### CMVpos3
variable.sites_CMVpos3 <- IdentifyVariants(CMVpos3_full, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos3_full, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
high.conf_CMVpos3 <- subset(
  variable.sites_CMVpos3, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

# Compute variant frequency per cell
CMVpos3_full <- AlleleFreq(
  object = CMVpos3_full,
  variants = high.conf_CMVpos3$variant,
  assay = "mito"
)
CMVpos3_full[["alleles"]]

DefaultAssay(CMVpos3_full) <- "alleles"


#### CMVpos4
# Calling of variable sites and selection of variants with high confidence
variable.sites_CMVpos4 <- IdentifyVariants(CMVpos4_full, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos4, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
high.conf_CMVpos4 <- subset(
  variable.sites_CMVpos4, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

CMVpos4_full <- AlleleFreq(
  object = CMVpos4_full,
  variants = c(high.conf_CMVpos4$variant),
  assay = "mito"
)


CMVpos4_full[["alleles"]]

DefaultAssay(CMVpos4_full) <- "alleles"


## Cluster on allele frequencies to define clonotypes
## Clustering parameters are chosen so that clonotypes are defined by ~ 1 mutation

### CMVpos2
CMVpos2_full <- FindClonotypes(CMVpos2_full, k = 4, metric = "euclidean")
CMVpos2_full <- FindClusters(CMVpos2_full, resolution = 1.3, group.singletons = F, algorithm = 1)


# Stash Clonotypes for quantification
CMVpos2_full$clonotype <- Idents(CMVpos2_full)

# Find significantly enriched mutations to identify meaningful clonotypes
sigmut_CMVpos2  <- FindAllMarkers(CMVpos2_full, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut_CMVpos2 <- sigmut_CMVpos2  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)


DoHeatmap(CMVpos2_full, features = sigmut_CMVpos2$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos2_full, expression = clonotype %in% topmut_CMVpos2$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()

# Look at distribution of highly frequent mutations
FeaturePlot(
  object = CMVpos2_full,
  features = c("152T>C", "5590G>A", "11157T>C", "1415G>A"),
  order = TRUE,
  cols = c("grey", "darkred"),
  max.cutoff = "q90"
)&NoLegend()&NoAxes()

# Remove clonotypes described solely by differences in 152T>C frequencies
sigmut_CMVpos2  <- FindAllMarkers(CMVpos2_full, only.pos = T, logfc.threshold = 0.1,
                                  min.pct = 0.6,
                                  features = rownames(CMVpos2_full)[!rownames(CMVpos2_full)%in%"152T>C"])
topmut_CMVpos2 <- sigmut_CMVpos2  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)

DoHeatmap(CMVpos2_full, features = sigmut_CMVpos2$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos2_full, expression = clonotype %in% topmut_CMVpos2$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()

# Filter clonotypes
filtered_clonotypes_CMVpos2 <- topmut_CMVpos2$cluster[!topmut_CMVpos2$cluster%in%c("singleton")]

# Plot 152T>C as well!
DoHeatmap(CMVpos2_full, features = c(topmut_CMVpos2$gene, "152T>C"), slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos2_full, expression = clonotype %in% filtered_clonotypes_CMVpos2), draw.lines = F,
          group.colors = large_pal)+
  scale_fill_gradientn(colours = c("white", "darkred"))+NoLegend()

### CMVpos3
CMVpos3_full <- FindClonotypes(CMVpos3_full, k = 4, metric = "euclidean")
CMVpos3_full <- FindClusters(CMVpos3_full, resolution = 1.4, group.singletons = F, algorithm = 1)


# Stash Clonotypes for quantification
CMVpos3_full$clonotype <- Idents(CMVpos3_full)

# Find significantly enriched mutations to identify meaningful clonotypes
sigmut_CMVpos3  <- FindAllMarkers(CMVpos3_full, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut_CMVpos3 <- sigmut_CMVpos3  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)


DoHeatmap(CMVpos3_full, features = sigmut_CMVpos3$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos3_full, expression = clonotype %in% topmut_CMVpos3$cluster & clonotype!= "singleton"))+
  scale_fill_viridis_c()+NoLegend()

# Filter clonotypes
filtered_clonotypes_CMVpos3 <- topmut_CMVpos3$cluster[!topmut_CMVpos3$cluster%in%c("singleton")]

DoHeatmap(CMVpos3_full, features = topmut_CMVpos3$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype", 
          cells = WhichCells(CMVpos3_full, expression = clonotype %in% filtered_clonotypes_CMVpos3), draw.lines = F,
          group.colors = large_pal)+
  scale_fill_gradientn(colours = c("white", "darkred"))+NoLegend()


#### CMVpos4
CMVpos4_full <- FindClonotypes(CMVpos4_full, k = 4, metric = "euclidean")
CMVpos4_full <- FindClusters(CMVpos4_full, resolution = 1.2, group.singletons = F, algorithm = 1)


# Stash Clonotypes for quantification
CMVpos4_full$clonotype <- Idents(CMVpos4_full)

# Find significantly enriched mutations to identify meaningful clonotypes
sigmut_CMVpos4  <- FindAllMarkers(CMVpos4_full, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut_CMVpos4 <- sigmut_CMVpos4  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)

DoHeatmap(CMVpos4_full, features = sigmut_CMVpos4$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos4_full, expression = clonotype %in% topmut_CMVpos4$cluster & clonotype!= "singleton"))+
  scale_fill_viridis_c()+NoLegend()

FeaturePlot(
  object = CMVpos4_full,
  features = c("16155A>G"),
  order = TRUE,
  cols = c("grey", "darkred"),
  max.cutoff = "q50"
)&NoLegend()&NoAxes()

# Remove clonotypes described solely by differences in 16155A>G frequencies
sigmut_CMVpos4  <- FindAllMarkers(CMVpos4_full, only.pos = T, logfc.threshold = 0.1,
                                  min.pct = 0.6,
                                  features = rownames(CMVpos4_full)[!rownames(CMVpos4_full)%in%"16155A>G"])
topmut_CMVpos4 <- sigmut_CMVpos4  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)


# Filter clonotypes
filtered_clonotypes_CMVpos4 <- topmut_CMVpos4$cluster[!topmut_CMVpos4$cluster%in%c("singleton")]

DoHeatmap(CMVpos4_full, features = sigmut_CMVpos4$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos4_full, expression = clonotype %in% filtered_clonotypes_CMVpos4), draw.lines = F, group.colors = large_pal)+
  scale_fill_gradientn(colours = c("white", "darkred"))+NoLegend()



### CMVpos1
# Calling of variable sites and selection of variants with high confidence
VlnPlot(CMVpos1_full, features = "mtDNA_depth", group.by = "experiment")
CMVpos1_mito <- subset(CMVpos1_full, subset = mtDNA_depth > 5)
variable.sites_CMVpos1 <- IdentifyVariants(CMVpos1_mito, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos1[variable.sites_CMVpos1$vmr>0.000001,], concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
# Variable features
high.conf_CMVpos1 <- subset(
  variable.sites_CMVpos1, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)


# Compute variant frequency per cell
CMVpos1_mito <- AlleleFreq(
  object = CMVpos1_mito,
  variants = c(high.conf_CMVpos1$variant),
  assay = "mito"
)

CMVpos1_mito[["alleles"]]

DefaultAssay(CMVpos1_mito) <- "alleles"

### CMVpos1_mito
CMVpos1_mito <- FindClonotypes(CMVpos1_mito, k = 4, metric = "euclidean")
CMVpos1_mito <- FindClusters(CMVpos1_mito, resolution = 1.2, group.singletons = F, algorithm = 1)


# Stash Clonotypes for quantification
CMVpos1_mito$clonotype <- Idents(CMVpos1_mito)

# Find significantly enriched mutations to identify meaningful clonotypes
sigmut_CMVpos1  <- FindAllMarkers(CMVpos1_mito, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut_CMVpos1 <- sigmut_CMVpos1  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)


DoHeatmap(CMVpos1_mito, features = sigmut_CMVpos1$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos1_mito, expression = clonotype %in% topmut_CMVpos1$cluster & clonotype!= "singleton"))+
  scale_fill_viridis_c()+NoLegend()

# There is still a small contamination of cells from CMVpos3 => remove
VlnPlot(CMVpos1_mito, features = c(homoplasmic_muts_df$variant[40:52]), ncol = 3)
CMVpos1_mito <- subset(CMVpos1_mito, idents = "32", invert = TRUE)

# Repeat genotyping to only include informative mutations
# Calling of variable sites and selection of variants with high confidence
variable.sites_CMVpos1 <- IdentifyVariants(CMVpos1_mito, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos1, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
# Variable features
high.conf_CMVpos1 <- subset(
  variable.sites_CMVpos1, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)


# Compute variant frequency per cell
CMVpos1_mito <- AlleleFreq(
  object = CMVpos1_mito,
  variants = c(high.conf_CMVpos1$variant),
  assay = "mito"
)

CMVpos1_mito[["alleles"]]

DefaultAssay(CMVpos1_mito) <- "alleles"

### CMVpos1_mito
CMVpos1_mito <- FindClonotypes(CMVpos1_mito, k = 4, metric = "euclidean")
CMVpos1_mito <- FindClusters(CMVpos1_mito, resolution = 1, group.singletons = F, algorithm = 1)


# Stash Clonotypes for quantification
CMVpos1_mito$clonotype <- Idents(CMVpos1_mito)

# Find significantly enriched mutations to identify meaningful clonotypes
sigmut_CMVpos1  <- FindAllMarkers(CMVpos1_mito, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut_CMVpos1 <- sigmut_CMVpos1  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)


DoHeatmap(CMVpos1_mito, features = sigmut_CMVpos1$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos1_mito, expression = clonotype %in% topmut_CMVpos1$cluster & clonotype!= "singleton"))+
  scale_fill_viridis_c()+NoLegend()

# Look at distribution highly frequent mutations

FeaturePlot(
  object = CMVpos1_mito,
  features = c("5320C>T", "13710A>G"),
  order = TRUE,
  cols = c("grey", "darkred"),
  max.cutoff = "q80"
)&NoLegend()&NoAxes()

# Remove clonotypes described solely by differences in 5320T>C frequencies
sigmut_CMVpos1  <- FindAllMarkers(CMVpos1_mito, only.pos = T, logfc.threshold = 0.1,
                                  min.pct = 0.6,
                                  features = rownames(CMVpos1_mito)[!rownames(CMVpos1_mito)%in%"5320C>T"])
topmut_CMVpos1 <- sigmut_CMVpos1  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)

DoHeatmap(CMVpos1_mito, features = sigmut_CMVpos1$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(CMVpos1_mito, expression = clonotype %in% topmut_CMVpos1$cluster & clonotype!= "singleton"))+
  scale_fill_viridis_c()+NoLegend()


filtered_clonotypes_CMVpos1 <- topmut_CMVpos1$cluster[!topmut_CMVpos1$cluster%in%c("singleton")]

DoHeatmap(CMVpos1_mito, features = topmut_CMVpos1$gene, slot = "data", disp.min = 0.01, disp.max = 0.1, group.by = "clonotype", 
          cells = WhichCells(CMVpos1_mito, expression = clonotype %in% filtered_clonotypes_CMVpos1), draw.lines = F,
          group.colors = large_pal)+
  scale_fill_gradientn(colours = c("white", "darkred"))+NoLegend()

# Reset idents to annotations
Idents(CMVpos1_mito) <- CMVpos1_full$annotation
Idents(CMVpos2_full) <- CMVpos2_full$annotation
Idents(CMVpos3_full) <- CMVpos3_full$annotation
Idents(CMVpos4_full) <- CMVpos4_full$annotation


# Export
saveRDS(CMVpos2_full, file = paste0(path, "CMVpos2_full_alleles"))
saveRDS(CMVpos3_full, file = paste0(path, "CMVpos3_full_alleles"))
saveRDS(CMVpos4_full, file = paste0(path, "CMVpos4_full_alleles"))
saveRDS(CMVpos1_mito, file = paste0(path, "CMVpos1_mito"))


CMVpos2_full <- readRDS(paste0(path, "CMVpos2_full_alleles"))
CMVpos3_full <- readRDS(paste0(path, "CMVpos3_full_alleles"))
CMVpos4_full <- readRDS(paste0(path, "CMVpos4_full_alleles"))
CMVpos1_mito <- readRDS(paste0(path, "CMVpos1_mito"))





#### Export clonotype dataframes for quantification ####
# Put into dataframes for chisquare test

All_NK_clone_CMVpos1 <- data.frame(
  Clusters = Idents(CMVpos1_mito),
  mito_cluster = paste0(CMVpos1_mito$clonotype, "CMVpos1")
)
Adaptive_clone_CMVpos1 <- data.frame(
  Clusters = Idents(subset(CMVpos1_mito, idents = c("Adaptive1", "Adaptive2", "Adaptive3",
                                                    "Adaptive4", "Adaptive5"))),
  mito_cluster = paste0(subset(CMVpos1_mito, idents = c("Adaptive1", "Adaptive2", "Adaptive3",
                                                        "Adaptive4", "Adaptive5"))$clonotype, "CMVpos1")
)
Conventional_clone_CMVpos1 <- data.frame(
  Clusters = Idents(subset(CMVpos1_mito, idents = c("CD56bright", "EarlyCD56dim", "CD56dim"))),
  mito_cluster = paste0(subset(CMVpos1_mito, idents = c("CD56bright","EarlyCD56dim", "CD56dim"))$clonotype, "CMVpos1")
)

All_NK_clone_CMVpos2 <- data.frame(
  Clusters = Idents(CMVpos2_full),
  mito_cluster = paste0(CMVpos2_full$clonotype, "CMVpos2")
)
Adaptive_clone_CMVpos2 <- data.frame(
  Clusters = Idents(subset(CMVpos2_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3"))),
  mito_cluster = paste0(subset(CMVpos2_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3"))$clonotype, "CMVpos2")
)
Conventional_clone_CMVpos2 <- data.frame(
  Clusters = Idents(subset(CMVpos2_full, idents = c("CD56bright", "EarlyCD56dim", "CD56dim"))),
  mito_cluster = paste0(subset(CMVpos2_full, idents = c("CD56bright", "EarlyCD56dim", "CD56dim"))$clonotype, "CMVpos2")
)


All_NK_clone_CMVpos3 <- data.frame(
  Clusters = Idents(CMVpos3_full),
  mito_cluster = paste0(CMVpos3_full$clonotype, "CMVpos3")
)
Adaptive_clone_CMVpos3 <- data.frame(
  Clusters = Idents(subset(CMVpos3_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3"))),
  mito_cluster = paste0(subset(CMVpos3_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3"))$clonotype, "CMVpos3")
)
Conventional_clone_CMVpos3 <- data.frame(
  Clusters = Idents(subset(CMVpos3_full, idents = c("CD56bright", "EarlyCD56dim", "CD56dim"))),
  mito_cluster = paste0(subset(CMVpos3_full, idents = c("CD56bright", "EarlyCD56dim", "CD56dim"))$clonotype, "CMVpos3")
)


All_NK_clone_CMVpos4 <- data.frame(
  Clusters = Idents(CMVpos4_full),
  mito_cluster = paste0(CMVpos4_full$clonotype, "CMVpos4")
)
Adaptive_clone_CMVpos4 <- data.frame(
  Clusters = Idents(subset(CMVpos4_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3"))),
  mito_cluster = paste0(subset(CMVpos4_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3"))$clonotype, "CMVpos4")
)
Conventional_clone_CMVpos4 <- data.frame(
  Clusters = Idents(subset(CMVpos4_full, idents = c("CD56dim", "CD56bright"))),
  mito_cluster = paste0(subset(CMVpos4_full, idents = c("CD56bright",  "CD56dim"))$clonotype, "CMVpos4")
)




# Include only clonotypes which are clearly identified by mutational profile
All_NK_clone_CMVpos1 <- All_NK_clone_CMVpos1 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos1, "CMVpos1"))
Adaptive_clone_CMVpos1 <- Adaptive_clone_CMVpos1 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos1, "CMVpos1"))
Conventional_clone_CMVpos1 <- Conventional_clone_CMVpos1 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos1, "CMVpos1"))

All_NK_clone_CMVpos2 <- All_NK_clone_CMVpos2 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos2, "CMVpos2"))
Adaptive_clone_CMVpos2 <- Adaptive_clone_CMVpos2 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos2, "CMVpos2"))
Conventional_clone_CMVpos2 <- Conventional_clone_CMVpos2 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos2, "CMVpos2"))

All_NK_clone_CMVpos3 <- All_NK_clone_CMVpos3 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos3, "CMVpos3"))
Adaptive_clone_CMVpos3 <- Adaptive_clone_CMVpos3 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos3, "CMVpos3"))
Conventional_clone_CMVpos3 <- Conventional_clone_CMVpos3 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos3, "CMVpos3"))

All_NK_clone_CMVpos4 <- All_NK_clone_CMVpos4 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos4, "CMVpos4"))
Adaptive_clone_CMVpos4 <- Adaptive_clone_CMVpos4 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos4, "CMVpos4"))
Conventional_clone_CMVpos4 <- Conventional_clone_CMVpos4 %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes_CMVpos4, "CMVpos4"))



# Export clonotype tables
write.csv(All_NK_clone_CMVpos1, file = paste0(path, "All_NK_clone_CMVpos1"))
write.csv(Adaptive_clone_CMVpos1, file = paste0(path, "Adaptive_clone_CMVpos1"))
write.csv(Conventional_clone_CMVpos1, file = paste0(path, "Conventional_clone_CMVpos1"))
          
write.csv(All_NK_clone_CMVpos2, file = paste0(path, "All_NK_clone_CMVpos2"))
write.csv(Adaptive_clone_CMVpos2, file = paste0(path, "Adaptive_clone_CMVpos2"))
write.csv(Conventional_clone_CMVpos2, file = paste0(path, "Conventional_clone_CMVpos2"))

write.csv(All_NK_clone_CMVpos3, file = paste0(path, "All_NK_clone_CMVpos3"))
write.csv(Adaptive_clone_CMVpos3, file = paste0(path, "Adaptive_clone_CMVpos3"))
write.csv(Conventional_clone_CMVpos3, file = paste0(path, "Conventional_clone_CMVpos3"))

write.csv(All_NK_clone_CMVpos4, file = paste0(path, "All_NK_clone_CMVpos4"))
write.csv(Adaptive_clone_CMVpos4, file = paste0(path, "Adaptive_clone_CMVpos4"))
write.csv(Conventional_clone_CMVpos4, file = paste0(path, "Conventional_clone_CMVpos4"))


