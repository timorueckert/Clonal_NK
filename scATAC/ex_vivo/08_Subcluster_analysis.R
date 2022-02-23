# Analysis of adaptive subclusters per donor

path <- "~/"

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- met.brewer(name="Isfahan1",n=7,type="discrete")
clusterpal <- clusterpal[c(4,5,7,1)]
NKG2Cpal <- c("#40BAD5","#120136")
adaptivepal <- c("#880E4F","#C21807",  "#FF2400", "#7C0A02", "#ED2939", "#D21F3C", "#FA8072", "#BF0A30", "#BF0A30")
adaptivepal <- c("#880E4F","#7E191B", "#C0424E", "#E0115F","#BF0A30",  "#D21F3C",  "#FA8072","#CA3433",  "#CD5C5C" )
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

# From GreenleafLab/ArchR
blueyellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

# Read in data
CMVpos2_full <- readRDS(paste0(path, "CMVpos2_full"))
CMVpos4_full <- readRDS(paste0(path, "CMVpos4_full"))
CMVpos3_full <- readRDS(paste0(path, "CMVpos3_full"))
CMVpos1_full <- readRDS(paste0(path, "CMVpos1_full"))

# Dimplots for each donor
pdf(file = "Plots/new/Fig5_Dimplots.pdf", paper = "a4r")
DimPlot(CMVpos2_full, cols = c(clusterpal[1:3], adaptivepal[1:3]))&NoLegend()
DimPlot(CMVpos4_full, cols = c(clusterpal[1:3], adaptivepal[4:6]))&NoLegend()
DimPlot(CMVpos3_full, cols = c(clusterpal[1:3], adaptivepal[7:9]))&NoLegend()
DimPlot(CMVpos1_full, cols = c(clusterpal[1:3], adaptivepal))&NoLegend()
DimPlot(CMVpos2_full, cols = c(clusterpal[1:3], adaptivepal[1:3]))
DimPlot(CMVpos4_full, cols = c(clusterpal[1:3], adaptivepal[6:8]))
DimPlot(CMVpos3_full, cols = c(clusterpal[1:3], adaptivepal[c(4,5,9)]))
DimPlot(CMVpos1_full, cols = c(clusterpal[1:3], adaptivepal))
dev.off()

#### Add motif data for each donor #####
##### Re-run chromvar analysis for FOS::JUNB on a per-donor basis #####
library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE, collection = "CORE")
)
pfm

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix_D2 <- CreateMotifMatrix(
  features = granges(CMVpos2_full),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)
motif.matrix_D4 <- CreateMotifMatrix(
  features = granges(CMVpos4_full),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)
motif.matrix_D3 <- CreateMotifMatrix(
  features = granges(CMVpos3_full),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)

motif.matrix_D1 <- CreateMotifMatrix(
  features = granges(CMVpos1_full),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)

# Create a new Motif object to store the results
motif_D2 <- CreateMotifObject(
  data = motif.matrix_D2,
  pwm = pfm
)

motif_D3 <- CreateMotifObject(
  data = motif.matrix_D3,
  pwm = pfm
)

motif_D4 <- CreateMotifObject(
  data = motif.matrix_D4,
  pwm = pfm
)

motif_D1 <- CreateMotifObject(
  data = motif.matrix_D1,
  pwm = pfm
)


# Add the Motif object to the assay
CMVpos2_full <- SetAssayData(
  object = CMVpos2_full,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif_D2
)

CMVpos4_full <- SetAssayData(
  object = CMVpos4_full,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif_D4
)

CMVpos3_full <- SetAssayData(
  object = CMVpos3_full,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif_D3
)

CMVpos1_full <- SetAssayData(
  object = CMVpos1_full,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif_D3
)


# Calculate regional stats
CMVpos2_full <- RegionStats(object = CMVpos2_full, genome = BSgenome.Hsapiens.UCSC.hg38)
CMVpos4_full <- RegionStats(object = CMVpos4_full, genome = BSgenome.Hsapiens.UCSC.hg38)
CMVpos3_full <- RegionStats(object = CMVpos3_full, genome = BSgenome.Hsapiens.UCSC.hg38)
CMVpos1_full <- RegionStats(object = CMVpos1_full, genome = BSgenome.Hsapiens.UCSC.hg38)

CMVpos2_full <- RunChromVAR(
  object = CMVpos2_full,
  genome = BSgenome.Hsapiens.UCSC.hg38)
CMVpos4_full <- RunChromVAR(
  object = CMVpos4_full,
  genome = BSgenome.Hsapiens.UCSC.hg38)
CMVpos3_full <- RunChromVAR(
  object = CMVpos3_full,
  genome = BSgenome.Hsapiens.UCSC.hg38)
CMVpos1_full <- RunChromVAR(
  object = CMVpos1_full,
  genome = BSgenome.Hsapiens.UCSC.hg38)


# Plot FOS::JUNB motif 
pdf(file = "Plots/new/SupplFig5_AP1_activity.pdf", paper = "a4r")
VlnPlot(CMVpos2_full, features = "MA1134.1", cols = c(clusterpal, adaptivepal[1:3]), pt.size = 0)+scale_y_continuous(limits = c(-6,6))
VlnPlot(CMVpos4_full, features = "MA1134.1", cols = c(clusterpal, adaptivepal[4:6]), pt.size = 0)+scale_y_continuous(limits = c(-6,6))
VlnPlot(CMVpos3_full, features = "MA1134.1", cols = c(clusterpal, adaptivepal[7:9]), pt.size = 0)+scale_y_continuous(limits = c(-6,6))
VlnPlot(CMVpos1_full, features = "MA1134.1", cols = c(clusterpal, adaptivepal), pt.size = 0)+scale_y_continuous(limits = c(-6,6))
dev.off()

# Reviewer asked for statistical significance here...
#CMVpos2
AP1_df_P2 <- FetchData(CMVpos2_full, vars = c("MA1134.1", "annotation"))
#Get dim and adaptive
AP1_df_P2 <- AP1_df_P2 %>% dplyr::filter(!annotation%in%c("CD56bright", "EarlyCD56dim"))
# test
test <- pairwise.wilcox.test(AP1_df_P2$chromvar_MA1134.1, AP1_df_P2$annotation, p.adjust="BH")
test$p.value

#CMVpos3
AP1_df_P3 <- FetchData(CMVpos3_full, vars = c("MA1134.1", "annotation"))
#Get dim and adaptive
AP1_df_P3 <- AP1_df_P3 %>% dplyr::filter(!annotation%in%c("CD56bright", "EarlyCD56dim"))
# test
test1 <- pairwise.wilcox.test(AP1_df_P3$chromvar_MA1134.1, AP1_df_P3$annotation, p.adjust="BH")
test1$p.value

#CMVpos4
AP1_df_P4 <- FetchData(CMVpos4_full, vars = c("MA1134.1", "annotation"))
#Get dim and adaptive
AP1_df_P4 <- AP1_df_P4 %>% dplyr::filter(!annotation%in%c("CD56bright", "EarlyCD56dim"))
# test
test2 <- pairwise.wilcox.test(AP1_df_P4$chromvar_MA1134.1, AP1_df_P4$annotation, p.adjust="BH")
test2$p.value

#CMVpos1
AP1_df_P1 <- FetchData(CMVpos1_full, vars = c("MA1134.1", "annotation"))
#Get dim and adaptive
AP1_df_P1 <- AP1_df_P1 %>% dplyr::filter(!annotation%in%c("CD56bright", "EarlyCD56dim"))
# test
test3 <- pairwise.wilcox.test(AP1_df_P1$chromvar_MA1134.1, AP1_df_P1$annotation, p.adjust="BH")
test3$p.value

##### Find overlapping adaptive signature ####
DefaultAssay(CMVpos2_full) <- "MACS2"
DefaultAssay(CMVpos4_full) <- "MACS2"
DefaultAssay(CMVpos3_full) <- "MACS2"
DefaultAssay(CMVpos1_full) <- "MACS2"

# Find DAPs between adaptive and CD56dim
adaptive_DAP_CMVpos2 <- FindMarkers(
  object = CMVpos2_full,
  ident.1 = c("Adaptive1", "Adaptive2", "Adaptive3"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

adaptive_DAP_CMVpos3 <- FindMarkers(
  object = CMVpos3_full,
  ident.1 = c("Adaptive1", "Adaptive2", "Adaptive3"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

adaptive_DAP_CMVpos4 <- FindMarkers(
  object = CMVpos4_full,
  ident.1 = c("Adaptive1", "Adaptive2", "Adaptive3"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

adaptive_DAP_CMVpos1 <- FindMarkers(
  object = CMVpos1_full,
  ident.1 = c("Adaptive1", "Adaptive2", "Adaptive3", "Adaptive4", "Adaptive5"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

# Filter for statistically significant hits
adaptive_DAP_CMVpos2$gene <- rownames(adaptive_DAP_CMVpos2)
adaptive_DAP_CMVpos3$gene <- rownames(adaptive_DAP_CMVpos3)
adaptive_DAP_CMVpos4$gene <- rownames(adaptive_DAP_CMVpos4)
adaptive_DAP_CMVpos1$gene <- rownames(adaptive_DAP_CMVpos1)


adaptive_DAP_CMVpos2 <- adaptive_DAP_CMVpos2  %>% dplyr::filter(p_val_adj < 0.05)
adaptive_DAP_CMVpos3 <- adaptive_DAP_CMVpos3  %>% dplyr::filter(p_val_adj < 0.05)
adaptive_DAP_CMVpos4 <- adaptive_DAP_CMVpos4  %>% dplyr::filter(p_val_adj < 0.05)
adaptive_DAP_CMVpos1 <- adaptive_DAP_CMVpos1  %>% dplyr::filter(p_val_adj < 0.05)

# Generate Euler Plot
library(eulerr)
# Find relationships

all_peaks <- unique(c(rownames(adaptive_DAP_CMVpos1), rownames(adaptive_DAP_CMVpos2),
                      rownames(adaptive_DAP_CMVpos3),rownames(adaptive_DAP_CMVpos4)))

# Create matrix of logicals representing whether peak is a DAP for each dataset

peaks_CMVpos1 <- all_peaks%in%rownames(adaptive_DAP_CMVpos1)
peaks_CMVpos2 <- all_peaks%in%rownames(adaptive_DAP_CMVpos2)
peaks_CMVpos3 <- all_peaks%in%rownames(adaptive_DAP_CMVpos3)
peaks_CMVpos4 <- all_peaks%in%rownames(adaptive_DAP_CMVpos4)

DAP_relationship <- data.frame("CMVpos1" = peaks_CMVpos1,
                               "CMVpos2" = peaks_CMVpos2,
                               "CMVpos3" = peaks_CMVpos3,
                               "CMVpos4" = peaks_CMVpos4)

rownames(DAP_relationship) <- all_peaks


# Plot as an Euler diagram
fit1 <- euler(DAP_relationship)

plot(fit1,
     quantities = list(type = c("counts", "percent")),
     fills = c(adaptivepal[1], adaptivepal[1], adaptivepal[1], adaptivepal[1]),
     lty = 1:3,
     labels = list(font = 4))


# Add annotation to table
temp <- ClosestFeature(CMVpos2_full, regions = rownames(DAP_relationship))
DAP_relationship$gene_name <- temp$gene_name



### Bright vs Dim ###
# Find DAPs between adaptive and CD56dim
bright_dim_DAP_CMVpos1 <- FindMarkers(
  object = CMVpos1_full,
  ident.1 = c("CD56bright"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

bright_dim_DAP_CMVpos2 <- FindMarkers(
  object = CMVpos2_full,
  ident.1 = c("CD56bright"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

bright_dim_DAP_CMVpos3 <- FindMarkers(
  object = CMVpos3_full,
  ident.1 = c("CD56bright"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

bright_dim_DAP_CMVpos4 <- FindMarkers(
  object = CMVpos4_full,
  ident.1 = c("CD56bright"),
  ident.2 = c("CD56dim"),
  test.use = 'LR',
  only.pos = F,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1
)

# Filter for statistically significant hits
bright_dim_DAP_CMVpos2$gene <- rownames(bright_dim_DAP_CMVpos2)
bright_dim_DAP_CMVpos3$gene <- rownames(bright_dim_DAP_CMVpos3)
bright_dim_DAP_CMVpos4$gene <- rownames(bright_dim_DAP_CMVpos4)
bright_dim_DAP_CMVpos1$gene <- rownames(bright_dim_DAP_CMVpos1)

bright_dim_DAP_CMVpos2 <- bright_dim_DAP_CMVpos2  %>% dplyr::filter(p_val_adj < 0.05)
bright_dim_DAP_CMVpos3 <- bright_dim_DAP_CMVpos3  %>% dplyr::filter(p_val_adj < 0.05)
bright_dim_DAP_CMVpos4 <- bright_dim_DAP_CMVpos4  %>% dplyr::filter(p_val_adj < 0.05)
bright_dim_DAP_CMVpos1 <- bright_dim_DAP_CMVpos1  %>% dplyr::filter(p_val_adj < 0.05)

# Generate Euler Plot
# Find relationships

bright_dim_all_peaks <- unique(c(rownames(bright_dim_DAP_CMVpos2), rownames(bright_dim_DAP_CMVpos1),
                                 rownames(bright_dim_DAP_CMVpos3), rownames(bright_dim_DAP_CMVpos4)))

# Create matrix of logicals representing whether peak is a DAP for each dataset

bright_dim_peaks_CMVpos2 <- bright_dim_all_peaks%in%rownames(bright_dim_DAP_CMVpos2)
bright_dim_peaks_CMVpos3 <- bright_dim_all_peaks%in%rownames(bright_dim_DAP_CMVpos3)
bright_dim_peaks_CMVpos4 <- bright_dim_all_peaks%in%rownames(bright_dim_DAP_CMVpos4)
bright_dim_peaks_CMVpos1 <- bright_dim_all_peaks%in%rownames(bright_dim_DAP_CMVpos1)

DAP_relationship_bright_dim <- data.frame("CMVpos2" = bright_dim_peaks_CMVpos2,
                                          "CMvpos3" = bright_dim_peaks_CMVpos3,
                                          "CMVpos4" = bright_dim_peaks_CMVpos4,
                                          "CMVpos1" = bright_dim_peaks_CMVpos1)

rownames(DAP_relationship_bright_dim) <- bright_dim_all_peaks 

# Plot as an Euler diagram
fit2 <- euler(DAP_relationship_bright_dim)

plot(fit2,
     quantities = list(type = c("counts", "percent")),
     fills = c(clusterpal[2], clusterpal[2], clusterpal[2], clusterpal[2]),
     lty = 1:3,
     labels = list(font = 4))



##### Look at adaptive subcluster-specific peaks #####

# Find DAPs between adaptive
sub_adaptive_DAP_CMVpos2 <- FindAllMarkers(
  object = subset(CMVpos2_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
  test.use = 'LR',
  only.pos = T,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1,
  min.pct = 0.05
)

sub_adaptive_DAP_CMVpos3 <- FindAllMarkers(
  object = subset(CMVpos3_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
  test.use = 'LR',
  only.pos = T,
  latent.vars = 'nCount_MACS2',
  min.pct = 0.05,
  logfc.threshold = 0.1
)

sub_adaptive_DAP_CMVpos4 <- FindAllMarkers(
  object = subset(CMVpos4_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
  test.use = 'LR',
  only.pos = T,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1,
  min.pct = 0.05
)

sub_adaptive_DAP_CMVpos1 <- FindAllMarkers(
  object = subset(CMVpos1_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3",
                                           "Adaptive4", "Adaptive5")),
  test.use = 'LR',
  only.pos = T,
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1,
  min.pct = 0.05
)

# For CMVpos3, cluster 1 and 2 seem to be rather similar, include also DAPs comparing these two with cluster 3, otherwise
# those peaks are excluded because they are shared between 1 and 2
sub_adaptive_DAP2_CMVpos3 <- FindMarkers(
  object = CMVpos3_full,
  ident.1 = c("Adaptive1", "Adaptive2"),
  ident.2 = "Adaptive3",
  test.use = 'LR',
  only.pos = T,
  latent.vars = 'nCount_MACS2',
  min.pct = 0.05,
  logfc.threshold = 0.1
)
sub_adaptive_DAP2_CMVpos3$gene <- rownames(sub_adaptive_DAP2_CMVpos3)

# Filter for statistically significant hits
sub_adaptive_DAP_CMVpos2 <- sub_adaptive_DAP_CMVpos2  %>% dplyr::filter(p_val_adj < 0.05)
sub_adaptive_DAP_CMVpos3 <- sub_adaptive_DAP_CMVpos3  %>% dplyr::filter(p_val_adj < 0.05)
sub_adaptive_DAP2_CMVpos3 <- sub_adaptive_DAP2_CMVpos3  %>% dplyr::filter(p_val_adj < 0.05)
sub_adaptive_DAP_CMVpos4 <- sub_adaptive_DAP_CMVpos4  %>% dplyr::filter(p_val_adj < 0.05)
sub_adaptive_DAP_CMVpos1 <- sub_adaptive_DAP_CMVpos1  %>% dplyr::filter(p_val_adj < 0.05)

# Combine peaks for CMVpos3
sub_adaptive_DAP2_CMVpos3$cluster <- "Adaptive1"
rownames(sub_adaptive_DAP2_CMVpos3) <- NULL
sub_adaptive_DAP_CMVpos3_comb <- rbind(sub_adaptive_DAP_CMVpos3, sub_adaptive_DAP2_CMVpos3)


# Generate Euler Plot
# Find relationships

sub_all_peaks <- unique(c(sub_adaptive_DAP_CMVpos2$gene, sub_adaptive_DAP_CMVpos3_comb$gene,
                          sub_adaptive_DAP_CMVpos4$gene, sub_adaptive_DAP_CMVpos1$gene))

# Create matrix of logicals representing whether peak is a DAP for each dataset

sub_CMVpos2 <- sub_all_peaks%in%sub_adaptive_DAP_CMVpos2$gene
sub_CMVpos3 <- sub_all_peaks%in%sub_adaptive_DAP_CMVpos3_comb$gene
sub_CMVpos4 <- sub_all_peaks%in%sub_adaptive_DAP_CMVpos4$gene
sub_CMVpos1 <- sub_all_peaks%in%sub_adaptive_DAP_CMVpos1$gene

subDAP_relationship <- data.frame("CMVpos2" = sub_CMVpos2,
                                  "CMVpos3" = sub_CMVpos3,
                                  "CMVpos4" = sub_CMVpos4,
                                  "CMVpos1" = sub_CMVpos1)

rownames(subDAP_relationship) <- sub_all_peaks

# Plot as an Euler diagram
fit3 <- euler(subDAP_relationship)

plot(fit3,
     quantities = list(type = c("counts", "percent")),
     fills = c(adaptivepal[2], adaptivepal[3], adaptivepal[4], adaptivepal[5]),
     lty = 1:3,
     labels = list(font = 4 ))

# Add annotation to table
temp <- ClosestFeature(CMVpos3_full, regions = rownames(subDAP_relationship))
subDAP_relationship$gene_name <- temp$gene_name

# Look for some representative, donor specific, subcluster specific peaks
CMVpos2_specDAP <- subDAP_relationship %>% dplyr::filter(CMVpos2 == TRUE,
                                                         CMVpos3 == FALSE,
                                                         CMVpos4 == FALSE,
                                                         CMVpos1 == FALSE)
CMVpos2_specDAP$gene <- rownames(CMVpos2_specDAP)
CMVpos2_specDAP <- merge(CMVpos2_specDAP, sub_adaptive_DAP_CMVpos2, by = "gene")

CMVpos3_specDAP <- subDAP_relationship %>% dplyr::filter(CMVpos2 == FALSE,
                                                         CMVpos3 == TRUE,
                                                         CMVpos4 == FALSE,
                                                         CMVpos1 == FALSE)
CMVpos3_specDAP$gene <- rownames(CMVpos3_specDAP)
CMVpos3_specDAP <- merge(CMVpos3_specDAP, sub_adaptive_DAP_CMVpos3, by = "gene")

CMVpos4_specDAP <- subDAP_relationship %>% dplyr::filter(CMVpos2 == FALSE,
                                                         CMVpos3 == FALSE,
                                                         CMVpos4 == TRUE,
                                                         CMVpos1 == FALSE)
CMVpos4_specDAP$gene <- rownames(CMVpos4_specDAP)
CMVpos4_specDAP <- merge(CMVpos4_specDAP, sub_adaptive_DAP_CMVpos4, by = "gene")

CMVpos1_specDAP <- subDAP_relationship %>% dplyr::filter(CMVpos2 == FALSE,
                                                         CMVpos3 == FALSE,
                                                         CMVpos4 == FALSE,
                                                         CMVpos1 == TRUE)
CMVpos1_specDAP$gene <- rownames(CMVpos1_specDAP)
CMVpos1_specDAP <- merge(CMVpos1_specDAP, sub_adaptive_DAP_CMVpos1, by = "gene")



##### Pseudobulk_analysis to identify convergent and divergent aspects of
##### adaptive subcluster across donors ####


##### CD56 dim vs Adaptive #####
DefaultAssay(CMVpos2_full) <- "MACS2"
DefaultAssay(CMVpos3_full) <- "MACS2"
DefaultAssay(CMVpos4_full) <- "MACS2"
DefaultAssay(CMVpos1_full) <- "MACS2"

# Remove all assays that are not required
stripassay <- function(seu){
  seu[["mito"]] <- NULL
  seu[["alleles"]] <- NULL
  seu[["HTO"]] <- NULL
  seu[["chromvar_names"]] <- NULL
  seu[["RNA_imputed"]] <- NULL
  return(seu)
}

CMVpos2_full <- stripassay(CMVpos2_full)
CMVpos3_full <- stripassay(CMVpos3_full)
CMVpos4_full <- stripassay(CMVpos4_full)
CMVpos1_full <- stripassay(CMVpos1_full)


# Merge all donors
seu_dim_adaptive <- merge(subset(CMVpos2_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3", "CD56dim")),
                          c(subset(CMVpos3_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3","CD56dim")),
                            subset(CMVpos4_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3","CD56dim")),
                            subset(CMVpos1_full, idents = c("Adaptive1", "Adaptive2","Adaptive3",
                                                       "Adaptive4", "Adaptive5", "CD56dim"))
                            )
)
# Normalize
seu_dim_adaptive <- RunTFIDF(seu_dim_adaptive)

seu_dim_adaptive$annotation2 <- paste0(seu_dim_adaptive$annotation, seu_dim_adaptive$donor)
Idents(seu_dim_adaptive) <- seu_dim_adaptive$annotation2



adaptive_peaks <- unique(c(adaptive_DAP_CMVpos2 %>% dplyr::filter (avg_log2FC <0) %>% arrange(desc(avg_log2FC)) %>% pull(gene),
                           adaptive_DAP_CMVpos3 %>% dplyr::filter (avg_log2FC <0) %>% arrange(desc(avg_log2FC)) %>% pull(gene),
                           adaptive_DAP_CMVpos4 %>% dplyr::filter (avg_log2FC <0) %>% arrange(desc(avg_log2FC)) %>% pull(gene),
                           adaptive_DAP_CMVpos1 %>% dplyr::filter (avg_log2FC <0) %>% arrange(desc(avg_log2FC)) %>% pull(gene),
                           adaptive_DAP_CMVpos2 %>% dplyr::filter (avg_log2FC >0) %>% arrange(avg_log2FC) %>% pull(gene),
                           adaptive_DAP_CMVpos3 %>% dplyr::filter (avg_log2FC >0) %>% arrange(avg_log2FC) %>% pull(gene),
                           adaptive_DAP_CMVpos4 %>% dplyr::filter (avg_log2FC >0) %>% arrange(avg_log2FC) %>% pull(gene),
                           adaptive_DAP_CMVpos1 %>% dplyr::filter (avg_log2FC >0) %>% arrange(avg_log2FC) %>% pull(gene)
)
)

# Keep only those shared between at least three donors
adaptive_shared_peaks <- adaptive_peaks[adaptive_peaks%in%rownames(DAP_relationship[rowSums(DAP_relationship)>2,])]

average_acc_adaptive <- AverageExpression(seu_dim_adaptive, features = adaptive_shared_peaks, assays = "MACS2", slot = "data")

average_acc_adaptive <- as.data.frame(average_acc_adaptive[["MACS2"]])

average_acc_adaptive$region <- rownames(average_acc_adaptive)

average_acc_tidy_adaptive <- gather(average_acc_adaptive, key = "Cluster" , value = "Normalized_Accessibility", -region)

average_acc_tidy_adaptive %>%  tidy_heatmap(
  rows = Cluster,
  columns = region,
  values = Normalized_Accessibility,
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  scale = "column",
  colors =  rev(mapal)
)


### Comparison of subclusters ####
rm(seu_dim_adaptive)

seu_adaptive <- merge(subset(CMVpos2_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
                      c(subset(CMVpos3_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
                      subset(CMVpos4_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
                      subset(CMVpos1_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3",
                                                      "Adaptive4", "Adaptive5")))
)
# Normalize
seu_adaptive <- RunTFIDF(seu_adaptive)

seu_adaptive$annotation2 <- paste0(seu_adaptive$annotation, seu_adaptive$donor)

Idents(seu_adaptive) <- seu_adaptive$annotation2

sub_adaptive_DAP_CMVpos2 <- sub_adaptive_DAP_CMVpos2 %>% group_by(cluster) %>% arrange(p_val_adj, .by_group = T)
sub_adaptive_DAP_CMVpos3_comb <- sub_adaptive_DAP_CMVpos3_comb %>% group_by(cluster) %>% arrange(p_val_adj, .by_group = T)
sub_adaptive_DAP_CMVpos4 <- sub_adaptive_DAP_CMVpos4 %>% group_by(cluster) %>%  arrange(p_val_adj, .by_group = T)
sub_adaptive_DAP_CMVpos1 <- sub_adaptive_DAP_CMVpos1 %>% group_by(cluster) %>%  arrange(p_val_adj, .by_group = T)




average_acc <- AverageExpression(seu_adaptive, features = unique(c(sub_adaptive_DAP_CMVpos1$gene,
                                                                   sub_adaptive_DAP_CMVpos2$gene,
                                                                   sub_adaptive_DAP_CMVpos4$gene,
                                                                   sub_adaptive_DAP_CMVpos3_comb$gene
)), assays = "MACS2", slot = "data")

 
average_acc <- as.data.frame(average_acc[["MACS2"]])

average_acc$region <- rownames(average_acc)

average_acc_tidy <- gather(average_acc, key = "Cluster" , value = "Normalized_Accessibility", -region)

average_acc_tidy %>%  tidy_heatmap(
  rows = Cluster,
  columns = region,
  values = Normalized_Accessibility,
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  scale = "column",
  colors =  cool_warm(500)
)

# Reorder  peaks so that they follow the hierarchical clustering of the subclusters

peaks_ordered <- c(sub_adaptive_DAP_CMVpos2 %>% dplyr::filter(cluster == "Adaptive1") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos4 %>% dplyr::filter(cluster == "Adaptive1") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos4 %>% dplyr::filter(cluster == "Adaptive3") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos1 %>% dplyr::filter(cluster == "Adaptive1") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos1 %>% dplyr::filter(cluster == "Adaptive2") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos1 %>% dplyr::filter(cluster == "Adaptive3") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos1 %>% dplyr::filter(cluster == "Adaptive4") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos1 %>% dplyr::filter(cluster == "Adaptive5") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos3 %>% dplyr::filter(cluster == "Adaptive3") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos3 %>% dplyr::filter(cluster == "Adaptive1") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos3 %>% dplyr::filter(cluster == "Adaptive2") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos4 %>% dplyr::filter(cluster == "Adaptive2") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos2 %>% dplyr::filter(cluster == "Adaptive2") %>% pull(gene),
                   sub_adaptive_DAP_CMVpos2 %>% dplyr::filter(cluster == "Adaptive3") %>% pull(gene)
)

# Keep only unique ones
peaks_ordered <- unique(peaks_ordered)
peaks_ordered <- peaks_ordered[peaks_ordered%in%unique(c(sub_adaptive_DAP_CMVpos2$gene,sub_adaptive_DAP_CMVpos3_comb$gene,
                                                         sub_adaptive_DAP_CMVpos4$gene, sub_adaptive_DAP_CMVpos1$gene))]

# Get average expression again in correct order
average_acc <- AverageExpression(seu_adaptive, features = peaks_ordered, assays = "MACS2", slot = "data")
average_acc <- as.data.frame(average_acc[["MACS2"]])

average_acc$region <- rownames(average_acc)

# Subcluster-specific peaks of Adaptive1 CMVpos2 are buried a bit in shared signature with 1/3 of CMVpos4 => reorder
# Find specific peaks
a1c2_peaks <- FindMarkers(seu_adaptive, ident.1 = "Adaptive1CMVpos2",
                          ident.2 = c("Adaptive1CMVpos4", "Adaptive3CMVpos4"),
                          test.use = 'LR',
                          latent.vars = 'nCount_MACS2',
                          logfc.threshold = 0.1,
                          min.pct = 0.05,
                          only.pos = T)

a1c2_peaks$region <- rownames(a1c2_peaks)

# Keep only those present in original, within-donor analysis
a1c2_peaks <- a1c2_peaks[a1c2_peaks$region%in%sub_adaptive_DAP_CMVpos2$gene&!a1c2_peaks$region%in%sub_adaptive_DAP_CMVpos4$gene,]

# Plot first shared peaks of group 1, then subcluster specific peaks
peaks_reordered <- unique(c(average_acc$region[average_acc$region%in%a1c2_peaks$region],
                            sub_adaptive_DAP_CMVpos4 %>% dplyr::filter(cluster == "Adaptive1") %>% pull(gene),
                     average_acc$region
                     ))



average_acc_ordered <- AverageExpression(seu_adaptive, features = peaks_reordered, assays = "MACS2", slot = "data")
average_acc_ordered <- as.data.frame(average_acc_ordered[["MACS2"]])

average_acc_ordered$region <- rownames(average_acc_ordered)

# Alt1
average_acc_tidy <- gather(average_acc_ordered, key = "Cluster" , value = "Normalized_Accessibility", -region)
# Alt2 (w/o re-ordering)
average_acc_tidy <- gather(average_acc, key = "Cluster" , value = "Normalized_Accessibility", -region)

average_acc_tidy %>%  tidy_heatmap(
  rows = Cluster,
  columns = region,
  values = Normalized_Accessibility,
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  scale = "column",
  colors =  rev(mapal)
)


#### Find markers between subcluster groups ####


seu_adaptive <- NormalizeData(
  object = seu_adaptive,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(seu_adaptive$nCount_RNA_atac)
)

# Find differential gene scores between the two groups
DefaultAssay(seu_adaptive) <- "RNA_atac"
group_scores <- FindMarkers(
  object = seu_adaptive,
  ident.1 = c("Adaptive3CMVpos4", "Adaptive1CMVpos2", "Adaptive1CMVpos4"),
  ident.2 = c("Adaptive2CMVpos2", "Adaptive2CMVpos4",
              "Adaptive3CMVpos2", "Adaptive1CMVpos3", "Adaptive2CMVpos3",
              "Adaptive1CMVpos1", "Adaptive2CMVpos1", "Adaptive3CMVpos1",
              "Adaptive4CMVpos1", "Adaptive5CMVpos1","Adaptive3CMVpos3"),
  logfc.threshold = 0.1
)

group_scores$gene <- rownames(group_scores)
group_scores <- group_scores %>% dplyr::filter(p_val_adj < 0.05)
group_genes <- c(group_scores %>% top_n(n= 100, wt = avg_log2FC) %>% arrange(desc(avg_log2FC)) %>% pull(gene),
                       group_scores %>% top_n(n= 100, wt = -avg_log2FC) %>% arrange(avg_log2FC) %>% pull(gene))

average_act <- AverageExpression(seu_adaptive, features = group_genes, assays = "RNA_atac", slot = "data")

average_act <- as.data.frame(average_act[["RNA_atac"]])


average_act$gene <- rownames(average_act)

average_act_tidy <- gather(average_act, key = "Cluster" , value = "Activity", -gene)

average_act_tidy %>%  tidy_heatmap(
  rows = Cluster,
  columns = gene,
  values = Activity,
  cluster_rows = TRUE,
  scale = "column",
  colors =  blueyellow,
  fontsize_col = 2
)


# Peaks
DefaultAssay(seu_adaptive) <- "MACS2"
group_peaks <- FindMarkers(
  object = seu_adaptive,
  ident.1 = c("Adaptive3CMVpos4", "Adaptive1CMVpos2", "Adaptive1CMVpos4"),
  ident.2 = c("Adaptive3CMVpos3", "Adaptive2CMVpos2", "Adaptive2CMVpos4",
              "Adaptive3CMVpos2", "Adaptive1CMVpos3", "Adaptive2CMVpos3",
              "Adaptive1CMVpos1", "Adaptive2CMVpos1", "Adaptive3CMVpos1",
              "Adaptive4CMVpos1", "Adaptive5CMVpos1"),
  test.use = 'LR',
  latent.vars = 'nCount_MACS2',
  logfc.threshold = 0.1,
  min.pct = 0.05
)

group_peaks$region <- rownames(group_peaks)

group_peaks <- group_peaks %>% dplyr::filter(p_val_adj < 0.05)
group_peaks <- group_peaks %>% arrange(avg_log2FC)

average_acc_groups <- AverageExpression(seu_adaptive, features = group_peaks$region, assays = "MACS2", slot = "data")

average_acc_groups <- as.data.frame(average_acc_groups[["MACS2"]])

average_acc_groups$region <- rownames(average_acc_groups)

average_acc_groups_tidy <- gather(average_acc_groups, key = "Cluster" , value = "Normalized_Accessibility", -region)

average_acc_groups_tidy %>%  tidy_heatmap(
  rows = Cluster,
  columns = region,
  values = Normalized_Accessibility,
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  scale = "column",
  colors =  rev(mapal)
)

temp <- ClosestFeature(CMVpos2_full, rownames(group_peaks))
group_peaks$closest_gene <- temp$gene_name


##### Comparison to conventional subsets #####
integrated <- readRDS(paste0(path, "integrated_processed"))


# Plot signature of Adaptive1 vs Adaptive2 for conventional subsets
DefaultAssay(integrated) <- "MACS2"
average_acc_bright_dim <- AverageExpression(integrated, features = group_peaks$region, assays = "MACS2")

average_acc_bright_dim <- as.data.frame(average_acc_bright_dim[["MACS2"]])

average_acc_bright_dim$region <- rownames(average_acc_bright_dim)

average_acc_bright_dim <- gather(average_acc_bright_dim, key = "Cluster" , value = "Normalized_Accessibility", -region)

average_acc_bright_dim$Cluster <- factor(average_acc_bright_dim$Cluster, levels = c("CD56dim", "EarlyCD56dim", "CD56bright", "Adaptive"))
average_acc_bright_dim %>% dplyr::filter(Cluster%in%c("CD56bright", "CD56dim", "EarlyCD56dim")) %>%  tidy_heatmap(
  rows = Cluster,
  columns = region,
  values = Normalized_Accessibility,
  #cluster_rows = TRUE,
  #cluster_cols = TRUE,
  scale = "column",
  colors =  rev(mapal)
)


# Get foldchanges of significant peaks
peakFC_earlydim <- FoldChange(integrated, ident.1 = "CD56dim", ident.2 = "EarlyCD56dim",
                               features = group_peaks$region)
peakFC_bright <- FoldChange(integrated, ident.1 = "CD56dim", ident.2 = "CD56bright",
                              features = group_peaks$region)
peakFC_adaptive <- FoldChange(seu_adaptive, ident.1 = c("Adaptive3CMVpos3", "Adaptive2CMVpos2", "Adaptive2CMVpos4",
                                                        "Adaptive3CMVpos2", "Adaptive1CMVpos3", "Adaptive2CMVpos3",
                                                        "Adaptive1CMVpos1", "Adaptive2CMVpos1", "Adaptive3CMVpos1",
                                                        "Adaptive4CMVpos1", "Adaptive5CMVpos1"),
                         ident.2 =c("Adaptive3CMVpos4", "Adaptive1CMVpos2", "Adaptive1CMVpos4"),
                         features = group_peaks$region)

# Put into one dataframe for plotting
peakFC_earlydim$region <- rownames(peakFC_earlydim)
peakFC_bright$region <- rownames(peakFC_bright)
peakFC_adaptive$region <- rownames(peakFC_adaptive)
colnames(peakFC_earlydim)[1] <- "log2FC_CD56dimvsEarlyCD56dim"
colnames(peakFC_bright)[1] <- "log2FC_CD56dimvsCD56bright"
colnames(peakFC_adaptive)[1] <- "log2FC_Adaptive1_vs_Adaptive2"
peakFC_earlydim <- peakFC_earlydim[,c(1,4)]
peakFC_adaptive <- peakFC_adaptive[,c(1,4)]
peakFC_bright <- peakFC_bright[,c(1,4)]

combined_FC_early <- merge(peakFC_earlydim, peakFC_adaptive, by = "region", all = TRUE)
# or for bright
combined_FC_bright <- merge(peakFC_bright, peakFC_adaptive, by = "region", all = TRUE)

# Add closest gene
temp <- ClosestFeature(integrated, regions = combined_FC_early$region)
combined_FC_early$closest_gene <- temp$gene_name

temp <- ClosestFeature(integrated, regions = combined_FC_bright$region)
combined_FC_bright$closest_gene <- temp$gene_name


genes_to_label <- c("PCNT", "TCF7", "GLI3", "IL7R", "SLC9A9", "BACH2", "RUNX2")
combined_FC_early%>%
  ggplot(aes(x = log2FC_Adaptive1_vs_Adaptive2, y = log2FC_CD56dimvsEarlyCD56dim))+
  geom_point(size = 1)+
  scale_color_manual(values = c("black", "darkred"))+
  coord_cartesian()+
  geom_text_repel(aes(x = log2FC_Adaptive1_vs_Adaptive2, y = log2FC_CD56dimvsEarlyCD56dim, label = closest_gene),
                   data = combined_FC_early[combined_FC_early$closest_gene%in%genes_to_label,],
                   size = 4, max.overlaps = Inf, box.padding = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_smooth(method = "lm")+
  xlab("Adaptive1vsAdaptive2 (Log2FC)")+
  ylab("CD56dim vs EarlyCD56dim (Log2FC)")+
  #scale_y_continuous(limits = c(-0.6, 0.5))+
  #scale_x_continuous(limits = c(-0.6,0.5))+
  theme_classic()

combined_FC_bright%>%
  ggplot(aes(x = log2FC_Adaptive1_vs_Adaptive2, y = log2FC_CD56dimvsCD56bright))+
  geom_point(size = 1)+
  scale_color_manual(values = c("black", "darkred"))+
  coord_cartesian()+
  geom_text_repel(aes(x = log2FC_Adaptive1_vs_Adaptive2, y = log2FC_CD56dimvsCD56bright, label = closest_gene),
                  data = combined_FC_bright[combined_FC_bright$closest_gene%in%genes_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_smooth(method = "lm")+
  xlab("Adaptive1vsAdaptive2 (Log2FC)")+
  ylab("CD56dim vs CD56bright(Log2FC)")+
  #scale_y_continuous(limits = c(-0.6, 0.5))+
  #scale_x_continuous(limits = c(-0.6,0.5))+
  theme_classic()


lmfull <- lm(log2FC_Adaptive1_vs_Adaptive2~log2FC_CD56dimvsEarlyCD56dim, data = combined_FC_early) #Create the linear regression
summary(lmfull) #Review the results

lmfull <- lm(log2FC_Adaptive1_vs_Adaptive2~log2FC_CD56dimvsCD56bright, data = combined_FC_bright) #Create the linear regression
summary(lmfull) #Review the results





# Plot representatives
levels(seu_adaptive) <- c("Adaptive1CMVpos2", "Adaptive1CMVpos4","Adaptive3CMVpos4",
                          "Adaptive1CMVpos1","Adaptive2CMVpos1","Adaptive3CMVpos1",
                          "Adaptive4CMVpos1","Adaptive5CMVpos1","Adaptive3CMVpos3",
                          "Adaptive1CMVpos3","Adaptive2CMVpos3","Adaptive2CMVpos4",
                          "Adaptive2CMVpos2","Adaptive3CMVpos2")
                          
                          


# Shared peaks
# Aim2 "chr1-159076355-159077686"
CoveragePlot(
  object = seu_adaptive,
  region = "chr1-159076355-159077686",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,scale.factor = 1e7
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))

# Cadm1 "chr11-115222052-115223633"
CoveragePlot(
  object = seu_adaptive,
  region = "chr11-115222052-115223633",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,scale.factor = 1e7
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))




# PCNT chr21-46423829-46426010
CoveragePlot(
  object = seu_adaptive,
  region = "chr21-46423829-46426010",
  extend.upstream = 5000,
  extend.downstream = 5000,
  links = F,scale.factor = 1e7
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))


# TCF7 chr5-134132043-134132878
CoveragePlot(
  object = seu_adaptive,
  region = "chr5-134132043-134132878",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,scale.factor = 1e7
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))


# Subcluster-specific peaks
# CMVpos2
CoveragePlot(
  object = seu_adaptive,
  region = "chr21-31087852-31088428",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))
# CMVpos4
CoveragePlot(
  object = seu_adaptive,
  region = "chr17-67398748-67399002",
  extend.upstream = 4000,
  extend.downstream = 4000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))
# CMVpos3
CoveragePlot(
  object = seu_adaptive,
  region = "chr16-81678354-81679431",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))
# CMVpos1
CoveragePlot(
  object = seu_adaptive,
  region = "chr16-88013040-88013660",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(adaptivepal, adaptivepal))



# Plot Surface phenotypes
DefaultAssay(CMVpos2_full) <- "ADT"
CMVpos2_plot <- sample(WhichCells(CMVpos2_full, expression = experiment == "mtASAP3"), size = 2000)
DefaultAssay(CMVpos4_full) <- "ADT"
CMVpos4_plot <- sample(Cells(CMVpos4_full), size = 2000)
DefaultAssay(CMVpos3_full) <- "ADT"
CMVpos3_plot <- sample(Cells(CMVpos3_full), size = 2000)
DefaultAssay(CMVpos1_full) <- "ADT"
CMVpos1_plot <- sample(Cells(CMVpos1_full), size = 2000)

pdf(file = "Plots/new/SupplFig5_ADT.pdf")
FeaturePlot(CMVpos2_full, cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", coord.fixed = T,
            ncol = 3, cells = CMVpos2_plot, pt.size = 0.1, order = TRUE,
            features = c("NKp30", "anti-PE", "anti-Biotin",
                         "KIR2DL1-S1-S3-S5", "KIR3DL1", "KIR2DL2-L3"))&NoLegend()&NoAxes()
FeaturePlot(CMVpos4_full, cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", coord.fixed = T,
            ncol = 3,  cells = CMVpos4_plot, pt.size = 0.1, order = TRUE,
            features = c("NKp30", "anti-PE", "anti-Biotin",
                         "KIR2DL1-S1-S3-S5", "KIR3DL1", "KIR2DL2-L3"))&NoLegend()&NoAxes()
FeaturePlot(CMVpos3_full, cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", coord.fixed = T,
            ncol = 3,  cells = CMVpos3_plot, pt.size = 0.1, order = TRUE,
            features = c("NKp30", "Anti-PE", "anti-Biotin",
                         "KIR2DL1-S1-S3-S5", "KIR3DL1", "KIR2DL2-L3"))&NoLegend()&NoAxes()
FeaturePlot(CMVpos1_full, cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", coord.fixed = T,
            ncol = 3,  cells = CMVpos1_plot, pt.size = 0.1, order = TRUE,
            features = c("NKp30", "Anti-PE", "anti-Biotin",
                         "KIR2DL1-S1-S3-S5", "KIR3DL1", "KIR2DL2-L3"))&NoLegend()&NoAxes()
dev.off()

# CD62L expression of adaptive subcluster groups
seu_adaptive_ASAP <- subset(seu_adaptive, subset = experiment%in%c("mtASAP2", "mtASAP3", "mtASAP1", "mtASAP5"))
VlnPlot(seu_adaptive_ASAP, features = c("CD62L"), cols = c(adaptivepal, adaptivepal), pt.size = 0, y.max = 2.5)&NoLegend()
DefaultAssay(seu_adaptive_ASAP) <- "ADT"

CD62_test <- FindMarkers(seu_adaptive_ASAP,
                         ident.1 = c("Adaptive3CMVpos4", "Adaptive1CMVpos2", "Adaptive1CMVpos4"),
                         ident.2 = c("Adaptive3CMVpos3", "Adaptive2CMVpos2", "Adaptive2CMVpos4",
                                     "Adaptive3CMVpos2", "Adaptive1CMVpos3", "Adaptive2CMVpos3",
                                     "Adaptive1CMVpos1", "Adaptive2CMVpos1", "Adaptive3CMVpos1",
                                     "Adaptive4CMVpos1", "Adaptive5CMVpos1"),
                         features = "CD62L")



### Differential motif activity between groups ###

# Add assay with motif names rather than IDs for more convenient plotting
# Get residuals for merged object
motif.matrix <- CreateMotifMatrix(
  features = granges(seu_adaptive),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)

# Create a new Motif object to store the results
motif_obj <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
seu_adaptive <- SetAssayData(
  object = seu_adaptive,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif_obj
)


# Calculate regional stats
seu_adaptive <- RegionStats(object = seu_adaptive, genome = BSgenome.Hsapiens.UCSC.hg38)

seu_adaptive <- RunChromVAR(
  object = seu_adaptive,
  genome = BSgenome.Hsapiens.UCSC.hg38)

chromvar_names <- GetAssayData(seu_adaptive, assay = "chromvar", slot = "data")
rownames(chromvar_names) <- ConvertMotifID(seu_adaptive, id = rownames(seu_adaptive[["chromvar"]]), assay = "MACS2")
seu_adaptive[["chromvar_names"]] <- CreateAssayObject(data = chromvar_names)

DefaultAssay(seu_adaptive) <- "chromvar_names"

group_motifacts <- FindMarkers(
  object = seu_adaptive,
  ident.1 = c("Adaptive3CMVpos4", "Adaptive1CMVpos2", "Adaptive1CMVpos4"),
  ident.2 = c("Adaptive3CMVpos3", "Adaptive2CMVpos2", "Adaptive2CMVpos4",
              "Adaptive3CMVpos2", "Adaptive1CMVpos3", "Adaptive2CMVpos3",
              "Adaptive1CMVpos1", "Adaptive2CMVpos1", "Adaptive3CMVpos1",
              "Adaptive4CMVpos1", "Adaptive5CMVpos1"),
  logfc.threshold = 0.1,
  min.pct = 0.05
)

VlnPlot(seu_adaptive, features = c("TCF7L2"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()&scale_y_continuous(limits = c(-3,5))
VlnPlot(seu_adaptive, features = c("LEF1"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()&scale_y_continuous(limits = c(-4,5))
VlnPlot(seu_adaptive, features = c("CTCF"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()&scale_y_continuous(limits = c(-4,5))

VlnPlot(seu_adaptive, features = c("ZEB1"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()
VlnPlot(seu_adaptive, features = c("RUNX2"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()&scale_y_continuous(limits = c(-4,5))
VlnPlot(seu_adaptive, features = c("FOS::JUNB"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()
VlnPlot(seu_adaptive, features = c("LEF1"), cols = c(adaptivepal, adaptivepal), pt.size = 0)&NoLegend()&scale_y_continuous(limits = c(-4,5))
