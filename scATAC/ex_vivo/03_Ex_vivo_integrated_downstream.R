# Downstream analysis of integrated ex vivo dataset:
# Peak, gene score, motif analyses; Separated analysis of CMV+/-

library(Seurat)
library(Signac)
library(ggplot2)
library(scico)
library(heatmaply)
library(patchwork)
library(ggrepel)
library(MetBrewer)
library(RColorBrewer)
library(dplyr)

path <- "~/"

set.seed(1234)

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- met.brewer(name="Isfahan1",n=7,type="discrete")
clusterpal <- clusterpal[c(4,5,7,1)]
NKG2Cpal <- c("#40BAD5","#120136")
adaptivepal <- c("#880E4F",  "#C2185B",  "#E91E63",  "#F06292", "#F8BBD0", "#4A148C", "#7B1FA2", "#9C27B0", "#BA68C8")
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
# From GreenleafLab/ArchR
blueyellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")


# Import data
integrated <- readRDS(paste0(path, "integrated_processed"))

##### Figure 1: Clustering and annotation using CITE-Seq #####
pdf(file = "Plots/Fig1_Dimplot_ATAC.pdf", paper = "a4")
DimPlot(integrated, cols = c(clusterpal[1:3], adaptivepal))
DimPlot(integrated, cols = c(clusterpal[1:3], adaptivepal), shuffle = T)&NoLegend()
dev.off()

DefaultAssay(integrated) <- "ADT"

# Plot only cells with CITE-Seq stainings
cells_plot <- WhichCells(integrated, expression = experiment%in%c("mtASAP1", "mtASAP2"))

# Downsample to avoid overplotting
cells_downsampled <- sample(cells_plot, size = 3000)

pdf(file = "Plots/Fig1_CITE_ATAC.pdf", paper = "a4")
FeaturePlot(
  object = integrated,
  features = c("CD56", "CD62L","CD16", "CD57", "NKp30", "anti-Biotin", "CD161", "CD2"),
  pt.size = 0.01,
  max.cutoff = 'q99',
  min.cutoff = 'q5',
  ncol = 4,
  cols = colorscale,
  cells = cells_downsampled,
  coord.fixed = TRUE
)&NoLegend()&NoAxes()
FeaturePlot(
  object = integrated,
  features = c("CD56", "CD62L","CD16", "CD57", "NKp30", "anti-Biotin", "CD161", "CD2"),
  pt.size = 0.01,
  max.cutoff = 'q99',
  min.cutoff = 'q5',
  ncol = 4,
  cols = colorscale,
  cells = cells_downsampled,
  coord.fixed = TRUE
)&NoAxes()
dev.off()


##### Gene activity markers #####

DefaultAssay(integrated) <- 'RNA_atac'

activity_markers <- FindAllMarkers(integrated, only.pos = TRUE)
CD56dim_act <- FindMarkers(integrated, ident.1 = c("CD56dim"),
                                              ident.2 = c("CD56bright"), only.pos = T)
CD56dim_act$gene <- rownames(CD56dim_act)

# Plot representative markers
levels(integrated) <- rev(levels(integrated))

pdf(file = "Plots/Fig1_ATAC_dotplot.pdf", paper = "a4r")
DotPlot(integrated, features = c(
                                 "TCF7", "RUNX2", "BACH2", "MAML3", "ZMAT4", "TCF4","MEF2C", "ZBTB16", "ZNF516", "GLI3",
                                   "BCL11B", "ZEB2","ARID5B", "ZBTB38", "JAKMIP1"),
        dot.min = 0.01)+scale_color_gradientn(colours = blueyellow)+theme(axis.text.x = element_text(angle = 45))
dev.off()


# Pseudo-bulk heatmap
library(tidyr)
library(tidyheatmap)

# Combine markers from FindAllMarkers and Bright vs Dim comparison as otherwise these are masked
# due to similarity between adaptive and dim
heatmap_acts <- unique(c(activity_markers[activity_markers$cluster%in%c("CD56bright", "EarlyCD56dim"),"gene"],
                         CD56dim_act$gene[!CD56dim_act$gene%in%activity_markers$gene],
                     activity_markers[activity_markers$cluster%in%c("CD56dim", "Adaptive"),"gene"]))

averaged_gene_act <- AverageExpression(integrated, features = heatmap_acts, slot = "data", assays = "RNA_atac")

averaged_gene_act <- as.data.frame(averaged_gene_act$RNA_atac)

averaged_gene_act$gene <- rownames(averaged_gene_act)

averaged_gene_act_tidy <- gather(averaged_gene_act, key = "group", value = "activity", -gene)

pdf(file = "Plots/Fig1_ATAC_heatmap_activity.pdf", paper = "a4")
averaged_gene_act_tidy  %>% tidy_heatmap(
  rows = gene,
  columns = group,
  values = activity,
  cluster_rows = F,
  cluster_cols = F,
  scale = "row",
  colors =  blueyellow)
dev.off()


# Pairwise comparisons
EarlyCD56dim_vs_bright_activity <- FindMarkers(integrated, ident.1 = c("EarlyCD56dim"),
                                      ident.2 = c("CD56bright"))


CD56dim_vs_EarlyCD56dim_activity <- FindMarkers(integrated, ident.1 = c("CD56dim"),
                                       ident.2 = c("EarlyCD56dim"))

CD56dim_vs_CD56bright_activity <- FindMarkers(integrated, ident.1 = c("CD56dim"),
                                     ident.2 = c("CD56bright"))

adaptive_vs_dim_activity <- FindMarkers(integrated, ident.1 = c("Adaptive"),
                               ident.2 = c("CD56dim"))

adaptive_vs_dim_activity <- adaptive_vs_dim_activity %>% filter(p_val_adj < 0.05)
CD56dim_vs_EarlyCD56dim_activity <- CD56dim_vs_EarlyCD56dim_activity %>% filter(p_val_adj < 0.05)
CD56dim_vs_CD56bright_activity <- CD56dim_vs_CD56bright_activity %>% filter(p_val_adj < 0.05)
EarlyCD56dim_vs_bright_activity <- EarlyCD56dim_vs_bright_activity %>% filter(p_val_adj < 0.05)


# Plots
adaptive_vs_dim_activity$gene <- rownames(adaptive_vs_dim_activity)
EarlyCD56dim_vs_bright_activity$gene <- rownames(EarlyCD56dim_vs_bright_activity)
CD56dim_vs_EarlyCD56dim_activity$gene <- rownames(CD56dim_vs_EarlyCD56dim_activity)
CD56dim_vs_CD56bright_activity$gene <- rownames(CD56dim_vs_CD56bright_activity)


pdf(file = "Plots/SupplFig1_ATAC_volcano.pdf")
genestolabel <- c("ZBTB38", "CIITA", "RPTOR", "ZEB2", "ZBTB16", "IKZF2", "SYK", "CADM1", "KLRC2", "AIM2", "FCER1G", "IL12RB2", "BCL11B", "ARID5B")
adaptive_vs_dim_activity %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304), label = gene),
                  data = adaptive_vs_dim_activity[adaptive_vs_dim_activity$gene%in%c(genestolabel),],
                  size = 5, box.padding = 0.7)+
  scale_y_continuous(limits = c(0, 305))+
  scale_x_continuous(limits = c(-1.3, 1.3))+
  xlab(label = "Adaptive vs CD56dim (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()


genestolabel <- c("TCF7", "RUNX2", "BACH2", "MAML3", "ZMAT4", "TCF4","MEF2C", "IKZF2","ZBTB16", "ZNF516", "GLI3",
                  "ZEB2",  "BCL11B", "ZBTB38", "IL7R")
EarlyCD56dim_vs_bright_activity %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304), label = gene),
                  data = EarlyCD56dim_vs_bright_activity[EarlyCD56dim_vs_bright_activity$gene%in%c(genestolabel),],
                  size = 5, box.padding = 0.5)+
  scale_y_continuous(limits = c(0, 305))+
  scale_x_continuous(limits = c(-1.3, 1.3))+
  xlab(label = "Early CD56dim vs CD56bright (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()


CD56dim_vs_EarlyCD56dim_activity %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304), label = gene),
                  data = CD56dim_vs_EarlyCD56dim_activity[CD56dim_vs_EarlyCD56dim_activity$gene%in%c(genestolabel),],
                  size = 5, box.padding = 0.7)+
  scale_y_continuous(limits = c(0, 305))+
  scale_x_continuous(limits = c(-1.3, 1.3))+
  xlab(label = "CD56dim vs Early CD56dim (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()



CD56dim_vs_CD56bright_activity %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-304), label = gene),
                  data = CD56dim_vs_CD56bright_activity[CD56dim_vs_CD56bright_activity$gene%in%c(genestolabel),],
                  size = 5, box.padding = 0.7)+
  scale_y_continuous(limits = c(0, 305))+
  scale_x_continuous(limits = c(-1.3, 1.3))+
  xlab(label = "CD56dim vs CD56bright (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()
dev.off()





# Add annotation of peaks as promoter/exon/intron/intergenic
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annots <- c('hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)
DefaultAssay(integrated) <- "MACS2"

# Global overview of peak distribution
all_features <- ClosestFeature(integrated, regions = rownames(integrated[["MACS2"]]), annotation = annotations)

# Counts types
all_features %>% count(type) -> Peak_type_numbers
# Show as pie chart
pdf(file = "SupplFig1_ATAC_peak_types.pdf")
Peak_type_numbers %>% 
  ggplot(aes(x = "", y = n, fill = type))+
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette = "Accent")+
  theme_void()
dev.off()

# Compute and plot TSSenrichment for combined dataset
integrated <- TSSEnrichment(integrated, fast = FALSE, assay = "MACS2")
TSSPlot(integrated, assay = "MACS2", group.by = "donor")


#### Finding DA peaks #####
DefaultAssay(integrated) <- 'MACS2'
Annotation(integrated) <- Annotation(integrated[["ATAC"]])

levels(integrated) <- c("CD56bright", "EarlyCD56dim", "CD56dim",  "Adaptive")

# Peaks unique to cluster
unique_peaks <- FindAllMarkers(integrated, only.pos = TRUE, logfc.threshold = 0.10, 
                               latent.vars = 'nCount_MACS2', test.use = "LR")

saveRDS(unique_peaks, file = paste0(path, "unique_peaks"))
unique_peaks <- readRDS(file = paste0(path, "unique_peaks"))

unique_peaks <- unique_peaks %>% dplyr::filter(p_val_adj < 0.05)



# Pairwise comparisons
CD56dim_vs_CD56bright <- FindMarkers(integrated, ident.1 = c("CD56dim"),
                                     ident.2 = c("CD56bright"), min.pct = 0.05, test.use = 'LR',logfc.threshold = 0.10, 
                                     latent.vars = 'nCount_MACS2')

adaptive_vs_dim <- FindMarkers(integrated, ident.1 = c("Adaptive"),
                               ident.2 = c("CD56dim"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.10,
                               latent.vars = 'nCount_MACS2', features = VariableFeatures(integrated))

CD56dim_peaks <- FindMarkers(integrated, ident.1 = c("CD56dim"),
                             ident.2 = c("CD56bright"), min.pct = 0.05, test.use = 'LR',logfc.threshold = 0.10, 
                             latent.vars = 'nCount_MACS2', only.pos = T)

saveRDS(CD56dim_vs_CD56bright, file = paste0(path, "CD56dim_vs_CD56bright"))
saveRDS(adaptive_vs_dim, file = paste0(path, "adaptive_vs_dim"))
saveRDS(CD56dim_peaks, file = paste0(path, "CD56dim_peaks"))


CD56dim_vs_CD56bright <- readRDS(paste0(path, "CD56dim_vs_CD56bright"))
adaptive_vs_dim <- readRDS(paste0(path, "adaptive_vs_dim"))
CD56dim_peaks <- readRDS(paste0(path, "CD56dim_peaks"))

# Filter peaks for p_adj < 0.05
CD56dim_vs_CD56bright<- CD56dim_vs_CD56bright %>% dplyr::filter(p_val_adj < 0.05)
adaptive_vs_dim <- adaptive_vs_dim %>% dplyr::filter(p_val_adj < 0.05)
CD56dim_peaks <- CD56dim_peaks%>% dplyr::filter(p_val_adj < 0.05)

# Add annotation
adaptive_vs_dim$region <- rownames(adaptive_vs_dim)
temp <- ClosestFeature(integrated, regions = adaptive_vs_dim$region)
adaptive_vs_dim$closest_gene <- temp$gene_name

CD56dim_vs_CD56bright$region <- rownames(CD56dim_vs_CD56bright)
temp <- ClosestFeature(integrated, regions = CD56dim_vs_CD56bright$region)
CD56dim_vs_CD56bright$closest_gene <- temp$gene_name

# Pseudo-bulk heatmap
CD56dim_peaks$region <- rownames(CD56dim_peaks)
heatmap_peaks <- unique(c(unique_peaks[unique_peaks$cluster%in%c("CD56bright", "EarlyCD56dim"),"gene"],
                          CD56dim_peaks$region[!CD56dim_peaks$region%in%unique_peaks$gene],
                          unique_peaks[unique_peaks$cluster%in%c("CD56dim", "Adaptive"),"gene"]))

averaged_peaks <- AverageExpression(integrated, features = heatmap_peaks, slot = "data", assays = "MACS2")

averaged_peaks <- as.data.frame(averaged_peaks$MACS2)

averaged_peaks$region <- rownames(averaged_peaks)

averaged_peaks_tidy <- gather(averaged_peaks, key = "group", value = "accessbility", -region)

pdf(file = "Plots/Fig1_ATAC_heatmap_peaks.pdf", paper = "a4")
averaged_peaks_tidy  %>% tidy_heatmap(
  rows = region,
  columns = group,
  values = accessbility,
  cluster_rows = F,
  cluster_cols = F,
  scale = "row",
  colors =  rev(mapal))
dev.off()


#### Enhancers defined by correlation with gene expression
DefaultAssay(integrated) <- 'MACS2'

# Plot Links
links_df <- as.data.frame(Links(integrated))

# Pick some interesting links to label

genestolabel <- c("KLRC2", "KLRC1", "RUNX2", "ZBTB38", "ZBTB16", "TCF7", "ZEB2", "PRF1", "MYC", "BCL11B", "FCER1G", "IL7R")


links_df %>% ggplot(aes(x = score, y = -log10(pvalue)))+
  geom_point(alpha = 0.5)+
  geom_text_repel(aes(x = score, y = -log10(pvalue), label = gene),
                  data = rbind(links_df[links_df$gene=="KLRC2"&links_df$pvalue<1e-6,],
                               links_df[links_df$gene=="KLRC1"&links_df$pvalue<1e-5,],
                               links_df[links_df$gene=="RUNX2"&links_df$pvalue<1e-6,],
                               links_df[links_df$gene=="ZBTB38"&links_df$pvalue<1e-4,],
                               links_df[links_df$gene=="ZBTB16"&links_df$pvalue<1e-5,],
                               links_df[links_df$gene=="TCF7"&links_df$pvalue<1e-5,],
                               links_df[links_df$gene=="ZEB2"&links_df$pvalue<1e-5,],
                               links_df[links_df$gene=="PRF1"&links_df$pvalue<1e-2,],
                               links_df[links_df$gene=="MYC"&links_df$pvalue<1e-5,],
                               links_df[links_df$gene=="BCL11B"&links_df$pvalue<1e-5,],
                               links_df[links_df$gene=="FCER1G"&links_df$pvalue<1e-6,],
                               links_df[links_df$gene=="IL7R"&links_df$pvalue<1e-5,]),
                  size = 5, box.padding = 0.7, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0,50))+
  xlab("Pearson correlation coefficient")+
  ylab("-Log10(pvalue)")+
  theme_classic()


# Plot representative links
# TCF7
ranges.show <- c(StringToGRanges("chr5-134117873-134118572"),
                StringToGRanges("chr5-134129430-134129958"),
                StringToGRanges("chr5-134131112-134132878"),
                StringToGRanges("chr5-134122344-134123036"))

ranges.show$color <- c(clusterpal[1],
                       clusterpal[1],
                       clusterpal[2],
                       clusterpal[2])
CoveragePlot(
  object = integrated,
  region = "TCF7",
  features = "TCF7",
  expression.assay = "RNA_imputed",
  extend.upstream = 1000,
  extend.downstream = 1000,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

ExpressionPlot(
  object = integrated,
  features = c("TCF7"),
  assay = "RNA_imputed"
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

LinkPlot(integrated,region = "chr5-134113680-134149228")+scale_color_gradient2(
  low = "cyan",
  mid = "grey",
  high = "purple",
  midpoint = 0,
  limits = c(-0.11, 0.33)
)


CoveragePlot(
  object = integrated,
  region = "chr12-10410000-10460000",
  features = "KLRC2",
  expression.assay = "RNA_imputed",
  extend.upstream = 1000,
  extend.downstream = 1000
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

CoveragePlot(
  object = integrated,
  region = "chr12-10410000-10460000",
  features = "KLRC1",
  expression.assay = "RNA_imputed",
  extend.upstream = 1000,
  extend.downstream = 1000
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

LinkPlot(integrated,region = "chr12-10410000-10460000")+scale_color_gradient2(
  low = "cyan",
  mid = "grey",
  high = "purple",
  midpoint = 0,
  limits = c(-0.11, 0.33)
)

ExpressionPlot(
  object = integrated,
  features = c("KLRC1"),
  assay = "RNA_imputed"
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))
ExpressionPlot(
  object = integrated,
  features = c("KLRC2"),
  assay = "RNA_imputed"
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))



CoveragePlot(
  object = integrated,
  region = "chr2-144450000-144550000",
  features = "ZEB2",
  expression.assay = "RNA_imputed",
  extend.upstream = 10000,
  extend.downstream = 80000
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))


LinkPlot(integrated,region = "chr2-144440000-144630000")+scale_color_gradient2(
  low = "cyan",
  mid = "grey",
  high = "purple",
  midpoint = 0,
  limits = c(-0.11, 0.33)
)

ExpressionPlot(
  object = integrated,
  features = c("ZEB2"),
  assay = "RNA_imputed",
  )&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

CoveragePlot(
  object = integrated,
  region = "chr6-45300000-45600000",
  features = "RUNX2",
  expression.assay = "RNA_imputed",
  extend.upstream = 0,
  extend.downstream = 0
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))
ExpressionPlot(
  object = integrated,
  features = c("RUNX2"),
  assay = "RNA_imputed",
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))
LinkPlot(integrated,region = "chr6-45300000-45600000")+scale_color_gradient2(
  low = "cyan",
  mid = "grey",
  high = "purple",
  midpoint = 0,
  limits = c(-0.11, 0.33)
)


# Look at the motifs in the plotted peaks
# Get motifs in these peak
get_motifs_in_peaks <- function(seu, peak){
  lapply(peak, function(x){
    motifs <- names(which(seu[["MACS2"]]@motifs@data[x,]==1))
    motif_names <- ConvertMotifID(seu, id = motifs)
    return(cbind(motifs, motif_names))
  })
}

# KLRC2: 
# chr12-10411032-10411786, chr12-10420034-10420741, chr12-10425313-10425872,
#  chr12-10429318-10429642, chr12-10430415-10431014,   chr12-10435351-10436206
links_df[links_df$gene=="KLRC2",]

KLRC2_motifs <- get_motifs_in_peaks(integrated,
                            c("chr12-10411032-10411786", "chr12-10420034-10420741", "chr12-10425313-10425872",
                              "chr12-10429318-10429642", "chr12-10430415-10431014", "chr12-10435351-10436206"))

# "chr12-10425313-10425872" (which also correlates inversely with KLRC1) and
# "chr12-10430415-10431014" contain AP-1 motifs

# KLRC1: chr12-10422214-10422452, chr12-10425313-10425872, chr12-10437653-10437874, chr12-10453028-10453564
links_df[links_df$gene=="KLRC1",]


KLRC1_motifs <- get_motifs_in_peaks(integrated,
                            c("chr12-10422214-10422452",
                             "chr12-10425313-10425872",
                              "chr12-10437653-10437874",
                              "chr12-10453028-10453564"))

# TCF7, Only CD56bright: chr5-134117873-134118572,, chr5-134129430-134129958
# CD56bright and early dim: chr5-134122344-134123036 (bothchr5-134132043-134132878 (both)
TCF7_motifs <- get_motifs_in_peaks(integrated,
                            c("chr5-134117873-134118572",
                              "chr5-134122344-134123036",
                              "chr5-134129430-134129958",
                              "chr5-134132043-134132878"))




###### Fig. 2: CMV+ vs CMV-  ######

# NKG2C plot on integrated
ASAP_cells <- WhichCells(integrated, expression =  experiment%in%c("mtASAP1", "mtASAP2"))

FeaturePlot(integrated, features = "Anti-PE", cells = ASAP_cells, cols = colorscale, min.cutoff = "q1",  max.cutoff = "q99", order =TRUE, pt.size = 0.4, coord.fixed = TRUE)&NoLegend()&
  ggtitle(label = NULL)

#### NKG2C+ vs NKG2C- analysis was moved to script 06 to 
#### include additional CMV- donors from mtASAP5

# Analysis of NKG2C+ from CMV+ alone
RidgePlot(subset(integrated, cells = WhichCells(integrated, expression = experiment%in%c("mtASAP1", "mtASAP2"))), features = "Anti-PE")
NKG2Cpos <- integrated[,(integrated@assays$ADT@data["Anti-PE",]>0.5)]

NKG2Cpos_CMVpos <- subset(NKG2Cpos, subset = serostatus == "CMVpos")
NKG2Cpos_CMVpos <- RunUMAP(NKG2Cpos_CMVpos, dims = 2:10, reduction = "integrated_lsi")

DefaultAssay(NKG2Cpos_CMVpos) <- "RNA_atac"
DimPlot(NKG2Cpos_CMVpos, cols = c(clusterpal[1:3], adaptivepal))&NoLegend()
cells_downsampled3 <- sample(Cells(NKG2Cpos_CMVpos), size = 7000)

pdf(file = "Plots/SupplFig2_ATAC_representatives.pdf", paper = "a4r")
FeaturePlot(NKG2Cpos_CMVpos, features = c("JAKMIP1", "ZBTB38", "ZNF516", "ZBTB16"), cols = blueyellow,
            order =TRUE, slot = "data", min.cutoff = "q1", max.cutoff = "q95", coord.fixed = T, cells = cells_downsampled3)&NoLegend()&NoAxes()
FeaturePlot(NKG2Cpos_CMVpos, features = c("JAKMIP1", "ZBTB38", "ZNF516", "ZBTB16"), cols = blueyellow,
            order =F, slot = "data", min.cutoff = "q1", max.cutoff = "q95", coord.fixed = T, cells = cells_downsampled3)&NoLegend()&NoAxes()
dev.off()


# Plot signature from full datasets

# Activities
averaged_gene_act_2C_CMVpos <- AverageExpression(NKG2Cpos_CMVpos, features = heatmap_acts, slot = "data", assays = "RNA_atac")

averaged_gene_act_2C_CMVpos <- as.data.frame(averaged_gene_act_2C_CMVpos$RNA_atac)

averaged_gene_act_2C_CMVpos$gene <- rownames(averaged_gene_act_2C_CMVpos)

averaged_gene_act_2C_CMVpos_tidy <- gather(averaged_gene_act_2C_CMVpos, key = "group", value = "activity", -gene)

pdf(file = "Plots/SupplFig2_ATAC_heatmap_activity.pdf", paper = "a4")
averaged_gene_act_2C_CMVpos_tidy  %>% tidy_heatmap(
  rows = gene,
  columns = group,
  values = activity,
  cluster_rows = F,
  cluster_cols = F,
  scale = "row",
  colors =  blueyellow)
dev.off()

# Peaks  
averaged_peaks_2C_CMVpos <- AverageExpression(NKG2Cpos_CMVpos, features = heatmap_peaks, slot = "data", assays = "MACS2")

averaged_peaks_2C_CMVpos <- as.data.frame(averaged_peaks_2C_CMVpos$MACS2)

averaged_peaks_2C_CMVpos$region <- rownames(averaged_peaks_2C_CMVpos)

averaged_peaks_2C_CMVpos_tidy <- gather(averaged_peaks_2C_CMVpos, key = "group", value = "accessbility", -region)

pdf(file = "Plots/SupplFig2_ATAC_heatmap_peaks.pdf", paper = "a4")
averaged_peaks_2C_CMVpos_tidy  %>% tidy_heatmap(
  rows = region,
  columns = group,
  values = accessbility,
  cluster_rows = F,
  cluster_cols = F,
  scale = "row",
  colors =  rev(mapal))
dev.off()




###### Fig. 3 Motif Analysis/AP-1 ######

# Find motifs with differential chromvar activity
DefaultAssay(integrated) <- 'chromvar_names'


differential.activity <- FindAllMarkers(
  object = integrated,
  only.pos = T,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

dim_vs_bright_motifact <- FindMarkers(
  object = integrated,
  ident.1 = "CD56dim",
  ident.2 = "CD56bright",
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0
)
dim_vs_bright_motifact$motif <- rownames(dim_vs_bright_motifact)

adaptive_vs_dim_motifact <- FindMarkers(
  object = integrated,
  ident.1 = "Adaptive",
  ident.2 = "CD56dim",
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold = 0
)
adaptive_vs_dim_motifact$motif <- rownames(adaptive_vs_dim_motifact)

top_activity <- differential.activity %>% group_by(cluster) %>% top_n(n = 25, wt = avg_diff)

# Pseudobulk heatmap
averaged_motif <- AverageExpression(integrated, features = top_activity$gene, slot = "data", assays = "chromvar_names", return.seurat = TRUE)

averaged_motif <- as.data.frame(GetAssayData(averaged_motif, slot = "data"))

averaged_motif$motif <- rownames(averaged_motif)

averaged_motif_tidy <- gather(averaged_motif, key = "group", value = "activity", -motif)


pdf(file = "Plots/Fig3_ATAC_heatmap_motif_activity.pdf", paper = "a4")
averaged_motif_tidy  %>% tidy_heatmap(
  rows = motif,
  columns = group,
  values = activity,
  cluster_rows = T,
  cluster_cols = F,
  scale = "none",
  colors =  rev(mapal),
  color_legend_min = -5,
  color_legend_max = 5
)
dev.off()


# Volcanoplot of pair-wise comparisons
# Bright vs Dim
motifs_to_label <- c("RUNX2", "TCF7", "REL", "CTCF", "LEF1",  "ZEB1", "Stat4", "TBX1")

dim_vs_bright_motifact %>% 
  ggplot(aes(x = avg_diff, y = -log10(p_val_adj+1e-320)))+
  geom_point()+
  scale_color_manual(values = c("black", "darkred"))+
  geom_text_repel(aes(x = avg_diff, y = -log10(p_val_adj+1e-320), label = motif),
                  data = dim_vs_bright_motifact[dim_vs_bright_motifact$motif%in%c(motifs_to_label),],
                  size = 5, box.padding = 1, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 330))+
  scale_x_continuous(limits = c(-4.5, 4.5))+
  xlab(label = "CD56dim vs CD56bright (Avg z-score diff)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

# Adaptive vs dim

motifs_to_label <- c("FOS::JUNB")

adaptive_vs_dim_motifact %>% 
  ggplot(aes(x = avg_diff, y = -log10(p_val_adj+1e-320)))+
  geom_point()+
  scale_color_manual(values = c("black", "darkred"))+
  geom_text_repel(aes(x = avg_diff, y = -log10(p_val_adj+1e-320), label = motif),
                  data = adaptive_vs_dim_motifact[adaptive_vs_dim_motifact$motif%in%c(motifs_to_label),],
                  size = 5, box.padding = 1, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 330))+
  scale_x_continuous(limits = c(-4.5, 4.5))+
  xlab(label = "Adaptive vs CD56dim (Avg z-score diff)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()


# AP1 representatives
cells_plot <- sample(Cells(integrated), size = 10000)
pdf(file = "Plots/Fig3_ATAC_AP1_featureplot.pdf", paper = "a4")
FeaturePlot(integrated, features = "FOS::JUNB", cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q95", coord.fixed = TRUE, cells = cells_plot)+NoLegend()
dev.off()

# Plot some bright/dim specific examples
pdf(file = "Plots/SupplFig3_ATAC_featureplots.pdf", paper = "a4")
FeaturePlot(integrated, features = "RUNX2", cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE, coord.fixed = TRUE, cells = cells_plot)+NoLegend()
FeaturePlot(integrated, features = "REL", cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE, coord.fixed = TRUE, cells = cells_plot)+NoLegend()
FeaturePlot(integrated, features = "CTCF", cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE, coord.fixed = TRUE, cells = cells_plot)+NoLegend()
FeaturePlot(integrated, features = "ZEB1", cols = colorscale,
            min.cutoff = "q1", max.cutoff = "q99", order = TRUE, coord.fixed = TRUE, cells = cells_plot)+NoLegend()
dev.off()



##### Which peaks contain motifs? #####

# Annotate peaks
Peaks_to_motif <- function(seu, motif, peaks){
  # Get peaks
  motif_peaks <- seu[["MACS2"]]@motifs@data[rownames(peaks), motif]==1
  motif_peaks <- names(which(motif_peaks))
  # Add FC etc
  motif_peaks <- peaks[motif_peaks,]
  # Assign genes
  motif_genes <- ClosestFeature(seu, regions =rownames(motif_peaks))
  motif_peaks$gene <- motif_genes$gene_name
  motif_peaks$region <- rownames(motif_peaks)
  return(motif_peaks)
}

# Adaptive vs dim: AP-1 enrichment
AP1_peaks <- Peaks_to_motif(integrated, motif = "MA1130.1", peaks = adaptive_vs_dim)
adaptive_vs_dim$AP1 <- FALSE
adaptive_vs_dim$AP1[rownames(adaptive_vs_dim)%in%rownames(AP1_peaks)] <- TRUE

layer1_adaptive <- adaptive_vs_dim[adaptive_vs_dim$AP1,]
layer2_adaptive <- adaptive_vs_dim[!adaptive_vs_dim$AP1,]

pdf(file = "Plots/Fig3_ATAC_peaks_AP1.pdf")
peakstolabel <- c( "AIM2", "CADM1", "IFNG", "RPTOR")
ggplot()+
  geom_point(data= layer2_adaptive, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-320)), color = "grey")+
  geom_point(data = layer1_adaptive, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-320)), color = "darkred")+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-320), label = closest_gene), color = "darkred",
                  data = adaptive_vs_dim[adaptive_vs_dim$closest_gene%in%c(peakstolabel)&adaptive_vs_dim$AP1,],
                  size = 5, box.padding = 0.9, max.overlaps = Inf, point.padding = 0.5)+
  #scale_y_continuous(limits = c(0, 40))+
  scale_x_continuous(limits = c(-0.6, 0.6))+
  xlab(label = "Adaptive vs CD56dim (log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()
dev.off()


# Enrichment of motifs in adaptive peaks
DefaultAssay(integrated) <- "MACS2"
adaptive_peaks <- adaptive_vs_dim %>% dplyr::filter(avg_log2FC >0)


enriched.motifs_adaptive <- FindMotifs(
  object = integrated,
  features = rownames(adaptive_peaks)
)


enriched.motifs_adaptive <- enriched.motifs_adaptive %>% dplyr::filter(pvalue < 0.05)
enriched.motifs_adaptive <- enriched.motifs_adaptive %>% arrange(pvalue) %>% mutate(rank = row_number())

# Plot ranked enrichment
AP1_to_label <- c("BATF3", "BATF", "BATF::JUN", "JUN(var.2)", "FOSL1::JUND",
                  "FOSL1::JUN", "FOS::JUN", "JUNB", "FOS::JUND", "FOSL1", "FOSL1::JUNB",
                  "JUND", "FOSL2::JUND", "FOSB::JUNB", "FOSL2::JUN", "FOSL2", "JUN::JUNB", "FOSL2::JUNB",
                  "FOS::JUNB", "NFE2", "JDP2", "FOS")

TFs_to_label <- c("FOS::JUN", "TBX1", "EOMES", "TBX21")

pdf(file = "Plots/Fig3_ATAC_motif_enrichment.pdf")
enriched.motifs_adaptive %>%  ggplot(aes(x = rank, y = -log10(pvalue)))+
  geom_point(size = 1)+
  geom_point(color = "darkred", size = 1, data =  enriched.motifs_adaptive[enriched.motifs_adaptive$motif.name%in%AP1_to_label,])+
  geom_text_repel(aes(x = rank, y = -log10(pvalue), label = motif.name),
                  data = enriched.motifs_adaptive[enriched.motifs_adaptive$motif.name%in%TFs_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 0.7)+
  xlab(label = "Rank Motif")+
  ylab(label = "-log10(pvalue)")+
  theme_classic()
dev.off()

# Enrichment of motifs in Bright peaks
DefaultAssay(integrated) <- "MACS2"
bright_peaks <- CD56dim_vs_CD56bright %>% dplyr::filter(avg_log2FC <0)


enriched.motifs_bright <- FindMotifs(
  object = integrated,
  features = rownames(bright_peaks)
)


enriched.motifs_bright <- enriched.motifs_bright %>% dplyr::filter(pvalue < 0.05)
enriched.motifs_bright <- enriched.motifs_bright %>% arrange(pvalue) %>% mutate(rank = row_number())

# Plot ranked enrichment
TFs_to_label <- c("ZEB1", "RUNX2", "TCF7", "RELA", "Stat4", "TBX1", "KLF4", "EOMES", "TBX21")

pdf(file = "Plots/SupplFig3_ATAC_motif_enrichment_bright.pdf")
enriched.motifs_bright %>%  ggplot(aes(x = rank, y = -log10(pvalue)))+
  geom_point(size = 1)+
  geom_text_repel(aes(x = rank, y = -log10(pvalue), label = motif.name),
                  data = enriched.motifs_bright[enriched.motifs_bright$motif.name%in%TFs_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 0.7)+
  xlab(label = "Rank Motif")+
  ylab(label = "-log10(pvalue)")+
  theme_classic()
dev.off()

# Enrichment of motifs in Dim peaks
DefaultAssay(integrated) <- "MACS2"
dim_peaks <- CD56dim_vs_CD56bright %>% dplyr::filter(avg_log2FC >0)


enriched.motifs_dim <- FindMotifs(
  object = integrated,
  features = rownames(dim_peaks)
)


enriched.motifs_dim <- enriched.motifs_dim %>% dplyr::filter(pvalue < 0.05)
enriched.motifs_dim <- enriched.motifs_dim %>% arrange(pvalue) %>% mutate(rank = row_number())

# Plot ranked enrichment
TFs_to_label <- c("ZNF460", "KLF4", "EOMES", "TBX1", "TBX21")

pdf(file = "Plots/SupplFig3_ATAC_motif_enrichment_dim.pdf")
enriched.motifs_dim %>%  ggplot(aes(x = rank, y = -log10(pvalue)))+
  geom_point(size = 1)+
  geom_text_repel(aes(x = rank, y = -log10(pvalue), label = motif.name),
                  data = enriched.motifs_dim[enriched.motifs_dim$motif.name%in%TFs_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 0.7)+
  xlab(label = "Rank Motif")+
  ylab(label = "-log10(pvalue)")+
  theme_classic()
dev.off()



# Plot some representative peaks
DefaultAssay(integrated) <- "MACS2"

# IFNG enhancers
# Adaptive/Dim Proximal: chr12-68482371-68483250
# Adaptive: Distal "chr12-68119138-68483879"
# IFNG peak of bright chr12-68175228-68176254
ranges.show <- StringToGRanges("chr12-68482371-68483879")
ranges.show$color <- "darkred"
p1 <- CoveragePlot(
  object = integrated,
  region = "chr12-68482371-68483250",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  region.highlight = ranges.show,
  ymax = 390
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

ranges.show <- c(StringToGRanges("chr12-68119138-68120343"),
                 StringToGRanges("chr12-68175228-68176254"))
ranges.show$color <- c("darkred", "green")
p2 <- CoveragePlot(
  object = integrated,
  region = "chr12-68119138-68176254",
  extend.upstream = 1000,
  extend.downstream = 1000,
  links = F,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

# Distal enhancer also contains AP-1 motif, which is however split
# into smaller, directly neighboring peak:
grep("chr12-68*",names(which(integrated[["MACS2"]]@motifs@data[, "MA1134.1"]==1)), value = TRUE)
# Coordinates: "chr12-68483580-68483879"

# Stat4 in IFNG locus: Contained in bright-specific peak!
grep("chr12-68*",names(which(integrated[["MACS2"]]@motifs@data[, "MA0518.1"]==1)), value = TRUE)



p1&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))
p2&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

# Aim2 peak
ranges.show <- StringToGRanges("chr1-159076355-159077686")
ranges.show$color <- "darkred"
CoveragePlot(
  object = integrated,
  region = "chr1-159076355-159077686",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

# Cadm1 peak
ranges.show <- StringToGRanges("chr11-115222052-115223633")
ranges.show$color <- "darkred"
CoveragePlot(
  object = integrated,
  region = "chr11-115222052-115223633",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))


# RPTOR chr17-80891280-80892607
ranges.show <- StringToGRanges("chr17-80891280-80892607")
ranges.show$color <- "darkred"
CoveragePlot(
  object = integrated,
  region = "chr17-80891280-80892607",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))