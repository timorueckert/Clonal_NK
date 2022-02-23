# Downstream analysis: Response patterns, surface features, peaks, motifs

library(Signac)
library(ggplot)
library(ggrepel)

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- brewer.pal(8, "Dark2")
pairedpal <- brewer.pal(12, "Paired")
NKG2Cpal <- c("#40BAD5","#120136")
clusterpal2 <- brewer.pal(8, "Set1")
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

# Import processed dataset
NK_act <- readRDS(paste0(path, "NK_act_dim"))

# Plot TSS plot per condition
TSSPlot(NK_act, assay = "peaks")

# Plot UMAP by condition, cluster and donor
pdf(file = "Plots/Fig4_Dimplots1.pdf", paper = "a4r")
DimPlot(NK_act, group.by = "condition", cols = pairedpal[c(3,4,7,8)], pt.size = 0.4)&NoLegend()
DimPlot(NK_act,  cols = pairedpal[c(1,5,2,6)], pt.size = 0.4, shuffle = TRUE)&NoLegend()
DimPlot(NK_act, group.by = "donor", pt.size = 0.4, cols = clusterpal2[c(2,5)])&NoLegend()
dev.off()

# Inspect 4-1BB and NKG2A expression across clusters
pdf(file = "Plots/Fig4_violin.pdf", paper = "a4r")
VlnPlot(NK_act, features = c("CD137"), pt.size = 0, cols = pairedpal[c(1,5,2,6)], log = T)
VlnPlot(NK_act, features = c("anti-Biotin"), pt.size = 0, cols = pairedpal[c(1,5,2,6)], log = T)
dev.off()


# Test for significance
NKG2A_df <- FetchData(NK_act, vars = c("anti-Biotin", "annotation"))
pairwise.wilcox.test(NKG2A_df$`ADT_anti-Biotin`,NKG2A_df$annotation, p.adjust="BH")

CD137_df <- FetchData(NK_act, vars = c("CD137", "annotation"))
pairwise.wilcox.test(CD137_df$ADT_CD137, CD137_df$annotation, p.adjust="BH")


##### DAP analysis & enriched motifs ####
# Find DAPs
DefaultAssay(NK_act) <- "MACS2"

peaks_both <- FindMarkers(NK_act, ident.1 = c("Peptide+Cytokines"),
                          ident.2 = c("Control"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.00,
                          latent.vars = 'nCount_MACS2')

peaks_peptide <- FindMarkers(NK_act, ident.1 = c("Peptide"),
                             ident.2 = c("Control"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.00,
                             latent.vars = 'nCount_MACS2')

peaks_cytokines <- FindMarkers(NK_act, ident.1 = c("Cytokines"),
                               ident.2 = c("Control"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.00,
                               latent.vars = 'nCount_MACS2')

peaks_synergism_peptide <- FindMarkers(NK_act, ident.1 = c("Peptide+Cytokines"),
                                       ident.2 = c("Peptide"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.00,
                                       latent.vars = 'nCount_MACS2')

peaks_synergism_cytokines <- FindMarkers(NK_act, ident.1 = c("Peptide+Cytokines"),
                                         ident.2 = c("Cytokines"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.00,
                                         latent.vars = 'nCount_MACS2')


# Add region as column
peaks_both$region <- rownames(peaks_both)
peaks_peptide$region <- rownames(peaks_peptide)
peaks_cytokines$region <- rownames(peaks_cytokines)

# Export
saveRDS(peaks_both, paste0(path, "peaks_both"))
saveRDS(peaks_peptide, paste0(path, "peaks_peptide"))
saveRDS(peaks_cytokines, paste0(path, "peaks_cytokines"))
saveRDS(peaks_synergism_peptide , paste0(path, "peaks_synergism_peptide "))
saveRDS(peaks_synergism_cytokines, paste0(path, "peaks_synergism_cytokines"))

# Import
peaks_both <- readRDS(paste0(path, "peaks_both"))
peaks_peptide <- readRDS(paste0(path, "peaks_peptide"))
peaks_cytokines <- readRDS(paste0(path, "peaks_cytokines"))
peaks_synergism_peptide <- readRDS(paste0(path, "peaks_synergism_peptide "))
peaks_synergism_cytokines <- readRDS(paste0(path, "peaks_synergism_cytokines"))



# Volcanoplots of pairwise comparisons

# Function to annotate peak tables with genes and motifs (Peaks_to_motif function is in ex vivo figures script)
add_peak_annotation <- function(peaks, motifs, seu){
  peaks$region <- rownames(peaks)
  annotation <- ClosestFeature(seu, regions = peaks$region)
  peaks$closest_gene <- annotation$gene_name
  motif_df <- data.frame(lapply(motifs, function(x){
    motif_peaks <- Peaks_to_motif(seu, motifs, peaks)
    motif_logical <- peaks$region%in%rownames(motif_peaks)
  }))
  colnames(motif_df) <- c(motifs)
  cbind(peaks, motif_df)
}

# Add Diffexpr column
peaks_both$diffexp <- "NO"
peaks_both$diffexp[peaks_both$p_val_adj < 0.05 & peaks_both$avg_log2FC > 0.1] <- "UP"
peaks_both$diffexp[peaks_both$p_val_adj < 0.05 & peaks_both$avg_log2FC < -0.1] <- "DOWN"

peaks_peptide$diffexp <- "NO"
peaks_peptide$diffexp[peaks_peptide$p_val_adj < 0.05 & peaks_peptide$avg_log2FC > 0.1] <- "UP"
peaks_peptide$diffexp[peaks_peptide$p_val_adj < 0.05 & peaks_peptide$avg_log2FC < -0.1] <- "DOWN"

peaks_cytokines$diffexp <- "NO"
peaks_cytokines$diffexp[peaks_cytokines$p_val_adj < 0.05 & peaks_cytokines$avg_log2FC > 0.1] <- "UP"
peaks_cytokines$diffexp[peaks_cytokines$p_val_adj < 0.05 & peaks_cytokines$avg_log2FC < -0.1] <- "DOWN"

peaks_synergism_peptide$diffexp <- "NO"
peaks_synergism_peptide$diffexp[peaks_synergism_peptide$p_val_adj < 0.05 & peaks_synergism_peptide$avg_log2FC > 0.1] <- "UP"
peaks_synergism_peptide$diffexp[peaks_synergism_peptide$p_val_adj < 0.05 & peaks_synergism_peptide$avg_log2FC < -0.1] <- "DOWN"

peaks_synergism_cytokines$diffexp <- "NO"
peaks_synergism_cytokines$diffexp[peaks_synergism_cytokines$p_val_adj < 0.05 & peaks_synergism_cytokines$avg_log2FC > 0.1] <- "UP"
peaks_synergism_cytokines$diffexp[peaks_synergism_cytokines$p_val_adj < 0.05 & peaks_synergism_cytokines$avg_log2FC < -0.1] <- "DOWN"


# Layer so that significant genes are on top
layer2_both <- peaks_both[peaks_both$diffexp%in%c("UP", "DOWN"),]
layer1_both <- peaks_both[peaks_both$diffexp%in%c("NO"),]

layer2_peptide <- peaks_peptide[peaks_peptide$diffexp%in%c("UP", "DOWN"),]
layer1_peptide <- peaks_peptide[peaks_peptide$diffexp%in%c("NO"),]

layer2_cytokines <- peaks_cytokines[peaks_cytokines$diffexp%in%c("UP", "DOWN"),]
layer1_cytokines <- peaks_cytokines[peaks_cytokines$diffexp%in%c("NO"),]

layer2_synergism_peptide <- peaks_synergism_peptide[peaks_synergism_peptide$diffexp%in%c("UP", "DOWN"),]
layer1_synergism_peptide <- peaks_synergism_peptide[peaks_synergism_peptide$diffexp%in%c("NO"),]

layer2_synergism_cytokines <- peaks_synergism_cytokines[peaks_synergism_cytokines$diffexp%in%c("UP", "DOWN"),]
layer1_synergism_cytokines <- peaks_synergism_cytokines[peaks_synergism_cytokines$diffexp%in%c("NO"),]

# LFL+IL-12/18
peaks_both <- add_peak_annotation(peaks_both,"MA1142.1", NK_act)
peakstolabel <- c("ZEB2", "IFNG", "ARID5B", "ZBTB16", "AIM2", "TNFRSF9", "CD2")
ggplot()+
  geom_point(data = layer1_both, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  geom_point(data = layer2_both, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = closest_gene),
                  data = peaks_both[peaks_both$closest_gene%in%c(peakstolabel)&peaks_both$p_val_adj < 0.05 & abs(peaks_both$avg_log2FC) > 0.1,],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_x_continuous(limits = c(-0.75, 0.75))+
  xlab(label = "LFL+IL-12+IL18 vs Control (log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

sum(!peaks_both$diffexp=="NO")


# Peptide
peaks_peptide <- add_peak_annotation(peaks_peptide,"MA1142.1", NK_act)
peakstolabel <- c("ZEB2", "IFNG", "ARID5B", "ZBTB16", "AIM2", "TNFRSF9", "CD2")

ggplot()+
  geom_point(data = layer1_peptide, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  geom_point(data = layer2_peptide, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = closest_gene),
                  data = peaks_peptide[peaks_peptide$closest_gene%in%c(peakstolabel)&peaks_peptide$p_val_adj < 0.05 & abs(peaks_peptide$avg_log2FC) > 0.1,],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_x_continuous(limits = c(-0.75, 0.75))+
  xlab(label = "LFL vs Control (log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

sum(!peaks_peptide$diffexp=="NO")

# Cytokines 
peaks_cytokines <- add_peak_annotation(peaks_cytokines,"MA1142.1", NK_act)

ggplot()+
  geom_point(data = layer1_cytokines, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  geom_point(data = layer2_cytokines, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = closest_gene),
                  data = peaks_cytokines[peaks_cytokines$closest_gene%in%c(peakstolabel)&peaks_cytokines$p_val_adj < 0.05 & abs(peaks_cytokines$avg_log2FC) > 0.1,],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_x_continuous(limits = c(-0.75, 0.75))+
  xlab(label = "IL-12+IL18 vs Control (log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

sum(!peaks_cytokines$diffexp=="NO")


# Synergism vs peptide
peaks_synergism_peptide <- add_peak_annotation(peaks_synergism_peptide,"MA1142.1", NK_act)
peakstolabel <- c("ZEB2", "IFNG", "ARID5B", "ZBTB16", "AIM2", "TNFRSF9", "CD2")
ggplot()+
  geom_point(data = layer1_synergism_peptide, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  geom_point(data = layer2_synergism_peptide, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = closest_gene),
                  data = peaks_synergism_peptide[peaks_synergism_peptide$closest_gene%in%c(peakstolabel)&peaks_synergism_peptide$p_val_adj < 0.05 & abs(peaks_synergism_peptide$avg_log2FC) > 0.1,],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_x_continuous(limits = c(-0.75, 0.75))+
  xlab(label = "LFL+IL-12+IL18 vs LFL (log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

sum(!peaks_synergism_peptide$diffexp=="NO")

# Synergism vs cytokines
peaks_synergism_cytokines <- add_peak_annotation(peaks_synergism_cytokines,"MA1142.1", NK_act)
peakstolabel <- c("ZEB2", "IFNG", "ARID5B", "ZBTB16", "AIM2", "TNFRSF9", "CD2")
ggplot()+
  geom_point(data = layer1_synergism_cytokines, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  geom_point(data = layer2_synergism_cytokines, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = closest_gene),
                  data = peaks_synergism_cytokines[peaks_synergism_cytokines$closest_gene%in%c(peakstolabel)&peaks_synergism_cytokines$p_val_adj < 0.05 & abs(peaks_synergism_cytokines$avg_log2FC) > 0.1,],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_x_continuous(limits = c(-0.75, 0.75))+
  xlab(label = "LFL+IL-12+IL18 vs IL-12+IL18 (log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

sum(!peaks_synergism_cytokines$diffexp=="NO")

# Filter for statistically significant
peaks_both <- peaks_both %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1)
peaks_peptide <- peaks_peptide %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1)
peaks_cytokines <- peaks_cytokines %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1)
peaks_synergism_peptide  <- peaks_synergism_peptide %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1)
peaks_synergism_cytokines <- peaks_synergism_cytokines %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1)


# Classification of induced peaks  as promoter/exon/intron/intergenig
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annots <- c('hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)
DefaultAssay(NK_act) <- "MACS2"

# Global overview of peak distribution
peaks_meta_both <- ClosestFeature(NK_act, regions = rownames(peaks_both), annotation = annotations)

# Counts types
peaks_meta_both %>% count(type) -> Peak_type_numbers
# Show as pie chart
pdf(file = "Plots/SupplFig4_ATAC_peak_types_activated.pdf")
Peak_type_numbers %>% 
  ggplot(aes(x = "", y = n, fill = type))+
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette = "Accent")+
  theme_void()
dev.off()


##### Plot euler diagrams of peaks in different conditions ####
library(eulerr)
all_peaks <- unique(c(peaks_both$region,
                      peaks_cytokines$region,
                      peaks_peptide$region
))

# Create matrix of logicals representing whether peak is a DAP for each dataset

peaks_both_log <-  all_peaks%in%peaks_both$region
peaks_peptide_log <- all_peaks%in%peaks_peptide$region
peaks_cytokines_log <- all_peaks%in%peaks_cytokines$region

DAP_relationship <- data.frame("Peptide_Cytokines" = peaks_both_log,
                               "Peptide" = peaks_peptide_log,
                               "Cytokines" = peaks_cytokines_log)

rownames(DAP_relationship) <- all_peaks

# Plot as an Euler diagram
fit2 <- euler(DAP_relationship)

plot(fit2,
     quantities = list(type = c("counts", "percent")),
     lty = 1:3,
     labels = list(font = 4),
     fills = c(pairedpal[6], pairedpal[5], pairedpal[2]))





#### Analysis of chromvar activity across conditions ####
# Differential motif activity
DefaultAssay(NK_act) <- "chromvar_names"
motifs_activated <- FindMarkers(NK_act, ident.1 = "Peptide+Cytokines", ident.2 = "Control", mean.fxn=rowMeans, fc.name="avg_diff")
motifs_cytokines <- FindMarkers(NK_act, ident.1 = "Cytokines", ident.2 = "Control", mean.fxn=rowMeans, fc.name="avg_diff")
motifs_synergism <- FindMarkers(NK_act, ident.1 = "Peptide+Cytokines", ident.2 = "Cytokines", mean.fxn=rowMeans, fc.name="avg_diff")
motifs_peptide <- FindMarkers(NK_act, ident.1 = "Peptide", ident.2 = "Control", mean.fxn=rowMeans, fc.name="avg_diff")
motifs_synergism_peptide <- FindMarkers(NK_act, ident.1 = "Peptide+Cytokines", ident.2 = "Peptide", mean.fxn=rowMeans, fc.name="avg_diff")
motifs_activated$motif <- rownames(motifs_activated)
motifs_cytokines$motif <- rownames(motifs_cytokines)
motifs_synergism$motif <- rownames(motifs_synergism)
motifs_peptide$motif <- rownames(motifs_peptide)
motifs_synergism_peptide$motif <- rownames(motifs_synergism_peptide)

motifs_activated$group <- "CytokinesPeptide"
motifs_cytokines$group <- "Cytokines"
motifs_peptide$group <- "Peptide"

# Filter for significant
motifs_activated <- motifs_activated %>% dplyr::filter(p_val_adj < 0.05)
motifs_cytokines <- motifs_cytokines %>% dplyr::filter(p_val_adj < 0.05)
motifs_synergism <- motifs_synergism %>% dplyr::filter(p_val_adj < 0.05)
motifs_peptide <- motifs_peptide %>% dplyr::filter(p_val_adj < 0.05)
motifs_synergism_peptide <- motifs_synergism_peptide %>% dplyr::filter(p_val_adj < 0.05)


# # Plot Heatmap
top_motifs <-  Reduce(function(...) merge(..., all = TRUE),
                      list(motifs_activated,
                           motifs_cytokines,
                           motifs_peptide))
top_motif_names <- top_motifs %>% group_by(group) %>% top_n(n = 50, wt = avg_diff)



# Averaged
library(tidyr)
library(tidyheatmap)

averaged_motif <- AverageExpression(NK_act, features = top_motif_names$motif, slot = "data", assays = "chromvar_names", return.seurat = TRUE)

averaged_motif <- as.data.frame(GetAssayData(averaged_motif, slot = "data"))

averaged_motif$motif <- rownames(averaged_motif)

averaged_motif_tidy <- gather(averaged_motif, key = "group", value = "activity", -motif)


pdf(file = "Plots/Fig4_ATAC_heatmap_motif_activity.pdf", paper = "a4")
averaged_motif_tidy  %>% tidy_heatmap(
  rows = motif,
  columns = group,
  values = activity,
  cluster_rows = T,
  cluster_cols = F,
  scale = "none",
  colors =  rev(mapal),
  color_legend_min = -4,
  color_legend_max = 4
)
dev.off()

# Plot representative
pdf(file = "Plots/SupplFig4_JUNB.pdf", paper = "a4r")
VlnPlot(NK_act, features = "JUNB", c(pairedpal[c(1,5,2,6)]), pt.size = 0)&NoLegend()
dev.off()


# Test for significance
JUNB_df <- FetchData(NK_act, vars = c("JUNB", "annotation"))
pairwise.wilcox.test(JUNB_df$JUNB, JUNB_df$annotation, p.adjust="BH")



# Representative peaks
# ZBTB16
ranges.show <- StringToGRanges("chr11-114073169-114074341")
ranges.show$color <- "darkred"
CoveragePlot(
  object = NK_act,
  region = "chr11-114073169-114074341",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(pairedpal[c(1,5,2,6)]))

# AIM2
ranges.show <- StringToGRanges("chr1-159076355-159077686")
ranges.show$color <- "darkred"
CoveragePlot(
  object = NK_act,
  region = "chr1-159076355-159077686",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(pairedpal[c(1,5,2,6)]))

# ZEB2 chr2-144357786-144359741; chr2-144413683-144414957, chr2-144405415-144406610
ranges.show <- c(StringToGRanges("chr2-144357786-144359741"),
                 StringToGRanges("chr2-144405415-144406610"),
                 StringToGRanges("chr2-144413683-144414957")
                 )
ranges.show$color <- "darkred"
CoveragePlot(
  object = NK_act,
  region = "chr2-144357786-144414957",
  extend.upstream = 1000,
  extend.downstream = 1000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(pairedpal[c(1,5,2,6)]))

# IFNG 
ranges.show <- c(StringToGRanges("chr12-68080332-68081978"),
                 StringToGRanges("chr12-68117462-68117953"),
                 StringToGRanges("chr12-68120996-68121533")
)
ranges.show$color <- "darkred"
CoveragePlot(
  object = NK_act,
  region = "chr12-68080332-68081978",
  extend.upstream = 1000,
  extend.downstream = 80000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(pairedpal[c(1,5,2,6)]))



# ARID5B
ranges.show <- c(StringToGRanges("chr10-61882367-61883992"),
                 StringToGRanges("chr10-61893982-61894619"),
  StringToGRanges("chr10-61942562-61944316"),
  StringToGRanges("chr10-61944623-61945010"),
                 StringToGRanges("chr10-61937440-61938505")
)
ranges.show$color <- "darkred"
CoveragePlot(
  object = NK_act,
  region = "chr10-61882367-61944316",
  extend.upstream = 5000,
  extend.downstream = 5000,
  links = F,
  scale.factor = 1e7,
  region.highlight = ranges.show
)&scale_fill_manual(values = c(pairedpal[c(1,5,2,6)]))
