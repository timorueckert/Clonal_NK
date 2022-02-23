# Downstream analysis: Differential gene expression, comparison of datasets
# separated by CMV serostatus

library(scico)
library(heatmaply)
library(patchwork)

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- met.brewer(name="Isfahan1",n=7,type="discrete")
clusterpal <- clusterpal[c(4,5,7,2)]
NKG2Cpal <- c("#40BAD5","#120136")
adaptivepal <- c("#880E4F",  "#C2185B",  "#E91E63",  "#F06292", "#F8BBD0", "#4A148C", "#7B1FA2", "#9C27B0", "#BA68C8")
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
# From GreenleafLab/ArchR
blueyellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60')
solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')


path <- "~/"

# Import data
integrated_scRNA <- readRDS(paste0(path, "integrated_scRNA"))
CMVpos_scRNA <- readRDS(paste0(path, "CMVpos_scRNA"))
CMVneg_scRNA <- readRDS(paste0(path, "CMVneg_scRNA"))



#### Figure 1
pdf(file = "Plots/Fig1_scRNA_Dimplot.pdf", paper = "a4")
DimPlot(integrated_scRNA, cols = c(clusterpal[1:3], adaptivepal[1], clusterpal[4]))
DimPlot(integrated_scRNA, cols = c(clusterpal[1:3], adaptivepal[1], clusterpal[4]), shuffle = T)&NoLegend()
dev.off()

cells_plot <- sample(Cells(integrated_scRNA), size = 3000)


DefaultAssay(integrated_scRNA) <- "ADT"
pdf(file = "Plots/Fig1_scRNA_CITE.pdf", paper = "a4")
FeaturePlot(
  object = integrated_scRNA,
  features = c("CD56-adt", "CD62L-adt", "CD16-adt", "CD57-adt", "NKp30-adt",  "Anti-Biotin-NKG2A-adt",  "CD161-adt","CD2-adt"),
  pt.size = 0.01,
  max.cutoff = 'q99',
  min.cutoff = 'q5',
  ncol = 4,
  cols = colorscale,
  cells = cells_plot,
  coord.fixed = TRUE
)&NoLegend()&NoAxes()
FeaturePlot(
  object = integrated_scRNA,
  features = c("CD56-adt", "CD62L-adt", "CD16-adt", "CD57-adt", "NKp30-adt",  "Anti-Biotin-NKG2A-adt",  "CD161-adt","CD2-adt"),
  pt.size = 0.01,
  max.cutoff = 'q99',
  min.cutoff = 'q5',
  ncol = 4,
  cols = colorscale,
  cells = cells_plot,
  coord.fixed = TRUE
)&NoAxes()
dev.off()


# DEGs: Combine markers from FindAllMarkers and Bright vs Dim comparison as otherwise these are masked
# due to similarity between adaptive and dim
DefaultAssay(integrated_scRNA) <- "RNA"

unique_genes <- FindAllMarkers(integrated_scRNA, only.pos = TRUE)

CD56dim_markers <- FindMarkers(integrated_scRNA, ident.1 = c("CD56dim"),
                               ident.2 = c("CD56bright"), only.pos = T)
CD56dim_markers$gene <- rownames(CD56dim_markers)

heatmap_markers <- unique(c(unique_genes[unique_genes$cluster%in%c("CD56bright", "EarlyCD56dim"),"gene"],
                     CD56dim_markers$gene[!CD56dim_markers$gene%in%unique_genes$gene],
                     unique_genes[unique_genes$cluster%in%c("CD56dim", "Adaptive"),"gene"]))

# Pseudo-bulk heatmap
library(tidyr)
library(tidyheatmap)
averaged_expr <- AverageExpression(integrated_scRNA, features = heatmap_markers, slot = "data", assays = "RNA")

averaged_expr <- as.data.frame(averaged_expr$RNA)

averaged_expr$gene <- rownames(averaged_expr)

averaged_expr_tidy <- gather(averaged_expr, key = "group", value = "activity", -gene)

pdf(file = "Plots/Fig1_RNA_heatmap.pdf", paper = "a4")
averaged_expr_tidy  %>% dplyr::filter(group != "Proliferating") %>%  tidy_heatmap(
  rows = gene,
  columns = group,
  values = activity,
  cluster_rows = F,
  cluster_cols = F,
  scale = "row",
  colors =  solarExtra)
dev.off()


levels(integrated_scRNA) <- rev(levels(integrated_scRNA))
pdf(file = "Plots/Fig1_scRNA_dotplot.pdf", paper = "a4r")
DotPlot(integrated_scRNA, features = c("GZMK", "XCL1", "IL7R", "TCF7",  "GPR183",
                                      "GZMB",  "PRF1",  "CX3CR1", "CD7", "FCER1G", "KLRB1",
                                      "KLRC2", "CD3E",  "PATL2", "ZBTB38"),
        dot.min = 0.01, idents = c("CD56bright", "EarlyCD56dim", "CD56dim", "Adaptive"))+scale_color_gradientn(colours = solarExtra)+theme(axis.text.x = element_text(angle = 45))
dev.off()
levels(integrated_scRNA) <- rev(levels(integrated_scRNA))

# Proliferating cluster
VlnPlot(integrated_scRNA, features = c("MKI67", "STMN1"), pt.size = 0, stack = T, fill.by ="ident" , cols = rev(c(clusterpal[1:3], adaptivepal[1], clusterpal[4])))&NoLegend()

#### Figure 2: CMVpos vs CMVneg ####

# Downsample CMVpos for plotting of NKG2C
cells_plot1 <- sample(Cells(CMVpos_scRNA), size = 10000)
cells_plot2 <- sample(Cells(CMVneg_scRNA), size = 5000)

pdf(file = "Plots/Fig2_scRNA_Dimplot_serostatus.pdf", paper = "a4")
DimPlot(CMVpos_scRNA, cols = c(clusterpal[1:4], adaptivepal))&NoLegend()&
  scale_y_continuous(limits = c(-5, 5), breaks = c(-4,-2,0,2,4,6))&
  scale_x_continuous(limits = c(-8, 6))
DimPlot(CMVneg_scRNA, cols = c(clusterpal[1:4], adaptivepal))&NoLegend()&
  scale_y_continuous(limits = c(-4.2, 6.2), breaks = c(-4,-2,0,2,4,6))

DimPlot(CMVpos_scRNA, group.by = "population", cols = c("grey",colorscale[7]), shuffle = TRUE, cells = cells_plot1)&NoLegend()&
  scale_y_continuous(limits = c(-5, 5), breaks = c(-4,-2,0,2,4,6))&
  scale_x_continuous(limits = c(-8, 6))&ggtitle(label = NULL)
DimPlot(CMVneg_scRNA, group.by = "population", cols = c("grey",colorscale[7]), shuffle = TRUE, pt.size = 0.1, cells = cells_plot2)&NoLegend()&
  scale_y_continuous(limits = c(-4.2, 6.2), breaks = c(-4,-2,0,2,4,6))&
  scale_x_continuous(limits = c(-8.5, 6))&ggtitle(label = NULL)

dev.off()

# DEGs between NKG2C+/-
# Downsample to equal numbers to enable fair comparison
table(CMVneg_scRNA$population)
table(CMVpos_scRNA$population)


downsampled_CMVpos <- c(sample(names(which(CMVpos_scRNA$population == "NKG2Cpos")), size = 6000),
                        sample(names(which(CMVpos_scRNA$population == "NKG2Cneg")), size = 6000))

downsampled_CMVneg <- c(sample(names(which(CMVneg_scRNA$population == "NKG2Cpos")), size = 6000),
                        sample(names(which(CMVneg_scRNA$population == "NKG2Cneg")), size = 6000))

# Export for reproducibility


CMVpos_down <- subset(CMVpos_scRNA, cells = downsampled_CMVpos)
CMVneg_down <- subset(CMVneg_scRNA, cells = downsampled_CMVneg)

Idents(CMVpos_down) <- CMVpos_down$population
Idents(CMVneg_down) <- CMVneg_down$population

# Include only genes that are detected with at least 100 counts in both donor groups
# to make the comparison fair
genes_tested <- intersect(names(which(rowSums(CMVpos_scRNA@assays$RNA@counts)>100)),
                          names(which(rowSums(CMVneg_scRNA@assays$RNA@counts)>100)))

markers_NKG2C_CMVpos <- FindMarkers(CMVpos_down, ident.1 = "NKG2Cpos", ident.2 = "NKG2Cneg", logfc.threshold = 0, features = genes_tested)
markers_NKG2C_CMVneg <- FindMarkers(CMVneg_down, ident.1 = "NKG2Cpos", ident.2 = "NKG2Cneg", logfc.threshold = 0, features = genes_tested)

# Add Diffexpr column
markers_NKG2C_CMVpos$diffexp <- "NO"
markers_NKG2C_CMVpos$diffexp[markers_NKG2C_CMVpos$p_val_adj < 0.05 & markers_NKG2C_CMVpos$avg_log2FC > 0.25] <- "UP"
markers_NKG2C_CMVpos$diffexp[markers_NKG2C_CMVpos$p_val_adj < 0.05 & markers_NKG2C_CMVpos$avg_log2FC < -0.25] <- "DOWN"

markers_NKG2C_CMVneg$diffexp <- "NO"
markers_NKG2C_CMVneg$diffexp[markers_NKG2C_CMVneg$p_val_adj < 0.05 & markers_NKG2C_CMVneg$avg_log2FC > 0.25] <- "UP"
markers_NKG2C_CMVneg$diffexp[markers_NKG2C_CMVneg$p_val_adj < 0.05 & markers_NKG2C_CMVneg$avg_log2FC < -0.25] <- "DOWN"

# Layer so that significant genes are on top
layer2_CMVpos <- markers_NKG2C_CMVpos[markers_NKG2C_CMVpos$diffexp%in%c("UP", "DOWN"),]
layer1_CMVpos <- markers_NKG2C_CMVneg[markers_NKG2C_CMVneg$diffexp%in%c("NO"),]

layer2_CMVneg <- markers_NKG2C_CMVneg[markers_NKG2C_CMVneg$diffexp%in%c("UP", "DOWN"),]
layer1_CMVneg <- markers_NKG2C_CMVneg[markers_NKG2C_CMVneg$diffexp%in%c("NO"),]


# Plot
genestolabel <- c("KLRC2", "KLRC1", "ZBTB16", "FCER1G", "CD2", "CD3E", "CD3D", "PATL2",
                  "LILRB1", "KLRB1", "IL32", "CD52", "CD7", "SYK", "SH2D1B", "IL2RB", "JAK1")
markers_NKG2C_CMVpos$gene <- rownames(markers_NKG2C_CMVpos)
ggplot()+
  geom_point(data = layer1_CMVpos, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  geom_point(data = layer2_CMVpos, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = markers_NKG2C_CMVpos[markers_NKG2C_CMVpos$gene%in%c(genestolabel)&!markers_NKG2C_CMVpos$diffexp=="NO",],
                  size = 5, box.padding = 0.5, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 301))+
  scale_x_continuous(limits = c(-2.65, 2.65))+
  xlab(label = "NKG2Cpos vs NKG2Cneg (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

markers_NKG2C_CMVneg$gene <- rownames(markers_NKG2C_CMVneg)
ggplot()+
  geom_point(data = layer1_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  geom_point(data = layer2_CMVneg, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), color = diffexp))+
  scale_color_manual(values = c("blue", "black", "red"))+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = markers_NKG2C_CMVneg[markers_NKG2C_CMVneg$gene%in%c(genestolabel)&!markers_NKG2C_CMVneg$diffexp=="NO",],
                  size = 5, box.padding = 0.7, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 301))+
  scale_x_continuous(limits = c(-2.65, 2.65))+
  xlab(label = "NKG2Cneg vs NKG2Cneg (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()





# No of significant genes
sum(markers_NKG2C_CMVpos$p_val_adj<0.05 & abs(markers_NKG2C_CMVpos$avg_log2FC)>0.25)
sum(markers_NKG2C_CMVneg$p_val_adj<0.05 & abs(markers_NKG2C_CMVneg$avg_log2FC)>0.25)


# Adaptive signature
integrated_scRNA_no_prol <- subset(integrated_scRNA, idents = "Proliferating", invert = TRUE)
adaptive_markers <- FindMarkers(integrated_scRNA_no_prol, ident.1 = "Adaptive", ident.2 = "CD56dim")


# Plot some representatives separately for CMVpos/CMVneg

cells_plot3 <- sample(Cells(CMVpos_scRNA), size = 5000)
cells_plot4 <- sample(Cells(CMVneg_scRNA), size = 4000)

pdf(file = "Plots/Fig2_scRNA_representatives.pdf")
FeaturePlot(CMVpos_scRNA, features = c("CD3E", "IL32", "FCER1G", "KLRB1"), coord.fixed = T, pt.size = 0.2,
            cols =  colorscale, slot = "data", min.cutoff = "q1", max.cutoff = "q95", cells = cells_plot3)&NoLegend()&NoAxes()
FeaturePlot(CMVneg_scRNA, features = c("CD3E", "IL32", "FCER1G", "KLRB1"), coord.fixed = T,pt.size = 0.2,
            cols =  colorscale, slot = "data", min.cutoff = "q1", max.cutoff = "q95", cells = cells_plot4)&NoLegend()&NoAxes()
FeaturePlot(CMVneg_scRNA, features = c("CD3E"), coord.fixed = T,pt.size = 0.2,
            cols =  colorscale, slot = "data", min.cutoff = "q1", max.cutoff = "q95", cells = cells_plot4)&NoAxes()
dev.off()


# CD56dim/bright and adaptive signature
# Find DEGs on full dataset and plot separately for CMVpos/neg
DefaultAssay(integrated_scRNA) <- "RNA"
bright_vs_dim <- FindMarkers(integrated_scRNA, ident.1 = "CD56bright", ident.2 = "CD56dim")
earlydim_vs_bright <- FindMarkers(integrated_scRNA, ident.1 = "EarlyCD56dim", ident.2 = "CD56bright")
dim_vs_early_dim <- FindMarkers(integrated_scRNA, ident.1 = "CD56dim", ident.2 = "EarlyCD56dim")
adaptive_vs_dim <- FindMarkers(integrated_scRNA, ident.1 = "Adaptive", ident.2 = "CD56dim")
all_markers <- FindAllMarkers(integrated_scRNA, only.pos = T)
top_markers <- all_markers %>% filter(cluster != "Proliferating") %>%  group_by(cluster) %>% top_n(25, wt = avg_log2FC)

# Filter for significant hits
adaptive_vs_dim <- adaptive_vs_dim %>% dplyr::filter(p_val_adj < 0.05)
bright_vs_dim <- bright_vs_dim %>% dplyr::filter(p_val_adj < 0.05)
earlydim_vs_bright <- earlydim_vs_bright %>% filter(p_val_adj < 0.05)
dim_vs_early_dim <- dim_vs_early_dim %>% filter(p_val_adj < 0.05)

bright_vs_dim$gene <- rownames(bright_vs_dim)
bright_vs_dim %>% dplyr::filter(avg_log2FC < 0) %>% top_n(n = -25, wt= avg_log2FC) %>% pull(gene) -> dim_markers
bright_vs_dim %>% dplyr::filter(avg_log2FC > 0) %>% top_n(n = 25, wt= avg_log2FC) %>% pull(gene) -> bright_markers

adaptive_vs_dim$gene <- rownames(adaptive_vs_dim)
adaptive_vs_dim %>% dplyr::filter(avg_log2FC > 0) %>% top_n(n = 25, wt= avg_log2FC) %>% pull(gene) -> adaptive_markers
adaptive_vs_dim %>% dplyr::filter(avg_log2FC < 0) %>% top_n(n = -25, wt= avg_log2FC) %>% pull(gene) -> conventional_markers
heatmap_markers <- unique(c(bright_markers, dim_markers,conventional_markers, adaptive_markers))

# Volcanoplots of DEGs
genestolabel <- c("KLRC2", "KLRC1", "ZBTB16", "FCER1G", "CD2", "CD3E", "CD3D", "ZBTB38", "PATL2",
                  "LILRB1", "KLRB1", "IL32", "CD52", "CD7", "ZBTB16", "HLA-DRB1")
adaptive_vs_dim$gene <- rownames(adaptive_vs_dim)

adaptive_vs_dim %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = adaptive_vs_dim[adaptive_vs_dim$gene%in%c(genestolabel),],
                  size = 5, box.padding = 0.8, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 301))+
  scale_x_continuous(limits = c(-3.5, 3.5))+
  xlab(label = "Adaptive vs CD56dim (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()


genestolabel <- c("GZMK", "XCL1", "IL7R", "TCF7", "FOS", "GPR183",
                  "GZMB",  "PRF1",  "CX3CR1", "ZBTB16")
bright_vs_dim$gene <- rownames(bright_vs_dim)

bright_vs_dim %>% ggplot(aes(x = -avg_log2FC, y = -log10(p_val_adj+1e-301)))+
  geom_point()+
  geom_text_repel(aes(x = -avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = bright_vs_dim[bright_vs_dim$gene%in%c(genestolabel),],
                  size = 5, box.padding = 1, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 301))+
  scale_x_continuous(limits = c(-3.5, 3.5))+
  xlab(label = "CD56dim vs CD56bright (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

earlydim_vs_bright$gene <- rownames(earlydim_vs_bright)

earlydim_vs_bright %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = earlydim_vs_bright[earlydim_vs_bright$gene%in%c(genestolabel),],
                  size = 5, box.padding = 1, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 301))+
  scale_x_continuous(limits = c(-3.5, 3.5))+
  xlab(label = "EarlyCD56dim vs CD56bright (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()

dim_vs_early_dim $gene <- rownames(dim_vs_early_dim) 

dim_vs_early_dim %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301)))+
  geom_point()+
  geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj+1e-301), label = gene),
                  data = dim_vs_early_dim[dim_vs_early_dim $gene%in%c(genestolabel),],
                  size = 5, box.padding = 1, max.overlaps = Inf)+
  scale_y_continuous(limits = c(0, 301))+
  scale_x_continuous(limits = c(-3.5, 3.5))+
  xlab(label = "CD56dim vs EarlyCD56dim (Log2FC)")+
  ylab(label = "-log10(pvalue)")+
  theme_classic(base_size = 14)&NoLegend()




# Naive NKG2C+ from CMV+
NKG2Cpos_CMVpos <- subset(CMVpos_scRNA, subset = population == "NKG2Cpos") 
NKG2Cpos_CMVpos <- RunUMAP(NKG2Cpos_CMVpos, dims = 1:10)
DimPlot(NKG2Cpos_CMVpos, cols = c(clusterpal[1:3], adaptivepal[1], clusterpal[4]))&NoLegend()&scale_y_continuous(breaks = c(-4,-2,0,2,4))&scale_x_continuous(breaks = c(-4,-2,0,2,4,6))

NKG2C_downsampled <- sample(Cells(NKG2Cpos_CMVpos), size = 7000)
FeaturePlot(NKG2Cpos_CMVpos, features = c("CD3E", "IL32", "FCER1G", "KLRB1"), coord.fixed = T,
            cols =  colorscale, slot = "data", min.cutoff = "q1", max.cutoff = "q95", cells = NKG2C_downsampled)&NoLegend()&NoAxes()

# Plot signatures from full dataset
# Pseudo-bulk heatmap
library(tidyr)
library(tidyheatmap)
averaged_expr_2Cpos <- AverageExpression(NKG2Cpos_CMVpos, features = heatmap_markers, slot = "data", assays = "RNA")

averaged_expr_2Cpos <- as.data.frame(averaged_expr_2Cpos$RNA)

averaged_expr_2Cpos$gene <- rownames(averaged_expr_2Cpos)

averaged_expr_2Cpos_tidy <- gather(averaged_expr_2Cpos, key = "group", value = "activity", -gene)

pdf(file = "Plots/SupplFig2_RNA_heatmap.pdf", paper = "a4")
averaged_expr_2Cpos_tidy  %>% dplyr::filter(group != "Proliferating") %>%  tidy_heatmap(
  rows = gene,
  columns = group,
  values = activity,
  cluster_rows = F,
  cluster_cols = F,
  scale = "row",
  colors =  solarExtra)
dev.off()





