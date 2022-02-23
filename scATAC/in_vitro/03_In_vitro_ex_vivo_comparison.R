# Projection of in vitro activated cells onto ex vivo dataset
# to analyze the signals involved in adaptive remodelling

rm(list = ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(future)


set.seed(1234)


# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
clusterpal <- brewer.pal(8, "Dark2")
NKG2Cpal <- c("#40BAD5","#120136")
pairedpal <- brewer.pal(12, "Paired")

# Import processed datasets
integrated_exvivo <- readRDS(paste0(path, "integrated_processed"))
NK_act <- readRDS(paste0(path, "NK_act_dim1"))

# Check datasets
DimPlot(integrated_exvivo)&DimPlot(NK_act)

# Quantify ex vivo peaks in in vitro dataset
DefaultAssay(integrated_exvivo ) <- "MACS2"
DefaultAssay(NK_act) <- "MACS2"

counts <- FeatureMatrix(
  fragments = Fragments(NK_act),
  features = granges(integrated_exvivo),
  cells = colnames(NK_act)
)

# create object
exvivo.assay <- CreateChromatinAssay(
  counts = counts,
  min.features = 1000,
  fragments = Fragments(NK_act)
)

NK_act_exvivo <- CreateSeuratObject(counts = exvivo.assay, assay = "MACS2")

# Add metadata from original analysis
NK_act_exvivo$donor  <- NK_act$donor
NK_act_exvivo$condition  <- NK_act$condition
NK_act_exvivo$cluster <- Idents(NK_act)
Idents(NK_act_exvivo) <- NK_act_exvivo$cluster

# compute LSI
NK_act_exvivo <- FindTopFeatures(NK_act_exvivo, min.cutoff = "q1")
NK_act_exvivo<- RunTFIDF(NK_act_exvivo)
NK_act_exvivo <- RunSVD(NK_act_exvivo)

# Calculate UMAP and look at donor distribution
NK_act_exvivo <- RunUMAP(NK_act_exvivo, dims = 2:30, reduction = "lsi")
DimPlot(NK_act_exvivo, group.by = "donor")

# Harmony integration on ex vivo assay in NK_act
library(harmony)
NK_act_exvivo <- RunHarmony(
  object = NK_act_exvivo,
  group.by.vars = 'donor',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE
)

NK_act_exvivo <- RunUMAP(object = NK_act_exvivo, reduction = 'harmony', dims = 2:30)

DimPlot(NK_act_exvivo, group.by = "donor")
DimPlot(NK_act_exvivo, group.by = "cluster")


# Add dataset-metadata
integrated_exvivo$dataset <- "ex vivo"
NK_act_exvivo$dataset <- "in vitro"

# compute UMAP and store the UMAP model
integrated_exvivo  <- RunUMAP(integrated_exvivo, reduction = "integrated_lsi", dims = 2:10, return.model = TRUE)

# Compute an uncorrected lsi for projection
integrated_exvivo <- RunTFIDF(integrated_exvivo)
integrated_exvivo <- RunSVD(integrated_exvivo)

# Create new Seurat object only containing required assay and reductions
NK_exvivo <- CreateSeuratObject(counts = integrated_exvivo[["MACS2"]], assay = "MACS2")
NK_exvivo[["lsi"]] <- integrated_exvivo[["lsi"]]
NK_exvivo[["umap"]] <- integrated_exvivo[["umap"]]
Idents(NK_exvivo) <- Idents(integrated_exvivo)

NK_exvivo
DimPlot(NK_exvivo)

# Subset to shared features
shared_features <- intersect(VariableFeatures(NK_exvivo), VariableFeatures(NK_act_exvivo))
NK_exvivo <- NK_exvivo[shared_features,]
NK_act_exvivo <- NK_act_exvivo[shared_features,]

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = NK_exvivo,
  query = NK_act_exvivo,
  reference.assay = "MACS2",
  query.assay = "MACS2",
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:10,
  verbose = TRUE,
  k.filter = NA
)

# Transfer labels from ex vivo to in vitro
label_predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = Idents(NK_exvivo),
  weight.reduction = "lsiproject"
)

# Export
saveRDS(label_predictions, paste0(path, "label_predictions"))
saveRDS(NK_act_exvivo, file = paste0(path, "NK_act_exvivo_assay"))

label_predictions <- readRDS(paste0(path, "label_predictions"))
NK_act_exvivo <- readRDS(paste0(path, "NK_act_exvivo_assay"))

# Transfer prediction score to NK_dim from main analysis
NK_act <- AddMetaData(NK_act, metadata = label_predictions)

pdf("Plots/Fig4_prediction_scores.pdf", paper = "a4r")
VlnPlot(NK_act, features = c("prediction.score.Adaptive"), cols =  pairedpal[c(1,5,2,6)], pt.size = 0)&NoLegend()
VlnPlot(NK_act, features = c("prediction.score.CD56dim"), cols =  pairedpal[c(1,5,2,6)], pt.size = 0)&NoLegend()
VlnPlot(NK_act, features = c("prediction.score.CD56bright"), cols =  pairedpal[c(1,5,2,6)], pt.size = 0.1)&NoLegend()
dev.off()


# Test for significance
label_df <- data.frame(cluster = Idents(NK_act),
                       adaptive = NK_act$prediction.score.Adaptive,
                       dim = NK_act$prediction.score.CD56dim,
                       bright = NK_act$prediction.score.CD56bright)


pairwise.wilcox.test(label_df$adaptive, label_df$cluster)
pairwise.wilcox.test(label_df$dim, label_df$cluster)



# Look for overlapping differential peaks
adaptive_peaks_ex_vivo <- FindMarkers(NK_exvivo, ident.1 = c("Adaptive"),
                                                         ident.2 = c("CD56dim"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.1,
                                                         latent.vars = 'nCount_MACS2')

Idents(NK_act_exvivo) <- NK_act_exvivo$cluster
peaks_in_vitro <- FindMarkers(NK_act_exvivo, ident.1 = c("Peptide+Cytokines"),
                                      ident.2 = c("Control"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.1,
                                      latent.vars = 'nCount_MACS2')


# Filter for statistically significant
adaptive_peaks_ex_vivo <- adaptive_peaks_ex_vivo %>% dplyr::filter(p_val_adj < 0.05)
peaks_in_vitro <- peaks_in_vitro %>% dplyr::filter(p_val_adj < 0.05)

saveRDS(adaptive_peaks_ex_vivo, file = paste0(path, "adaptive_peaks_ex_vivo"))
saveRDS(peaks_in_vitro, file = paste0(path, "peaks_in_vitro"))

adaptive_peaks_ex_vivo <- readRDS(paste0(path, "adaptive_peaks_ex_vivo"))
peaks_in_vitro <- readRDS(paste0(path, "peaks_in_vitro"))

# Plot together
adaptive_FC <- FoldChange(NK_exvivo, ident.1 = "Adaptive", ident.2 = "CD56dim",
                          features = unique(c(rownames(adaptive_peaks_ex_vivo), rownames(peaks_in_vitro))))
invitro_FC <- FoldChange(NK_act_exvivo, ident.1 = "Peptide+Cytokines", ident.2 = "Control",
                         features = unique(c(rownames(adaptive_peaks_ex_vivo), rownames(peaks_in_vitro))))

# Put into one dataframe for plotting
adaptive_FC$region <- rownames(adaptive_FC)
invitro_FC$region <- rownames(invitro_FC)

colnames(adaptive_FC)[1] <- "log2FC_Adaptive_vs_CD56dim"
colnames(invitro_FC)[1] <- "log2FC_PeptideCytokines_vs_Control"
adaptive_FC <- adaptive_FC[,c(1,4)]
invitro_FC <- invitro_FC[,c(1,4)]

combined_FC <- merge(adaptive_FC, invitro_FC, by = "region", all = TRUE)


# Add closest gene
temp <- ClosestFeature(NK_exvivo, regions = combined_FC$region)
combined_FC$closest_gene <- temp$gene_name

# HIghlight AP-1 peaks
rownames(combined_FC) <- combined_FC$region
AP1_peaks <- Peaks_to_motif(NK_exvivo, "MA1142.1", combined_FC)

# Plot, label some genes characteristic for adaptive
genes_to_label <- c("ZBTB16", "AIM2",  "IFNG")

combined_FC %>%
  ggplot(aes(x = log2FC_PeptideCytokines_vs_Control, y = log2FC_Adaptive_vs_CD56dim))+
  geom_point(color = "grey49")+
  geom_point(data = combined_FC[combined_FC$region%in%rownames(AP1_peaks),], color = "darkred")+
  coord_cartesian()+
  geom_smooth(method = "lm")+
  geom_text_repel(aes(x = log2FC_PeptideCytokines_vs_Control, y = log2FC_Adaptive_vs_CD56dim, label = closest_gene),
                  data = combined_FC[combined_FC$closest_gene%in%genes_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  xlab("LFL+IL-12+IL-18 vs Control (log2FC)")+
  ylab("Adaptive vs CD56dim (log2FC)")+
  theme_classic()

# Linear regression
lmfull <- lm(log2FC_PeptideCytokines_vs_Control~log2FC_Adaptive_vs_CD56dim, data = combined_FC)
summary(lmfull)

# Plot peaks statistically significant in both!
combined_FC <- combined_FC[combined_FC$region%in%intersect(rownames(adaptive_peaks_ex_vivo), rownames(peaks_in_vitro)),]

combined_FC %>%
  ggplot(aes(x = log2FC_PeptideCytokines_vs_Control, y = log2FC_Adaptive_vs_CD56dim))+
  geom_point(color = "grey49")+
  geom_point(data = combined_FC[combined_FC$region%in%rownames(AP1_peaks),], color = "darkred")+
  coord_cartesian()+
  geom_smooth(method = "lm")+
  geom_text_repel(aes(x = log2FC_PeptideCytokines_vs_Control, y = log2FC_Adaptive_vs_CD56dim, label = closest_gene),
                  data = combined_FC[combined_FC$closest_gene%in%genes_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  xlab("LFL+IL-12+IL-18 vs Control (log2FC)")+
  ylab("Adaptive vs CD56dim (log2FC)")+
  theme_classic()

# Linear regression
lmfull <- lm(log2FC_PeptideCytokines_vs_Control~log2FC_Adaptive_vs_CD56dim, data = combined_FC)
summary(lmfull)

# Test for motif enrichment in induced peaks shared with ex vivo signature
enriched.motifs_invitro <- FindMotifs(
  object = NK_exvivo,
  features = combined_FC %>%
    pull(region)
)

# enriched.motifs_invitro  <- enriched.motifs_invitro  %>% dplyr::filter(pvalue < 0.05)
enriched.motifs_invitro  <- enriched.motifs_invitro  %>% arrange(pvalue) %>% mutate(rank = row_number())

AP1tolabel <- c("JUNB", "BATF3", "BATF::JUN", "JUND", "BATF", "FOSL1::JUND",
                "JDP2", "FOSL1",  "FOS", "FOSL1::JUNB", "FOS::JUNB",
                "JUN::JUNB", "JUN(var.2)", "FOSL1::JUN", "FOS::JUN",
                "FOS::JUND", "FOSL2::JUN", "FOSL2::JUNB", "FOSL2::JUND",
                "FOSB::JUNB", "FOS::JUNB", "FOSB::JUNB(var.2)","FOSL2", "BACH2","MAF::NFE2" )

TFs_to_label <- c("RELA", "NFATC2", "REL", "Stat4", "FOSL2")

enriched.motifs_invitro %>%  ggplot(aes(x = rank, y = -log10(pvalue)))+
  geom_point(size = 1)+
  geom_point(color = "darkred", size = 1, data =  enriched.motifs_invitro[enriched.motifs_invitro$motif.name%in%AP1tolabel,])+
  geom_text_repel(aes(x = rank, y = -log10(pvalue), label = motif.name),
                  data = enriched.motifs_invitro[enriched.motifs_invitro$motif.name%in%TFs_to_label,],
                  size = 4, max.overlaps = Inf, box.padding = 0.7)+
  xlab(label = "Motif rank")+
  ylab(label = "-log10(pvalue)")+
  theme_classic()


