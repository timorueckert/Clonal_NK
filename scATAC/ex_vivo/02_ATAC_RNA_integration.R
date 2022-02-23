# Integration of ATAC and scRNAseq

rm(list = ls())
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(SeuratWrappers)
library(Matrix)
library(data.table)
library(future)
library(GenomicRanges)
library(viridis)
library(dplyr)
library(harmony)


set.seed(1234)
path <- "~/"
path2 <- "~/"


# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- brewer.pal(8, "Dark2")
NKG2Cpal <- c("#40BAD5","#120136")

# Import processed datasets
seu_ATAC <-readRDS(paste0(path, "integrated_processed"))
seu_RNA <- readRDS(paste0(path2, "integrated_scRNA"))

# Look at datasets before integration
DimPlot(seu_ATAC)
DimPlot(seu_RNA)

# Add technology to metadata
seu_ATAC$tech <- "ATAC"
seu_RNA$tech <- "RNA"


DefaultAssay(seu_RNA) <- "integrated"

# Calculate transfer anchors for imputation; Proliferating cells are not resolved in ATAC, exclude for integration 
transfer.anchors <- FindTransferAnchors(seu_RNA, query = seu_ATAC, features = VariableFeatures(object = seu_RNA), 
                                        reference.assay = "integrated", query.assay = "RNA_atac", reduction = "cca")



# Get data for all genes from full dataset
refdata <- GetAssayData(seu_RNA, assay = "RNA", slot = "data")

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seu_ATAC[["integrated_lsi"]],
                           dims = 2:10)


# Add imputed data to seu_ATAC object
seu_ATAC[["RNA_imputed"]] <- imputation


DefaultAssay(seu_ATAC) <- "RNA_imputed"
seu_ATAC <- FindVariableFeatures(seu_ATAC, nfeatures = 3000)

# Dotplot of integrated assay
DotPlot(seu_ATAC, features = rev(c("GZMK", "XCL1", "IL7R", "TCF7", "FOS", "GPR183",
                                           "GZMB",  "PRF1",  "CX3CR1",   "FCER1G", "KLRB1",
                                           "KLRC2", "CD3E", "IL32", "CD52", "GZMH", "PATL2", "ZBTB38")),
        cols = c(colorscale[1], colorscale[7]),dot.min = 0.01)+scale_color_gradientn(colours = cool_warm(500))+coord_flip()

FeaturePlot(seu_ATAC, features = c("TCF7"), cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", slot = "data")

# Find links between peaks and expression
# first compute the GC content for each peak
DefaultAssay(seu_ATAC) <- "RNA_imputed"
seu_ATAC <- FindVariableFeatures(seu_ATAC, nfeatures = 5000)
DefaultAssay(seu_ATAC) <- "MACS2"
seu_ATAC <- RegionStats(seu_ATAC, genome = BSgenome.Hsapiens.UCSC.hg38)


seu_ATAC <- LinkPeaks(
  object = seu_ATAC,
  peak.assay = "MACS2",
  expression.assay = "RNA_imputed",
  genes.use = rownames(seu_ATAC[["RNA_imputed"]]),
)

saveRDS(seu_ATAC, paste0(path, "integrated_processed"))


