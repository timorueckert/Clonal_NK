# Upstream analysis of scATAC seq data of human NK cells from
# CMV+/- individuals 

rm(list = ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(SeuratWrappers)
library(motifmatchr)
library(Matrix)
library(data.table)
library(future)
library(GenomicRanges)
library(viridis)
library(dplyr)


set.seed(1234)

path <- "~/"

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- brewer.pal(8, "Dark2")
NKG2Cpal <- c("#40BAD5","#120136")

##### Load and join peaks to harmonize datasets  #####
# read in peak sets
CMVpos2_peaks <- read.table(
  file = paste0(path, "data/CMVpos2/CMVpos2/outs/peaks.bed"),
  col.names = c("chr", "start", "end")
)

CMVneg2_peaks <- read.table(
  file = paste0(path, "data/CMVneg2/CMVneg2/outs/peaks.bed"),
  col.names = c("chr", "start", "end")
)

mtASAP1_peaks <- read.table(
  file = paste0(path, "data/mtASAP1/mtASAP1_aggr/outs/peaks.bed"),
  col.names = c("chr", "start", "end")
)

mtASAP2_peaks <- read.table(
  file = paste0(path, "data/mtASAP2/mtASAP2_aggr/outs/peaks.bed"),
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
CMVpos2_peaks<- makeGRangesFromDataFrame(CMVpos2_peaks)
CMVneg2_peaks<- makeGRangesFromDataFrame(CMVneg2_peaks)
mtASAP1_peaks<- makeGRangesFromDataFrame(mtASAP1_peaks)
mtASAP2_peaks<- makeGRangesFromDataFrame(mtASAP2_peaks)

# Reduce peaks
combined.peaks <- reduce(x = c(CMVpos2_peaks, CMVneg2_peaks, mtASAP1_peaks, mtASAP2_peaks))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Read in metadata
CMVpos2_metadata <- read.csv(
  file = paste0(path, "data/CMVpos2/CMVpos2/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

CMVneg2_metadata <- read.csv(
  file = paste0(path, "data/CMVneg2/CMVneg2/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

mtASAP1_metadata <- read.csv(
  file = paste0(path, "data/mtASAP1/mtASAP1_aggr/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

mtASAP2_metadata <- read.csv(
  file = paste0(path, "data/mtASAP2/mtASAP2_aggr/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

# perform an initial filtering of low count cells
CMVpos2_metadata <- CMVpos2_metadata[CMVpos2_metadata$passed_filters > 500, ]
CMVneg2_metadata <- CMVneg2_metadata[CMVneg2_metadata$passed_filters > 500, ]
mtASAP1_metadata <- mtASAP1_metadata[mtASAP1_metadata$passed_filters > 500, ]
mtASAP2_metadata <- mtASAP2_metadata[mtASAP2_metadata$passed_filters > 500, ]

# Import fragments
CMVpos2_frag <- CreateFragmentObject(
  path = paste0(path, "data/CMVpos2/CMVpos2/outs/fragments.tsv.gz"),
  cells = rownames(CMVpos2_metadata)
)

CMVneg2_frag <- CreateFragmentObject(
  path = paste0(path, "data/CMVneg2/CMVneg2/outs/fragments.tsv.gz"),
  cells = rownames(CMVneg2_metadata)
)

mtASAP1_frag <- CreateFragmentObject(
  path = paste0(path, "data/mtASAP1/mtASAP1_aggr/outs/fragments.tsv.gz"),
  cells = rownames(mtASAP1_metadata)
)

mtASAP2_frag <- CreateFragmentObject(
  path = paste0(path, "data/mtASAP2/mtASAP2_aggr/outs/fragments.tsv.gz"),
  cells = rownames(mtASAP2_metadata)
)

# Quantify Peaks
CMVpos2.counts <- FeatureMatrix(
  fragments = CMVpos2_frag,
  features = combined.peaks,
  cells = rownames(CMVpos2_metadata)
)

CMVneg2.counts <- FeatureMatrix(
  fragments = CMVneg2_frag,
  features = combined.peaks,
  cells = rownames(CMVneg2_metadata)
)

mtASAP1.counts <- FeatureMatrix(
  fragments = mtASAP1_frag,
  features = combined.peaks,
  cells = rownames(mtASAP1_metadata)
)

mtASAP2.counts <- FeatureMatrix(
  fragments = mtASAP2_frag,
  features = combined.peaks,
  cells = rownames(mtASAP2_metadata)
)



# Connection to UCSF servers fails, seq_info was downloaded manually
seq_info <- readRDS(paste0(path, "seq_info"))


# Create Objects
CMVpos2_assay <- CreateChromatinAssay(CMVpos2.counts, fragments = CMVpos2_frag, genome = seq_info)
CMVpos2 <- CreateSeuratObject(CMVpos2_assay, assay = "ATAC", meta.data = CMVpos2_metadata)

CMVneg2_assay <- CreateChromatinAssay(CMVneg2.counts, fragments = CMVneg2_frag, genome = seq_info)
CMVneg2 <- CreateSeuratObject(CMVneg2_assay, assay = "ATAC", meta.data = CMVneg2_metadata)

mtASAP1_assay <- CreateChromatinAssay(mtASAP1.counts, fragments = mtASAP1_frag, genome = seq_info)
mtASAP1 <- CreateSeuratObject(mtASAP1_assay, assay = "ATAC", meta.data = mtASAP1_metadata)

mtASAP2_assay <- CreateChromatinAssay(mtASAP2.counts, fragments = mtASAP2_frag, genome = seq_info)
mtASAP2 <- CreateSeuratObject(mtASAP2_assay, assay = "ATAC", meta.data = mtASAP2_metadata)

Human_NK <- list(CMVpos2, CMVneg2, mtASAP1, mtASAP2)

rm(list = "CMVpos2", "CMVneg2", "mtASAP1", "mtASAP2")
rm(list = "CMVpos2_assay", "CMVneg2_assay", "mtASAP1_assay", "mtASAP2_assay")
rm(list = "CMVpos2_metadata", "CMVneg2_metadata", "mtASAP1_metadata", "mtASAP2_metadata")
rm(list = "CMVpos2_frag", "CMVneg2_frag", "mtASAP1_frag", "mtASAP2_frag")

Human_NK[[1]]
Human_NK[[2]]
Human_NK[[3]]
Human_NK[[4]]


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Human_NK <- lapply(Human_NK, function(x){
  Annotation(x) <- annotations
  return(x)
})

#### QC ####

Human_NK <- lapply(Human_NK, function(x){
  x <- NucleosomeSignal(object = x)
  return(x)
})

Human_NK <- lapply(Human_NK, function(x){
  x <- TSSEnrichment(object = x, fast = FALSE)
  return(x)
})


# add blacklist ratio and fraction of reads in peaks
Human_NK <- lapply(Human_NK, function(x){
  x$pct_reads_in_peaks <- x$peak_region_fragments / x$passed_filters * 100
  x$blacklist_ratio <- x$blacklist_region_fragments / x$peak_region_fragments
  return(x)}
  )

# Plot QC metrics
lapply(Human_NK, function(x){
  VlnPlot(
    object =x,
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  )
})


# Remove outliers
Human_NK[[1]] <- subset(
  x = Human_NK[[1]],
  subset = peak_region_fragments < 13000 &
    peak_region_fragments > 2000 &
    pct_reads_in_peaks > 60 &
    blacklist_ratio < 0.0001 &
    nucleosome_signal < 1 &
    TSS.enrichment > 3
)

Human_NK[[2]] <- subset(
  x = Human_NK[[2]],
  subset = peak_region_fragments < 13000 &
    peak_region_fragments > 2000 &
    pct_reads_in_peaks > 70 &
    blacklist_ratio < 0.0001 &
    nucleosome_signal < 1 &
    TSS.enrichment > 4
)


Human_NK[[3]] <- subset(
  x = Human_NK[[3]],
  subset = peak_region_fragments < 10000 &
    peak_region_fragments > 500 &
    pct_reads_in_peaks > 55 &
    blacklist_ratio < 0.0001 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 3
)

Human_NK[[4]] <- subset(
  x = Human_NK[[4]], 
  subset = peak_region_fragments < 11000 &
    peak_region_fragments > 500 &
    pct_reads_in_peaks > 60 &
    blacklist_ratio < 0.0001 &
    nucleosome_signal < 1.1 &
    TSS.enrichment > 2.5)

##### Addition of CITE-Seq data to ASAP datasets #####
# Import data
import_kite_counts <- function(path){
  mtx <- fread(paste0(path,"featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(path,"featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- (fread(paste0(path,"featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

mtASAP1_LibA_HTO <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/HTO_A_featurecounts/"))
mtASAP1_LibB_HTO <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/HTO_B_featurecounts/"))
mtASAP1_LibC_HTO <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/HTO_C_featurecounts/"))
mtASAP1_LibD_HTO <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/HTO_D_featurecounts/"))

mtASAP1_LibA_ADT <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/ADT_A_featurecounts/"))
mtASAP1_LibB_ADT <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/ADT_B_featurecounts/"))
mtASAP1_LibC_ADT <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/ADT_C_featurecounts/"))
mtASAP1_LibD_ADT <- import_kite_counts(paste0(path, "data/mtASAP1/mtASAP1_CITE/ADT_D_featurecounts/"))


mtASAP2_LibA_HTO <- import_kite_counts(paste0(path, "data/mtASAP2/mtASAP2_CITE/HTO_A2_featurecounts/"))
mtASAP2_LibB_HTO <- import_kite_counts(paste0(path, "data/mtASAP2/mtASAP2_CITE/HTO_B2_featurecounts/"))

mtASAP2_LibA_ADT <- import_kite_counts(paste0(path, "data/mtASAP2/mtASAP2_CITE/ADT_A2_featurecounts/"))
mtASAP2_LibB_ADT <- import_kite_counts(paste0(path, "data/mtASAP2/mtASAP2_CITE/ADT_B2_featurecounts/"))

# Harmonize cellnames
colnames(mtASAP1_LibA_HTO) <- paste0(colnames(mtASAP1_LibA_HTO), "-1")
table(colnames(mtASAP1_LibA_HTO)%in%colnames(Human_NK[[3]]))
mtASAP1_LibA_HTO <- mtASAP1_LibA_HTO[,colnames(mtASAP1_LibA_HTO)%in%colnames(Human_NK[[3]])]
colnames(mtASAP1_LibA_ADT) <- paste0(colnames(mtASAP1_LibA_ADT), "-1")
table(colnames(mtASAP1_LibA_ADT)%in%colnames(Human_NK[[3]]))
mtASAP1_LibA_ADT <- mtASAP1_LibA_ADT[,colnames(mtASAP1_LibA_ADT)%in%colnames(Human_NK[[3]])]

colnames(mtASAP1_LibB_HTO) <- paste0(colnames(mtASAP1_LibB_HTO), "-2")
table(colnames(mtASAP1_LibB_HTO)%in%colnames(Human_NK[[3]]))
mtASAP1_LibB_HTO <- mtASAP1_LibB_HTO[,colnames(mtASAP1_LibB_HTO)%in%colnames(Human_NK[[3]])]
colnames(mtASAP1_LibB_ADT) <- paste0(colnames(mtASAP1_LibB_ADT), "-2")
table(colnames(mtASAP1_LibB_ADT)%in%colnames(Human_NK[[3]]))
mtASAP1_LibB_ADT <- mtASAP1_LibB_ADT[,colnames(mtASAP1_LibB_ADT)%in%colnames(Human_NK[[3]])]

colnames(mtASAP1_LibC_HTO) <- paste0(colnames(mtASAP1_LibC_HTO), "-3")
table(colnames(mtASAP1_LibC_HTO)%in%colnames(Human_NK[[3]]))
mtASAP1_LibC_HTO <- mtASAP1_LibC_HTO[,colnames(mtASAP1_LibC_HTO)%in%colnames(Human_NK[[3]])]
colnames(mtASAP1_LibC_ADT) <- paste0(colnames(mtASAP1_LibC_ADT), "-3")
table(colnames(mtASAP1_LibC_ADT)%in%colnames(Human_NK[[3]]))
mtASAP1_LibC_ADT <- mtASAP1_LibC_ADT[,colnames(mtASAP1_LibC_ADT)%in%colnames(Human_NK[[3]])]

colnames(mtASAP1_LibD_HTO) <- paste0(colnames(mtASAP1_LibD_HTO), "-4")
table(colnames(mtASAP1_LibD_HTO)%in%colnames(Human_NK[[3]]))
mtASAP1_LibD_HTO <- mtASAP1_LibD_HTO[,colnames(mtASAP1_LibD_HTO)%in%colnames(Human_NK[[3]])]
colnames(mtASAP1_LibD_ADT) <- paste0(colnames(mtASAP1_LibD_ADT), "-4")
table(colnames(mtASAP1_LibD_ADT)%in%colnames(Human_NK[[3]]))
mtASAP1_LibD_ADT <- mtASAP1_LibD_ADT[,colnames(mtASAP1_LibD_ADT)%in%colnames(Human_NK[[3]])]

colnames(mtASAP2_LibA_HTO) <- paste0(colnames(mtASAP2_LibA_HTO), "-1")
table(colnames(mtASAP2_LibA_HTO)%in%colnames(Human_NK[[4]]))
mtASAP2_LibA_HTO <- mtASAP2_LibA_HTO[,colnames(mtASAP2_LibA_HTO)%in%colnames(Human_NK[[4]])]
colnames(mtASAP2_LibA_ADT) <- paste0(colnames(mtASAP2_LibA_ADT), "-1")
table(colnames(mtASAP2_LibA_ADT)%in%colnames(Human_NK[[4]]))
mtASAP2_LibA_ADT <- mtASAP2_LibA_ADT[,colnames(mtASAP2_LibA_ADT)%in%colnames(Human_NK[[4]])]

colnames(mtASAP2_LibB_HTO) <- paste0(colnames(mtASAP2_LibB_HTO), "-2")
table(colnames(mtASAP2_LibB_HTO)%in%colnames(Human_NK[[4]]))
mtASAP2_LibB_HTO <- mtASAP2_LibB_HTO[,colnames(mtASAP2_LibB_HTO)%in%colnames(Human_NK[[4]])]
colnames(mtASAP2_LibB_ADT) <- paste0(colnames(mtASAP2_LibB_ADT), "-2")
table(colnames(mtASAP2_LibB_ADT)%in%colnames(Human_NK[[4]]))
mtASAP2_LibB_ADT <- mtASAP2_LibB_ADT[,colnames(mtASAP2_LibB_ADT)%in%colnames(Human_NK[[4]])]

# One sampe of mtASAP2 was mislabelled during counting, CMVneg4 is in fact CMVpos4 => correct
row.names(mtASAP2_LibA_HTO) <- c("CMVpos3", "CMVpos4")
row.names(mtASAP2_LibB_HTO) <- c("CMVpos3", "CMVpos4")


# Build assays
mtASAP1_LibA_HTO <- CreateAssayObject(counts = mtASAP1_LibA_HTO, min.features = 1)
mtASAP1_LibB_HTO <- CreateAssayObject(counts = mtASAP1_LibB_HTO, min.features = 1)
mtASAP1_LibC_HTO <- CreateAssayObject(counts = mtASAP1_LibC_HTO, min.features = 1)
mtASAP1_LibD_HTO <- CreateAssayObject(counts = mtASAP1_LibD_HTO, min.features = 1)

mtASAP2_LibA_HTO <- CreateAssayObject(counts = mtASAP2_LibA_HTO, min.features = 1)
mtASAP2_LibB_HTO <- CreateAssayObject(counts = mtASAP2_LibB_HTO, min.features = 1)



# Merge
mtASAP1_HTO <- CreateAssayObject(counts= cbind(mtASAP1_LibA_HTO@counts, mtASAP1_LibB_HTO@counts,
                                               mtASAP1_LibC_HTO@counts,mtASAP1_LibD_HTO@counts))
mtASAP2_HTO <- CreateAssayObject(counts= cbind(mtASAP2_LibA_HTO@counts, mtASAP2_LibB_HTO@counts))

mtASAP1_HTO <- NormalizeData(mtASAP1_HTO, normalization.method = "CLR")
mtASAP2_HTO <- NormalizeData(mtASAP2_HTO, normalization.method = "CLR")

# Include only cells which are present in both modalities
Human_NK[[3]] <- subset(Human_NK[[3]], cells = Cells(mtASAP1_HTO))
Human_NK[[4]] <- subset(Human_NK[[4]], cells = Cells(mtASAP2_HTO))


# Add HTO to Objects
Human_NK[[3]][["HTO"]] <- mtASAP1_HTO
Human_NK[[4]][["HTO"]] <- mtASAP2_HTO

# Demultiplex
Human_NK[[3]] <- HTODemux(Human_NK[[3]], assay = "HTO", positive.quantile = 0.95)
table(Human_NK[[3]]$HTO_classification.global)
Human_NK[[3]] <- subset(Human_NK[[3]], HTO_classification.global == "Singlet")

Human_NK[[4]] <- HTODemux(Human_NK[[4]], assay = "HTO", positive.quantile = 0.95)
table(Human_NK[[4]]$HTO_classification.global)
Human_NK[[4]] <- subset(Human_NK[[4]], HTO_classification.global == "Singlet")

table(Human_NK[[3]]$HTO_maxID)
table(Human_NK[[4]]$HTO_maxID)

RidgePlot(Human_NK[[3]], assay = "HTO", 
          features = rownames(Human_NK[[3]][["HTO"]]), ncol = 2)

RidgePlot(Human_NK[[4]], assay = "HTO", 
          features = rownames(Human_NK[[4]][["HTO"]]))

# Add ADT data
mtASAP1_LibA_ADT <- CreateAssayObject(counts = mtASAP1_LibA_ADT)
mtASAP1_LibB_ADT <- CreateAssayObject(counts = mtASAP1_LibB_ADT)
mtASAP1_LibC_ADT <- CreateAssayObject(counts = mtASAP1_LibC_ADT)
mtASAP1_LibD_ADT <- CreateAssayObject(counts = mtASAP1_LibD_ADT)

mtASAP2_LibA_ADT <- CreateAssayObject(counts = mtASAP2_LibA_ADT)
mtASAP2_LibB_ADT <- CreateAssayObject(counts = mtASAP2_LibB_ADT)

# Normalize per Library
mtASAP1_LibA_ADT <- NormalizeData(mtASAP1_LibA_ADT, normalization.method = "CLR")
mtASAP1_LibB_ADT <- NormalizeData(mtASAP1_LibB_ADT, normalization.method = "CLR")
mtASAP1_LibC_ADT <- NormalizeData(mtASAP1_LibC_ADT, normalization.method = "CLR")
mtASAP1_LibD_ADT <- NormalizeData(mtASAP1_LibD_ADT, normalization.method = "CLR")

mtASAP2_LibA_ADT <- NormalizeData(mtASAP2_LibA_ADT, normalization.method = "CLR")
mtASAP2_LibB_ADT <- NormalizeData(mtASAP2_LibB_ADT, normalization.method = "CLR")

# Merge, scale and add to Objects
mtASAP1_ADT <- CreateAssayObject(data=cbind(mtASAP1_LibA_ADT@data, mtASAP1_LibB_ADT@data,
                                            mtASAP1_LibC_ADT@data, mtASAP1_LibD_ADT@data))

mtASAP2_ADT <- CreateAssayObject(data=cbind(mtASAP2_LibA_ADT@data, mtASAP2_LibB_ADT@data))

mtASAP1_ADT <- ScaleData(mtASAP1_ADT)
mtASAP2_ADT <- ScaleData(mtASAP2_ADT)

mtASAP1_ADT <- subset(mtASAP1_ADT, cells= Cells(Human_NK[[3]]))
mtASAP2_ADT <- subset(mtASAP2_ADT, cells= Cells(Human_NK[[4]]))

Human_NK[[3]][["ADT"]] <- mtASAP1_ADT
Human_NK[[4]][["ADT"]] <- mtASAP2_ADT


##### Addition of mitochondrial data from mgatk output #####
# mtASAP1 had almost no mitochondrial recovery, so we only import
# mgatk output of the rest

# load mgatk output
CMVpos2_mito.data <- ReadMGATK(dir = paste0(path, "data/CMVpos2/CMVpos2_mgatk/final"))

CMVneg2_mito.data1 <- ReadMGATK(dir = paste0(path, "data/CMVneg2/CMVneg2_1_mgatk/final"))
CMVneg2_mito.data2 <- ReadMGATK(dir = paste0(path, "data/CMVneg2/CMVneg2_2_mgatk/final"))


mtASAP2_LibA_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP2/mtASAP2_LibA_mgatk/final"))
mtASAP2_LibB_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP2/mtASAP2_LibB_mgatk/final"))

# Export refallele for later analysis
mito.refallele <- CMVpos2_mito.data$refallele
saveRDS(mito.refallele, file = paste0(path, "mito.refallele"))

# Correct cellname suffixes
CMVneg2_mito.data2$counts@Dimnames[[2]] <- gsub("1", "2", CMVneg2_mito.data2$counts@Dimnames[[2]])
rownames(CMVneg2_mito.data2$depth) <- gsub("1", "2", rownames(CMVneg2_mito.data2$depth))


mtASAP2_LibB_mito.data$counts@Dimnames[[2]] <- gsub("1", "2", mtASAP2_LibB_mito.data$counts@Dimnames[[2]])
rownames(mtASAP2_LibB_mito.data$depth) <- gsub("1", "2", rownames(mtASAP2_LibB_mito.data$depth))


# Merge
CMVneg2_mitocounts_aggr <- cbind(CMVneg2_mito.data1$counts, CMVneg2_mito.data2$counts)
CMVneg2_mitodepth_aggr <- rbind(CMVneg2_mito.data1$depth, CMVneg2_mito.data2$depth)


mtASAP2_mitocounts_aggr <- cbind(mtASAP2_LibA_mito.data$counts, mtASAP2_LibB_mito.data$counts)
mtASAP2_mitodepth_aggr <- rbind(mtASAP2_LibA_mito.data$depth, mtASAP2_LibB_mito.data$depth)


# Create Assays, subset to cells present in both objects
CMVpos2_mito <- CreateAssayObject(counts = CMVpos2_mito.data$counts)
CMVneg2_mito <- CreateAssayObject(counts = CMVneg2_mitocounts_aggr)
mtASAP2_mito <- CreateAssayObject(counts = mtASAP2_mitocounts_aggr)


CMVpos2_mito <- subset(CMVpos2_mito, cells = Cells(Human_NK[[1]])[Cells(Human_NK[[1]])%in%Cells(CMVpos2_mito)])
CMVneg2_mito <- subset(CMVneg2_mito, cells = Cells(Human_NK[[2]])[Cells(Human_NK[[2]])%in%Cells(CMVneg2_mito)])
mtASAP2_mito <- subset(mtASAP2_mito, cells = Cells(Human_NK[[4]])[Cells(Human_NK[[4]])%in%Cells(mtASAP2_mito)])

Human_NK[[1]] <- subset(Human_NK[[1]], cells = Cells(CMVpos2_mito))
Human_NK[[2]] <- subset(Human_NK[[2]], cells = Cells(CMVneg2_mito))
Human_NK[[4]] <- subset(Human_NK[[4]], cells = Cells(mtASAP2_mito))


# add assay and metadata to the seurat object
Human_NK[[1]][["mito"]] <- CMVpos2_mito
Human_NK[[2]][["mito"]] <- CMVneg2_mito
Human_NK[[4]][["mito"]] <- mtASAP2_mito

Human_NK[[1]] <- AddMetaData(Human_NK[[1]], metadata = (CMVpos2_mito.data$depth), col.name = "mtDNA_depth")
Human_NK[[2]] <- AddMetaData(Human_NK[[2]], metadata = (CMVneg2_mitodepth_aggr), col.name = "mtDNA_depth")
Human_NK[[4]] <- AddMetaData(Human_NK[[4]], metadata = (mtASAP2_mitodepth_aggr), col.name = "mtDNA_depth")

# Add Experiment ID to dataset
Human_NK[[1]]$experiment <- "CMVpos2"
Human_NK[[2]]$experiment <- "CMVneg2"
Human_NK[[3]]$experiment <- "mtASAP1"
Human_NK[[4]]$experiment <- "mtASAP2"

# Add Donor ID to dataset
Human_NK[[1]]$donor <- "CMVpos2"
Human_NK[[2]]$donor <- "CMVneg2"
Human_NK[[3]]$donor <- Human_NK[[3]]$HTO_maxID
Human_NK[[4]]$donor <- Human_NK[[4]]$HTO_maxID

# Add CMVstatus to dataset
Human_NK[[1]]$serostatus <- "CMVpos"
Human_NK[[2]]$serostatus <- "CMVneg"
Human_NK[[4]]$serostatus <- "CMVpos"

mtASAP1_serostatus <- character(length=ncol(Human_NK[[3]]))
mtASAP1_serostatus[Human_NK[[3]]$HTO_maxID%in%c("CMVpos1", "CMVpos3", "CMVpos4")] <- "CMVpos"
mtASAP1_serostatus[Human_NK[[3]]$HTO_maxID%in%c("CMVneg1")] <- "CMVneg"
Human_NK[[3]]$serostatus <- mtASAP1_serostatus


# Export Full Datasets
saveRDS(Human_NK[[1]], file = "CMVpos2")
saveRDS(Human_NK[[2]], file = "CMVneg2")
saveRDS(Human_NK[[3]], file = "mtASAP1")
saveRDS(Human_NK[[4]], file = "mtASAP2")

# Re-import
CMVpos2 <- readRDS(paste0(path, "CMVpos2"))
CMVneg2 <- readRDS(paste0(path, "CMVneg2"))
mtASAP1 <- readRDS(paste0(path, "mtASAP1"))
mtASAP2 <- readRDS(paste0(path, "mtASAP2"))

Human_NK <- list(CMVpos2, CMVneg2, mtASAP1, mtASAP2)
rm(CMVpos2, CMVneg2, mtASAP1, mtASAP2)

# Normalize each dataset individually
Human_NK[[1]] <- RunTFIDF(Human_NK[[1]])
Human_NK[[2]] <- RunTFIDF(Human_NK[[2]])
Human_NK[[3]] <- RunTFIDF(Human_NK[[3]])
Human_NK[[4]] <- RunTFIDF(Human_NK[[4]])

# Find Top Features
Human_NK <- lapply(Human_NK, function(x){
  x <- FindTopFeatures(x, min.cutoff = '10')
  return(x)
})


# Merge datasets
Human_NK_merged <- merge(
  x = Human_NK[[1]],
  y = list(Human_NK[[2]], Human_NK[[3]], Human_NK[[4]]),
  add.cell.ids = c("CMVpos2", "CMVneg2", "mtASAP1", "mtASAP2")
)

Human_NK_merged
Human_NK_merged[["ATAC"]]





# Dim-Reduction
Human_NK_merged <- FindTopFeatures(Human_NK_merged, min.cutoff = '10')
Human_NK_merged <- RunSVD(Human_NK_merged)

# Check for correlation with sequencing depth
DepthCor(Human_NK_merged)

Human_NK_merged <- RunUMAP(object = Human_NK_merged, reduction = 'lsi', dims = 2:10)

DimPlot(Human_NK_merged, group.by = "donor")
DimPlot(Human_NK_merged, group.by = "serostatus")
DimPlot(Human_NK_merged, group.by = "experiment")

# Cluster at high resolution to Exclude ILCs
Human_NK_merged <- FindNeighbors(object = Human_NK_merged, dims = 2:10, reduction = 'lsi')
Human_NK_merged <- FindClusters(object = Human_NK_merged, verbose = FALSE, algorithm = 3, resolution = 1)

DefaultAssay(Human_NK_merged) <- "ADT"
FeaturePlot(Human_NK_merged, features = c("CD127", "CD117", "CD56", "CD16"), cols = colorscale, order = TRUE)
DimPlot(Human_NK_merged, label = TRUE) #Cluster 16 are ILCs

Human_NK_merged <- subset(Human_NK_merged, idents = "16", invert = TRUE)

# Plot QC of final merged dataset
VlnPlot(
  object =Human_NK_merged,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5,
  group.by = "orig.ident"
)





## Quick integration with Harmony to perform cluster-specific 
# peak calling for the merged dataset
library(harmony)
DefaultAssay(Human_NK_merged) <- "ATAC"


Human_NK_merged <- RunHarmony(
  object = Human_NK_merged,
  group.by.vars = 'experiment',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)



# re-compute the UMAP using corrected LSI embeddings
Human_NK_merged <- RunUMAP(Human_NK_merged, dims = 2:10, reduction = 'harmony', reduction.key = "hUMAP", reduction.name = "hUMAP")

p1 <- DimPlot(Human_NK_merged, group.by = "donor", reduction = "umap", shuffle = T, pt.size = 0.1)
p2 <- DimPlot(Human_NK_merged, group.by = "donor", reduction = "hUMAP", shuffle = T, pt.size = 0.1)
p3 <- DimPlot(Human_NK_merged, group.by = "experiment", reduction = "umap", cols = clusterpal, shuffle = T, pt.size = 0.1)
p4 <- DimPlot(Human_NK_merged, group.by = "experiment", reduction = "hUMAP", cols = clusterpal, shuffle = T, pt.size = 0.1)
p5 <- DimPlot(Human_NK_merged, group.by = "serostatus", reduction = "umap", cols = clusterpal, shuffle = T, pt.size = 0.1)
p6 <- DimPlot(Human_NK_merged, group.by = "serostatus", reduction = "hUMAP", cols = clusterpal, shuffle = T, pt.size = 0.1)

p1 + p2 + p3 + p4 + p5 + p6

#### Clustering
Human_NK_merged <- FindNeighbors(object = Human_NK_merged, reduction = 'harmony', dims = 2:10)
Human_NK_merged <- FindClusters(object = Human_NK_merged, verbose = FALSE, algorithm = 3, resolution = 0.2)
DimPlot(Human_NK_merged, reduction = "hUMAP")

# Merge adaptive clusters for peak calling to avoid donor specific effects
Human_NK_merged <- RenameIdents(Human_NK_merged, "2" = "1")



# Export merged object
saveRDS(Human_NK_merged, file = "Human_NK_merged")

Human_NK_merged <- readRDS("Human_NK_merged")

# Stash ident as metadata
Human_NK_merged$coarse <- Idents(Human_NK_merged)

# Call peaks on coarse clusters to recover fine peak resolution lost during reduction for merging
# Downsample to 10000 cells for peak calling to avoid calling of very low represented peaks
Downsampled <- subset(Human_NK_merged, cells = sample(Cells(Human_NK_merged), size = 10000))

DimPlot(Human_NK_merged)

peaks <- CallPeaks(Human_NK_merged,  macs2.path = "/opt/anaconda3/bin/macs2", group.by = "coarse")
saveRDS(peaks, file = "peaks_per_cluster")
peaks <- readRDS("peaks_per_cluster")


CoveragePlot(
  object = Human_NK_merged,
  region = "BCL11B",
  extend.upstream = 10000,
  extend.downstream = 10000,
  links = F,
  ranges = temp
)

New_counts <- FeatureMatrix(
  fragments = Fragments(Human_NK_merged),
  features = peaks,
  cells = colnames(Human_NK_merged)
)


MACS2_peaks <- CreateChromatinAssay(New_counts, fragments = Fragments(Human_NK_merged), genome = seq_info)

Human_NK_merged[["MACS2"]] <- MACS2_peaks 
DefaultAssay(Human_NK_merged) <- "MACS2"

# Remove old assay
Human_NK_merged[["ATAC"]] <- NULL

# Normalize new assay
Human_NK_merged <- RunTFIDF(Human_NK_merged)

# Dim-Reduction
Human_NK_merged <- FindTopFeatures(Human_NK_merged, min.cutoff = '20')
Human_NK_merged <- RunSVD(Human_NK_merged)


# Check for correlation with sequencing depth
DepthCor(Human_NK_merged)

Human_NK_merged <- RunUMAP(object = Human_NK_merged, reduction = 'lsi', dims = 2:10)

DimPlot(Human_NK_merged, group.by = "donor")
DimPlot(Human_NK_merged, group.by = "serostatus")
DimPlot(Human_NK_merged, group.by = "experiment")


##### Anchor-based integration of donors #####
# Split donors up
CMVpos1 <- subset(Human_NK_merged, donor == "CMVpos1")
CMVpos2 <- subset(Human_NK_merged, donor == "CMVpos2")
CMVpos3 <- subset(Human_NK_merged, donor == "CMVpos3")
CMVpos4 <- subset(Human_NK_merged, donor == "CMVpos4")
CMVneg1 <- subset(Human_NK_merged, donor == "CMVneg1")
CMVneg2 <- subset(Human_NK_merged, donor == "CMVneg2")

Single_donors <- list(CMVpos1, CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)
rm(CMVpos1, CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)

# Process single donors
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

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = Human_NK_merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)


integrated <- readRDS("integrated")

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:10,n.neighbors = 40, return.model = TRUE)
p1 <- DimPlot(integrated, group.by = "donor")
p2 <- DimPlot(integrated, group.by = "experiment")
p3 <- DimPlot(integrated, group.by = "serostatus")
p1&p2&p3

#### Clustering
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:10, k.param = 40)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = 0.2)



# Inspect CITE-Seq
DefaultAssay(integrated) <- "ADT"

cells_plot <- Cells(subset(integrated, experiment%in%c("mtASAP1", "mtASAP2")))

cells_downsampled <- sample(cells_plot, size = 20000)

p1 <- FeaturePlot(
  object = integrated,
  features = c("CD56", "CD16", "NKp30", "CD57", "anti-Biotin", "CD2"),
  pt.size = 0.01,
  max.cutoff = 'q99',
  min.cutoff = 'q1',
  ncol = 3,
  cols = colorscale,
  combine = FALSE,
  #reduction = "hUMAP",
  order = TRUE,
  cells = cells_downsampled
)

CombinePlots(lapply(p1, function(x) x+NoLegend()+NoAxes()))

# Subcluster bright to resolve early cd56dim
temp <- subset(integrated, idents ="2")
temp <- FindNeighbors(object = temp, reduction = 'integrated_lsi', dims = 2:10, k.param = 40)
temp <- FindClusters(object = temp, verbose = FALSE, algorithm = 3, resolution = 0.2)
DimPlot(temp)
Idents(integrated, cells = WhichCells(temp, idents =0)) <- "EarlyCD56dim"
DimPlot(integrated)

# Annotate clusters
integrated <- RenameIdents(integrated, "0" = "CD56dim", "1" = "Adaptive", "2" = "CD56bright")
levels(integrated) <- c("CD56bright", "EarlyCD56dim", "CD56dim", "Adaptive")
integrated$annotation <- Idents(integrated)
  
DimPlot(integrated, cols = clusterpal)

saveRDS(integrated, "integrated")



##### Finding DA peaks #####
DefaultAssay(integrated) <- 'MACS2'
Annotation(integrated) <- Annotation(integrated[["ATAC"]])


##### Export peaks for motif analysis with Homer #####
# Generate granges from Seurat subset to peaks of interest
adaptive_peaks <- FindMarkers(integrated, ident.1 = c("Adaptive"),
                              ident.2 = c("CD56dim"), min.pct = 0.05, test.use = 'LR', logfc.threshold = 0.1,
                              latent.vars = 'nCount_MACS2', only.pos = TRUE,assay = "MACS2")


adaptive_granges <- granges(integrated[rownames(adaptive_peaks),])

# Generate dataframe with Homer-compatible layout
adaptive_peaks_df <- data.frame(
  ID = 1:length(rownames(adaptive_peaks)),
  chrm = seqnames(adaptive_granges),
  start = start(ranges(adaptive_granges)),
  end = end(ranges(adaptive_granges)),
  strand = "+"
)


write.table(adaptive_peaks_df, file = "adaptive_peaks.txt",
            quote = F, sep="\t", row.names = F, col.names = F)



##### Gene Activity Matrix: Proximity Model #####
gene.activities <- GeneActivity(integrated)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
integrated[['RNA_atac']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA_atac)
)




##### TF motif Analysis #####
DefaultAssay(integrated) <- 'MACS2'
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE, collection = "CORE")
)

# Add missing STAT motifs from mouse database
mouse_pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Mus musculus", all_versions = FALSE, collection = "CORE")
)

Stat_ids <- c("MA0518.1","MA0519.1", "MA1623.1", "MA1624.1", "MA0520.1", "MA1625.1")

names(Stat_ids) <- c("Stat4", "Stat5a::Stat5b", "Stat2", "Stat5a", "Stat6", "Stat5b")

pfm[Stat_ids] <- mouse_pfm[Stat_ids]

pfm


motif.positions <- matchMotifs(
  pwms = pfm,
  subject = granges(integrated),
  out = 'positions',
  genome = 'hg38'
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(integrated),
  pwm = pfm,
  genome = 'hg38',
  use.counts = FALSE
)

# Create a new Motif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm,
  positions = motif.positions
)

# Add the Motif object to the assay
integrated <- SetAssayData(
  object = integrated,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif
)
integrated[["MACS2"]]

# Calculate regional stats
integrated <- RegionStats(object = integrated, genome = BSgenome.Hsapiens.UCSC.hg38)


##### Calculate motif activities #####
integrated <- RunChromVAR(
  object = integrated,
  genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(integrated) <- 'chromvar'



# Add assay with motif names rather than IDs for more convenient plotting##
chromvar_names <- GetAssayData(integrated, assay = "chromvar", slot = "data")
rownames(chromvar_names) <- c(ConvertMotifID(integrated, id = rownames(integrated[["chromvar"]])[1:633], assay = "MACS2"),
                              names(Stat_ids))
integrated[["chromvar_names"]] <- CreateAssayObject(data = chromvar_names)


# Export processed object
saveRDS(integrated, file = "integrated_processed")
integrated <- readRDS("integrated_processed")


# Subset CMVpos/neg
CMVpos <- subset(integrated, donor%in%c("CMVpos1", "CMVpos2", "CMVpos3", "CMVpos4"))
CMVneg <- subset(integrated, donor%in%c("CMVneg1", "CMVneg2"))

# Re-calculate UMAP
CMVpos <- RunUMAP(CMVpos, dims = 2:10, reduction = 'integrated_lsi')
CMVneg <- RunUMAP(CMVneg, dims = 2:10, reduction = 'integrated_lsi')

# Check labels from clustering on full dataset
DimPlot(CMVpos)
DimPlot(CMVneg) 
# Cells that cluster as adaptive do not form a separate cluster here! => likely some conventional dim 
# which share part of the signature

# Re-cluster CMVneg
# CMVpos <- FindNeighbors(CMVpos, dims = 2:10, reduction = 'integrated_lsi')
# CMVpos <- FindClusters(CMVpos, resolution = 0.1, algorithm = 3)

CMVneg <- FindNeighbors(CMVneg, dims = 2:10, reduction = 'integrated_lsi')
CMVneg <- FindClusters(CMVneg, resolution = 0.1, algorithm = 3)
FeaturePlot(CMVneg, features = c("CD56", "CD16", "CD127", "CD117"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE)
# Subcluster "1" to resolve Early dim (as in original clustering)
temp <- subset(CMVneg, idents = "1")
temp <- FindNeighbors(temp, dims = 2:10, reduction = 'integrated_lsi')
temp <- FindClusters(temp, resolution = 0.1, algorithm = 3)
DimPlot(temp)
Idents(CMVneg, cells = WhichCells(temp, idents = "0")) <- "EarlyCD56dim"
CMVneg <- RenameIdents(CMVneg, "1" = "CD56bright", "0" = "CD56dim")
levels(CMVneg) <- c("CD56bright", "EarlyCD56dim", "CD56dim")
DimPlot(CMVneg)

# Export 
saveRDS(CMVpos, file = "CMVpos_processed")
saveRDS(CMVneg, file = "CMVneg_processed")


