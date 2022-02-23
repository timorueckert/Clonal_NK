# Upstream analysis of scRNA of human NK cells ex vivo: QC, normalization,
# scaling, clustering and dim-reduction

rm(list=ls())

library(dplyr)
library(parallel)
library(Seurat)
library(scater)
library(scran)
library(purrr)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(data.table)
library(stringr)

path <- "~/"

# Enable parallelization via mclapply
ncores <- detectCores()-1

# Sample names
my_samples <- c("CMVpos2", "CMVpos3", "CMVpos4", "CMVneg1", "CMVneg2", "CMVpos1_v3", "CMVpos5_v3")

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
clusterpal <- brewer.pal(8, "Dark2")
clusterpal <- c(clusterpal[1:3],clusterpal[5:8])
NKG2Cpal <- c("#40BAD5","#120136")
divpal <- brewer.pal(n = 9, name = "RdBu")
adaptivepal <- c("#880E4F",  "#C2185B",  "#E91E63",  "#F06292", "#F8BBD0", "#4A148C", "#7B1FA2", "#9C27B0", "#BA68C8")

viriscale <- viridis(9)

#colorscale <- plasma(9)


##### Load the datasets #####
CMVpos2.data <- Read10X(data.dir = paste0(path, "data_CMVpos2"), strip.suffix = T)
CMVpos3.data <- Read10X(data.dir = paste0(path, "data_CMVpos3"), strip.suffix = T)
CMVpos4.data <- Read10X(data.dir = paste0(path, "data_CMVpos4"), strip.suffix = T)
CMVneg1.data <- Read10X(data.dir = paste0(path, "data_CMVneg1"), strip.suffix = T)
CMVneg2.data <- Read10X(data.dir = paste0(path, "data_CMVneg2"), strip.suffix = T)
LibA.data <- Read10X(data.dir = paste0(path, "data_CMVpos1_5/LibA/filtered_feature_bc_matrix"))
LibB.data <- Read10X(data.dir = paste0(path, "data_CMVpos1_5/LibB/filtered_feature_bc_matrix"), strip.suffix = T))

# Add suffix '-2' to cellnames of LibB to allow merging
colnames(LibB.data) <- paste0(colnames(LibB.data), "-2")

# Initialize the Seurat objects with the raw (non-normalized data).
# Keep all features expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected features
CMVpos2 <- CreateSeuratObject(counts = CMVpos2.data, min.cells = 3, min.features = 200, project = "CMVpos2")
CMVpos3 <- CreateSeuratObject(counts = CMVpos3.data, min.cells = 3, min.features = 200, project = "CMVpos3")
CMVpos4 <- CreateSeuratObject(counts = CMVpos4.data, min.cells = 3, min.features = 200, project = "CMVpos4")
CMVneg1 <- CreateSeuratObject(counts = CMVneg1.data, min.cells = 3, min.features = 200, project = "CMVneg1")
CMVneg2 <- CreateSeuratObject(counts = CMVneg2.data, min.cells = 3, min.features = 200, project = "CMVneg2")
LibA <- CreateSeuratObject(counts = LibA.data, min.cells = 3, min.features = 200, project = "CMVpos1_5")
LibB <- CreateSeuratObject(counts = LibB.data, min.cells = 3, min.features = 200, project = "CMVpos1_5")

rm(CMVpos2.data, CMVpos3.data, CMVpos4.data, CMVneg1.data, CMVneg2.data, LibA.data, LibB.data)
CMVpos2
CMVpos3
CMVpos4
CMVneg1
CMVneg2
LibA
LibB

##### QC #####
# store mitochondrial percentage in object meta data
CMVpos2 <- PercentageFeatureSet(CMVpos2, pattern = "^MT-", col.name = "percent.mt")
CMVpos3 <- PercentageFeatureSet(CMVpos3, pattern = "^MT-", col.name = "percent.mt")
CMVpos4 <- PercentageFeatureSet(CMVpos4, pattern = "^MT-", col.name = "percent.mt")
CMVneg1 <- PercentageFeatureSet(CMVneg1, pattern = "^MT-", col.name = "percent.mt")
CMVneg2 <- PercentageFeatureSet(CMVneg2, pattern = "^MT-", col.name = "percent.mt")
CMVneg2 <- PercentageFeatureSet(CMVneg2, pattern = "^MT-", col.name = "percent.mt")
LibA <- PercentageFeatureSet(LibA, pattern = "^MT-", col.name = "percent.mt")
LibB <- PercentageFeatureSet(LibB, pattern = "^MT-", col.name = "percent.mt")

# Inspect Metadata
VlnPlot(CMVpos2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVpos3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVpos4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVneg1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVneg2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LibA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LibB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(CMVpos2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVpos2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVpos3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVpos3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVpos4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVpos4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVneg1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVneg1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVneg2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVneg2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(LibA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LibA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(LibB, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LibB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


# Filter based on QC metrics
CMVpos2 <- subset(CMVpos2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 6)
CMVpos3 <- subset(CMVpos3, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVpos4 <- subset(CMVpos4, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg1 <- subset(CMVneg1, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg2 <- subset(CMVneg2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
LibA <- subset(LibA, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 10)
LibB <- subset(LibB, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 10)


# Put Seurat objects into a list
my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, LibA, LibB)
#my_SeuratList <- list(LibA, LibB)
rm(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, LibA, LibB)

# Visualize markers of contaminating cells
my_SeuratList %>% 
  map(function(x) VlnPlot(object = x, features = c("HBA1", "HBA2", "IGJ"), slot = "counts"))

# Filter erythrocytes and b cells out
my_SeuratList %>% 
  map(function(x) subset(x, subset= HBA1 == 0 & HBA2 == 0 & IGJ == 0)) -> my_SeuratList
my_SeuratList

##### Scran normalization #####
scran_normalize <- function(seu){
  sce <- as.SingleCellExperiment(seu)
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = clusters)
  summary(sizeFactors(sce))
  sce <- logNormCounts(sce, log = FALSE)
  
  # Add normalized values and sizefactors to Seurat container
  lognormcounts <- log(sce@assays@data$normcounts+1)
  lognormcounts <- as(lognormcounts, "dgCMatrix")
  dimnames(lognormcounts) <- dimnames(seu[["RNA"]]@counts)
  
  sizefactors <- sizeFactors(sce)
  names(sizefactors) <- names(seu$nCount_RNA)
  
  seu <- AddMetaData(seu, sizefactors, col.name = "sizefactors")
  seu[["RNA"]]@data <- lognormcounts
  rm(sce)
  return(seu)
}

mclapply(my_SeuratList, mc.cores = ncores, function(x) scran_normalize(x)) -> my_SeuratList

# Visualize sizefactors
FeatureScatter(my_SeuratList[[2]], feature1 = "sizefactors", feature2 = "nCount_RNA")
FeatureScatter(my_SeuratList[[2]], feature1 = "sizefactors", feature2 = "nFeature_RNA")



# Add CITE-Seq data for LibA, LibB separately,for those it was counted with kite, so take these out for now
LibA <- my_SeuratList[[6]]
LibB <- my_SeuratList[[7]]

CMVpos1_5 <- merge(LibA, LibB, merge.data = TRUE, project = "CMVpos1_5")

my_SeuratList <- c(my_SeuratList[[1]], my_SeuratList[[2]], my_SeuratList[[3]],
                   my_SeuratList[[4]], my_SeuratList[[5]])


##### Add normalized CITE-Seq data and exclude doublets #####
# Define marker/hashtag panel and assign populations in the same order

my_markers <- c("CD62L-adt", "CD57-adt", "CD161-adt", "Anti_Biotin_NKG2A-adt", "CD56-adt",
                "CD16-adt", "CD2-adt", "KLRG1-adt", "CXCR4-adt", "KIR2DL2_L3-adt", "KIR2DL1_S135-adt",
                "KIR3DL1-adt", "LAG3-adt", "NKp30-adt")
my_hashtags <- c("Hashtag1-adt", "Hashtag2-adt")
my_populations <- c("NKG2Cpos", "NKG2Cneg") 

Add_normCITEseq <- function(my_Seurat, name = character()){
  
  # Read in ADT/HTO data
  my_ADT <- as.sparse((x = read.delim(file = paste0(paste0(path, "data_"), name,
                                                    "/CITE-Seq/ATCACGAT.umi.txt"),
                                      header = TRUE, row.names = 1)))
  my_ADT <- my_ADT[rownames(my_ADT)%in% my_markers,] 
  my_HTO <- as.sparse((x = read.delim(file = paste0(paste0(path, "data_"), name,
                                                    "/CITE-Seq/ATTACTCG.umi.txt"),
                                      header = TRUE, row.names = 1)))
  my_HTO <- my_HTO[rownames(my_HTO)%in% my_hashtags,] 
  my_HTO
  
  # Define HTOs
  rownames(my_HTO) <- my_populations
  
  # Trim data to cells present in transcriptomes and ADT/HTO
  cellfilter.adt <- colnames(my_ADT) %in% colnames(my_Seurat)
  cellfilter.transcriptomes <- colnames(my_Seurat)[colnames(my_Seurat) %in% colnames(my_ADT)]
  my_ADT <- my_ADT[,cellfilter.adt]
  my_Seurat <- subset(my_Seurat, cells = cellfilter.transcriptomes)
  cellfilter.hto <- colnames(my_HTO) %in% colnames(my_Seurat)
  my_HTO <- my_HTO[,cellfilter.hto]
  
  # Add ADT/HTO to Seurat
  my_Seurat[["ADT"]] <- CreateAssayObject(counts = my_ADT)
  my_Seurat[["HTO"]] <- CreateAssayObject(counts = my_HTO)
  
  # Normalize
  my_Seurat <- NormalizeData(object = my_Seurat, assay = "ADT", normalization.method = "CLR")
  my_Seurat <- NormalizeData(object = my_Seurat, assay = "HTO", normalization.method = "CLR")
  my_Seurat <- ScaleData(object = my_Seurat, assay = "ADT")
  
  # Remove doublets
  my_Seurat <- HTODemux(my_Seurat, assay = "HTO", positive.quantile = 0.99)
  print(table(my_Seurat$HTO_classification.global))
  my_Seurat <- my_Seurat[,my_Seurat$HTO_classification.global=="Singlet"]
  
  return(my_Seurat)
} 

my_SeuratList %>% 
  map2(my_samples[1:5], function(x, y) Add_normCITEseq(x, name = y)) -> my_SeuratList

##### Addition of CITE-Seq data quantified with bustools to LibA/LibB#####
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

LibA_HTO <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibA_HTO/featurecounts/"))
LibB_HTO <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibB_HTO/featurecounts/"))

LibA_ADT <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibA_ADT/featurecounts/"))
LibB_ADT <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibB_ADT/featurecounts/"))




# Harmonize cellnames
colnames(LibA_HTO) <- paste0(colnames(LibA_HTO), "-1")
table(colnames(LibA_HTO)%in%colnames(LibA))
LibA_HTO <- LibA_HTO[,colnames(LibA_HTO)%in%colnames(LibA)]

colnames(LibA_ADT) <- paste0(colnames(LibA_ADT), "-1")
table(colnames(LibA_ADT)%in%colnames(LibA))
LibA_ADT <- LibA_ADT[,colnames(LibA_ADT)%in%colnames(LibA)]

colnames(LibB_HTO) <- paste0(colnames(LibB_HTO), "-2")
table(colnames(LibB_HTO)%in%colnames(LibB))
LibB_HTO <- LibB_HTO[,colnames(LibB_HTO)%in%colnames(LibB)]

colnames(LibB_ADT) <- paste0(colnames(LibB_ADT), "-2")
table(colnames(LibB_ADT)%in%colnames(LibB))
LibB_ADT <- LibB_ADT[,colnames(LibB_ADT)%in%colnames(LibB)]

# Harmonize rownames with other objects
rownames(LibA_ADT) <- paste0(rownames(LibA_ADT), "-adt")
rownames(LibA_ADT)[rownames(LibA_ADT) == "anti-Biotin-adt"] <- "Anti-Biotin-NKG2A-adt"
rownames(LibA_ADT)[rownames(LibA_ADT) == "KIR2DL1-S1-S3-S5-adt"] <- "KIR2DL1-S135-adt"
rownames(LibA_ADT)[rownames(LibA_ADT) == "LAG-3-adt"] <- "LAG3-adt"

rownames(LibB_ADT) <- paste0(rownames(LibB_ADT), "-adt")
rownames(LibB_ADT)[rownames(LibB_ADT) == "anti-Biotin-adt"] <- "Anti-Biotin-NKG2A-adt"
rownames(LibB_ADT)[rownames(LibB_ADT) == "KIR2DL1-S1-S3-S5-adt"] <- "KIR2DL1-S135-adt"
rownames(LibB_ADT)[rownames(LibB_ADT) == "LAG-3-adt"] <- "LAG3-adt"

# One hashtag is mislabelled (CMVneg1 instead of CMVpos1), correct
rownames(LibA_HTO)[4] <- "CMVpos1_NKG2Cneg"
rownames(LibB_HTO)[4] <- "CMVpos1_NKG2Cneg"


# Create Assay Objects
LibA_HTO_assay <- CreateAssayObject(counts = LibA_HTO, min.features = 1)
LibB_HTO_assay <- CreateAssayObject(counts = LibB_HTO, min.features = 1)




# Merge
CMVpos1_5_HTO <- CreateAssayObject(counts= cbind(LibA_HTO_assay@counts, LibB_HTO_assay@counts))

# Include only cells which are present in both modalities
CMVpos1_5 <- subset(CMVpos1_5, cells = Cells(CMVpos1_5_HTO))

# Add HTO to objects
CMVpos1_5[["HTO"]] <- CMVpos1_5_HTO



# There is a problem with hashtag 1 => Allmost all cells are marked by it and the higher counts
# in the actual samples are not sufficient to get precise demultiplexing.
# Luckily, all others work nicely, so we can correct for this by selecting all cells
# which are negative for all other hashtags (thresholds set based on violin plots) and
# setting the value for hashtag 1 to 100, all others to 0. This way we can still use HTODemux
VlnPlot(CMVpos1_5, features = rownames(CMVpos1_5[["HTO"]][1:2]), log = TRUE, ncol =2)
VlnPlot(CMVpos1_5, features = rownames(CMVpos1_5[["HTO"]][3:4]), log = TRUE, ncol =2)
table(LibB_HTO[2,]<10&LibB_HTO[3,]<30&LibB_HTO[4,]<30)
table(LibA_HTO[2,]<10&LibA_HTO[3,]<30&LibA_HTO[4,]<30)

# Add a normal distribution around 8 to be able to visualize
LibA_HTO_corrected <- LibA_HTO
LibB_HTO_corrected <- LibB_HTO

LibA_HTO_corrected[1,] <- round(rnorm(ncol(LibA_HTO), mean = 5, sd = 1))
LibB_HTO_corrected[1,] <- round(rnorm(ncol(LibB_HTO), mean = 5, sd = 1))

LibA_HTO_corrected[1,LibA_HTO[2,]<10&LibA_HTO[3,]<30&LibA_HTO[4,]<30] <- round(rnorm(ncol(LibA_HTO_corrected[,LibA_HTO[2,]<10&LibA_HTO[3,]<30&LibA_HTO[4,]<30]), mean = 100, sd = 15))
LibB_HTO_corrected[1,LibB_HTO[2,]<10&LibB_HTO[3,]<30&LibB_HTO[4,]<30] <- round(rnorm(ncol(LibB_HTO_corrected[,LibB_HTO[2,]<10&LibB_HTO[3,]<30&LibB_HTO[4,]<30]), mean = 100, sd = 15))

# Create assay
# Create Assay Objects
LibA_HTO_assay_corrected <- CreateAssayObject(counts = LibA_HTO_corrected, min.features = 1)
LibB_HTO_assay_corrected <- CreateAssayObject(counts = LibB_HTO_corrected, min.features = 1)

# Merge 
CMVpos1_5_HTO_corrected <- CreateAssayObject(counts= cbind(LibA_HTO_assay_corrected@counts, LibB_HTO_assay_corrected@counts))

# Add corrected assay to Seurat object
CMVpos1_5[["HTO_corrected"]] <- CMVpos1_5_HTO_corrected



# Normalize
CMVpos1_5 <- NormalizeData(CMVpos1_5, assay = "HTO_corrected", normalization.method = "CLR")

# Demultiplex
CMVpos1_5<- HTODemux(CMVpos1_5, assay = "HTO_corrected", positive.quantile = 0.995,)
table(CMVpos1_5$HTO_corrected_classification.global)
CMVpos1_5 <- subset(CMVpos1_5, HTO_corrected_classification.global == "Singlet")

table(CMVpos1_5$HTO_corrected_maxID)

# Check how the distribution looks
RidgePlot(CMVpos1_5, assay = "HTO_corrected", 
          features = rownames(CMVpos1_5[["HTO_corrected"]]), ncol = 2, log = T)



# Add ADT data
LibA_ADT <- CreateAssayObject(counts = LibA_ADT)
LibB_ADT<- CreateAssayObject(counts = LibB_ADT)


# Normalize per Library
LibA_ADT <- NormalizeData(LibA_ADT, normalization.method = "CLR")
LibB_ADT <- NormalizeData(LibB_ADT, normalization.method = "CLR")




# Merge, scale and add to Objects
CMVpos1_5_ADT <- CreateAssayObject(data=cbind(LibA_ADT@data, LibB_ADT@data))
CMVpos1_5_ADT <- subset(CMVpos1_5_ADT, cells = Cells(CMVpos1_5))


CMVpos1_5[["ADT"]] <- CMVpos1_5_ADT

# Demultiplex CMVpos1/5 from LibA/LibB
CMVpos1_v3 <- subset(CMVpos1_5, `HTO_corrected_maxID`%in%c("CMVpos1-NKG2Cpos", "CMVpos1-NKG2Cneg"))
CMVpos5_v3 <- subset(CMVpos1_5, `HTO_corrected_maxID`%in%c("CMVpos5-NKG2Cpos", "CMVpos5-NKG2Cneg"))


my_SeuratList <- c(my_SeuratList, CMVpos1_v3, CMVpos5_v3)

# Export QCed datasets
my_SeuratList %>% 
  map2(my_samples, function(x, y) saveRDS(x, file = paste0(y, "_QC_ADT_HTO")))

# Import QCed datasets
CMVpos2 <- readRDS(file = "CMVpos2_QC_ADT_HTO")
CMVpos3 <- readRDS(file = "CMVpos3_QC_ADT_HTO")
CMVpos4 <- readRDS(file = "CMVpos4_QC_ADT_HTO")
CMVneg1 <- readRDS(file = "CMVneg1_QC_ADT_HTO")
CMVneg2 <- readRDS(file = "CMVneg2_QC_ADT_HTO")
CMVpos1_v3 <- readRDS(file = "CMVpos1_v3_QC_ADT_HTO")
CMVpos5_v3 <- readRDS(file = "CMVpos5_v3_QC_ADT_HTO")

my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVpos1_v3, CMVpos5_v3)
rm(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVpos1_v3, CMVpos5_v3)


##### Variable Features, Scaling #####
# Find variable features
my_SeuratList %>% 
  lapply(function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)) -> my_SeuratList




# Scale Data
lapply(my_SeuratList, function(x) ScaleData(x, features = rownames(x))) -> my_SeuratList

# Run PCA
lapply(my_SeuratList, function(x) RunPCA(x, verbose = FALSE)) -> my_SeuratList

# Inspect PCs in Elbowplot and DimHeatmaps
pdf(file = "PCA_results.pdf", title = "PCA Results", paper = "a4")
my_SeuratList %>% 
  map(function(x){
    DimHeatmap(x, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
    DimHeatmap(x, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
    DimHeatmap(x, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)
  })
my_SeuratList %>% 
  map(function(x) ElbowPlot(x))
dev.off()

# RunUMAP
my_SeuratList[[1]] <- RunUMAP(my_SeuratList[[1]], dims = 1:13, verbose = FALSE)
my_SeuratList[[2]] <- RunUMAP(my_SeuratList[[2]], dims = 1:14, verbose = FALSE)
my_SeuratList[[3]] <- RunUMAP(my_SeuratList[[3]], dims = 1:13, verbose = FALSE)
my_SeuratList[[4]] <- RunUMAP(my_SeuratList[[4]], dims = 1:10, verbose = FALSE)
my_SeuratList[[5]] <- RunUMAP(my_SeuratList[[5]], dims = 1:10, verbose = FALSE)
my_SeuratList[[6]] <- RunUMAP(my_SeuratList[[6]], dims = 1:16, verbose = FALSE)
my_SeuratList[[7]] <- RunUMAP(my_SeuratList[[7]], dims = 1:12, verbose = FALSE)

##### Clustering #####
# High resolution in some donors (i.e. overclustering)
# to identify contaminating ILCs (are in proximity to bright)
my_SeuratList[[1]] <- FindNeighbors(my_SeuratList[[1]], dims = 1:13, verbose = FALSE)
my_SeuratList[[1]] <- FindClusters(my_SeuratList[[1]], verbose = FALSE, resolution = 2, n.start = 100)

my_SeuratList[[2]] <- FindNeighbors(my_SeuratList[[2]], dims = 1:14, verbose = FALSE)
my_SeuratList[[2]] <- FindClusters(my_SeuratList[[2]], verbose = FALSE, resolution = 1, n.start = 100)

my_SeuratList[[3]] <- FindNeighbors(my_SeuratList[[3]], dims = 1:13, verbose = FALSE)
my_SeuratList[[3]] <- FindClusters(my_SeuratList[[3]], verbose = FALSE, resolution = 0.6, n.start = 100)

my_SeuratList[[4]] <- FindNeighbors(my_SeuratList[[4]], dims = 1:10, verbose = FALSE)
my_SeuratList[[4]] <- FindClusters(my_SeuratList[[4]], verbose = FALSE, resolution = 0.6, n.start = 100)

my_SeuratList[[5]] <- FindNeighbors(my_SeuratList[[5]], dims = 1:10, verbose = FALSE)
my_SeuratList[[5]] <- FindClusters(my_SeuratList[[5]], verbose = FALSE, resolution = 0.4, n.start = 100)

my_SeuratList[[6]] <- FindNeighbors(my_SeuratList[[6]], dims = 1:16, verbose = FALSE)
my_SeuratList[[6]] <- FindClusters(my_SeuratList[[6]], verbose = FALSE, resolution = 0.6, n.start = 100)

my_SeuratList[[7]] <- FindNeighbors(my_SeuratList[[7]], dims = 1:12, verbose = FALSE)
my_SeuratList[[7]] <- FindClusters(my_SeuratList[[7]], verbose = FALSE, resolution = 0.6, n.start = 100)

# Plot UMAPs with clustering and expression of ILC/NK markers
pdf(file = "UMAP_all_donors_ILCs.pdf", paper = "a4")
my_SeuratList %>% 
  map(function(x) DimPlot(x, label = T))
my_SeuratList %>% 
  map(function(x) FeaturePlot(x, features = c("CD56-adt", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
                              min.cutoff = "q1", max.cutoff = "q99", order = TRUE))
dev.off()



##### Removal of contaminating cells (mostly ILC2s: High IL7R, GATA3, IL2RA, CD40LG, negative for CD56) #####
my_SeuratList[[1]] <- subset(my_SeuratList[[1]], idents = "16", invert = TRUE) # Cluster 16 are ILC2
my_SeuratList[[2]] <- subset(my_SeuratList[[2]], idents = "9", invert = TRUE) # Cluster 9 are ILC2
my_SeuratList[[3]] <- subset(my_SeuratList[[3]], idents = "6", invert = TRUE) # Cluster 6 are ILC2
my_SeuratList[[4]] <- subset(my_SeuratList[[4]], idents = "7", invert = TRUE) # Cluster 7 are ILC2
my_SeuratList[[5]] <- subset(my_SeuratList[[5]], idents = "6", invert = TRUE) # Cluster 6 are ILC2
my_SeuratList[[6]] <- subset(my_SeuratList[[6]], idents = "7", invert = TRUE) # Cluster 7 are ILC2
my_SeuratList[[7]] <- subset(my_SeuratList[[7]], idents = "8", invert = TRUE) # Cluster 8 are ILC2


##### Re-run Dimred and Clustering #####
# Run PCA
my_SeuratList %>% 
  lapply(function(x) RunPCA(x, verbose = FALSE)) -> my_SeuratList

# Inspect PCs in Elbowplot and DimHeatmaps
pdf(file = "PCA_results_no_ILCs.pdf", title = "PCA Results", paper = "a4")
my_SeuratList %>% 
  map(function(x){
    DimHeatmap(x, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
    DimHeatmap(x, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
  })
my_SeuratList %>% 
  map(function(x) ElbowPlot(x))
dev.off()

# RunUMAP
my_SeuratList[[1]] <- RunUMAP(my_SeuratList[[1]], dims = 1:11, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[2]] <- RunUMAP(my_SeuratList[[2]], dims = 1:14, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[3]] <- RunUMAP(my_SeuratList[[3]], dims = 1:15, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[4]] <- RunUMAP(my_SeuratList[[4]], dims = 1:10, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[5]] <- RunUMAP(my_SeuratList[[5]], dims = 1:13, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[6]] <- RunUMAP(my_SeuratList[[6]], dims = 1:14, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[7]] <- RunUMAP(my_SeuratList[[7]], dims = 1:14, verbose = FALSE, n.neighbors = 20)





# Add donor annotation
CMVpos2$donor <- "CMVpos2"
CMVpos3$donor <- "CMVpos3"
CMVpos4$donor <- "CMVpos4"
CMVneg1$donor <- "CMVneg1"
CMVneg2$donor <- "CMVneg2"
CMVpos1_v3$donor <- "CMVpos1"
CMVpos5_v3$donor <- "CMVpos5"


# Put Seurat objects into a list
my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVpos1_v3, CMVpos5_v3)
rm(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVpos1_v3, CMVpos5_v3)

# Assign serostatus
my_SeuratList[[1]]$serostatus <- "CMVpos"
my_SeuratList[[2]]$serostatus <- "CMVpos"
my_SeuratList[[3]]$serostatus <- "CMVpos"
my_SeuratList[[4]]$serostatus <- "CMVneg"
my_SeuratList[[5]]$serostatus <- "CMVneg"
my_SeuratList[[6]]$serostatus <- "CMVpos"
my_SeuratList[[7]]$serostatus <- "CMVpos"


##### Integration of datasets #####
# Find variable features, scale
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 128000 * 1024^2) # for 128 Gb RAM

# Find overlapping variable features of CMVpos
shared_variable_features_CMVpos <- Reduce(intersect, list(VariableFeatures(my_SeuratList[[1]]),
                       VariableFeatures(my_SeuratList[[2]]),
                       VariableFeatures(my_SeuratList[[3]]),
                       VariableFeatures(my_SeuratList[[6]]),
                       VariableFeatures(my_SeuratList[[7]])))

shared_features <- Reduce(intersect, list(rownames(my_SeuratList[[1]]),
                                          rownames(my_SeuratList[[2]]),
                                          rownames(my_SeuratList[[3]]),
                                          rownames(my_SeuratList[[4]]),
                                          rownames(my_SeuratList[[5]]),
                                          rownames(my_SeuratList[[6]]),
                                          rownames(my_SeuratList[[7]])
                                          ))

all_features <- Reduce(unique, list(rownames(my_SeuratList[[1]]),
                                       rownames(my_SeuratList[[2]]),
                                       rownames(my_SeuratList[[3]]),
                                       rownames(my_SeuratList[[4]]),
                                       rownames(my_SeuratList[[5]]),
                                       rownames(my_SeuratList[[6]]),
                                    rownames(my_SeuratList[[7]])
                                    ))


# Integrate CMVpos
integration_anchors_CMVpos <- FindIntegrationAnchors(c(my_SeuratList[[1]], my_SeuratList[[2]], my_SeuratList[[3]], my_SeuratList[[6]], my_SeuratList[[7]]),
                                                     dims = 1:10, anchor.features = shared_variable_features_CMVpos)
CMVpos <- IntegrateData(integration_anchors_CMVpos, dims = 1:10, features.to.integrate = all_features)
DefaultAssay(CMVpos) <- "integrated"

# Process integrated dataset
CMVpos <- ScaleData(CMVpos, vars.to.regress = "nCount_RNA")


CMVpos <- RunPCA(CMVpos)

CMVpos <- RunUMAP(CMVpos, dims = 1:10, n.neighbors = 20, min.dist = 0.2)
CMVpos <- FindNeighbors(CMVpos, dims = 1:10)
CMVpos <- FindClusters(CMVpos, resolution = 0.2)

# Subcluster to resolve early dim
temp <- subset(CMVpos, idents = "2")
temp <- FindNeighbors(temp, dims = 1:10)
temp <- FindClusters(temp, resolution = 0.2)
DimPlot(temp)
Idents(CMVpos, cells = WhichCells(temp, idents = "0")) <- "EarlyCD56dim"
DimPlot(CMVpos)

CMVpos <- RenameIdents(CMVpos, "2" = "CD56bright", "1" = "CD56dim", "0" = "Adaptive", "3" = "Proliferating")
levels(CMVpos) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating",  "Adaptive")
CMVpos$annotation <- Idents(CMVpos)

DimPlot(CMVpos, cols = c(clusterpal[1:4], adaptivepal))&NoLegend()

# Name populations based on hashtags
CMVpos$population <- "NKG2Cneg"
CMVpos$population[CMVpos$HTO_maxID=="NKG2Cpos"] <- "NKG2Cpos"
CMVpos$population[CMVpos$HTO_corrected_maxID%in%c("CMVpos1-NKG2Cpos","CMVpos5-NKG2Cpos")] <- "NKG2Cpos"

DimPlot(CMVpos, group.by = "population", cols = c("grey",colorscale[7]), shuffle = TRUE)&NoLegend()

# Integration of CMVneg (use the same anchor features as for CMVpos to enable merging afterwards)
integration_anchors_CMVneg <- FindIntegrationAnchors(c(my_SeuratList[[4]], my_SeuratList[[5]]),
                                                     dims = 1:10, anchor.features = shared_variable_features_CMVpos)
CMVneg <- IntegrateData(integration_anchors_CMVneg, dims = 1:10, features.to.integrate = all_features)

# Process integrated dataset
DefaultAssay(CMVneg) <- "integrated"
CMVneg <- ScaleData(CMVneg, vars.to.regress = "nCount_RNA")
CMVneg <- RunPCA(CMVneg)
CMVneg <- RunUMAP(CMVneg, dims = 1:10)
DimPlot(CMVneg)

CMVneg <- FindNeighbors(CMVneg, dims = 1:10)
CMVneg <- FindClusters(CMVneg, resolution = 0.2)

# Subcluster to resolve early dim
temp <- subset(CMVneg, idents = "1")
temp <- FindNeighbors(temp, dims = 1:10)
temp <- FindClusters(temp, resolution = 0.2)
DimPlot(temp)
Idents(CMVneg, cells = WhichCells(temp, idents = "1")) <- "EarlyCD56dim"
DimPlot(CMVneg)

CMVneg <- RenameIdents(CMVneg, "1" = "CD56bright", "0" = "CD56dim", "2" = "Proliferating")
levels(CMVneg) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating")
CMVneg$annotation <- Idents(CMVneg)

# Assign populations also here
CMVneg$population <- "NKG2Cneg"
CMVneg$population[CMVneg$HTO_maxID=="NKG2Cpos"] <- "NKG2Cpos"


DimPlot(CMVneg, group.by = "population", cols = c("grey",colorscale[7]), shuffle = TRUE)

# Normalize RNA assay for full dataset
DefaultAssay(CMVpos) <- "RNA"
DefaultAssay(CMVneg) <- "RNA"
CMVpos <- NormalizeData(CMVpos)
CMVneg <- NormalizeData(CMVneg)
CMVpos <- ScaleData(CMVpos)
CMVneg <- ScaleData(CMVneg)

# Export integrated datasets
saveRDS(CMVpos, file = "CMVpos_scRNA")
saveRDS(CMVneg, file = "CMVneg_scRNA")

# Merge the integrated datasets
integrated_scRNA <- merge(CMVpos, y = CMVneg)
DefaultAssay(integrated_scRNA) <- "integrated"

# Scale, regress out depth effects
VariableFeatures(integrated_scRNA) <- intersect(VariableFeatures(CMVpos), rownames(CMVneg))

integrated_scRNA <- ScaleData(integrated_scRNA, features = VariableFeatures(integrated_scRNA), vars.to.regress = "nCount_RNA")

integrated_scRNA <- RunPCA(integrated_scRNA)

integrated_scRNA <- RunUMAP(integrated_scRNA, dims = c(1:10), min.dist = 0.2)

integrated_scRNA <- FindNeighbors(integrated_scRNA, dims = c(1:10))
integrated_scRNA <- FindClusters(integrated_scRNA, resolution = 0.2)
DimPlot(integrated_scRNA, cols = c(clusterpal[1:4], adaptivepal), pt.size = 0.1)

# Subcluster to resolve early CD56dim
temp <- subset(integrated_scRNA, idents = "2")
temp <- FindNeighbors(temp, dims = c(1:10))
temp <- FindClusters(temp, resolution = 0.2)
DimPlot(temp, pt.size = 0.8, label = T)
Idents(integrated_scRNA, cells =  WhichCells(temp, idents = "1")) <- "EarlyCD56dim"

integrated_scRNA <- RenameIdents(integrated_scRNA, "2" = "CD56bright", "0" = "CD56dim", "1" = "Adaptive", "3" = "Proliferating")
levels(integrated_scRNA) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating",  "Adaptive")
integrated_scRNA$annotation <- Idents(integrated_scRNA)

DimPlot(integrated_scRNA, cols = c(clusterpal[1:4], adaptivepal))


DimPlot(integrated_scRNA, group.by = "serostatus", cols = c("#47BCFF","#FF8A47"), pt.size = 0.6, shuffle = T)
DimPlot(integrated_scRNA, group.by = "serostatus", cols = c("#47BCFF","#FF8A47"), pt.size = 0.6, shuffle = T)+NoLegend()


# For downstream analysis, switch back to RNA assay
# and re-normalize to avoid depth effects
DefaultAssay(integrated_scRNA) <- "RNA"
integrated_scRNA <- NormalizeData(integrated_scRNA)
integrated_scRNA <- ScaleData(integrated_scRNA)

# Export processed dataset
saveRDS(integrated_scRNA, file = "integrated_scRNA")






