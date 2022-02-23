# Upstream analysis of mtASAP3 experiment

path <- "~/"
# Import data
counts <- Read10X_h5(filename = paste0(path, "data/mtASAP3/mtASAP3_aggr/outs/filtered_peak_bc_matrix.h5"))


metadata <- read.csv(
  file = paste0(path, "data/mtASAP3/mtASAP3_aggr/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

seq_info <- readRDS(paste0(path, "seq_info"))

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seq_info,
  fragments = paste0(path, "data/mtASAP3/mtASAP3_aggr/outs/fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200,
  verbose = TRUE
)


mtASAP3 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


mtASAP3

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(mtASAP3) <- annotations

# compute nucleosome signal score per cell
mtASAP3 <- NucleosomeSignal(object = mtASAP3)

# compute TSS enrichment score per cell
mtASAP3 <- TSSEnrichment(object = mtASAP3, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
mtASAP3$pct_reads_in_peaks <- mtASAP3$peak_region_fragments / mtASAP3$passed_filters * 100
mtASAP3$blacklist_ratio <- mtASAP3$blacklist_region_fragments / mtASAP3$peak_region_fragments

mtASAP3$high.tss <- ifelse(mtASAP3$TSS.enrichment > 1.2, 'High', 'Low')
TSSPlot(mtASAP3, group.by = 'high.tss') + NoLegend()

# Check nucleosome banding pattern
mtASAP3$nucleosome_group <- ifelse(mtASAP3$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = mtASAP3, group.by = 'nucleosome_group')

# Plot QC metrics
VlnPlot(
  object = mtASAP3,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# Remove outliers
mtASAP3 <- subset(
  x = mtASAP3,
  subset = peak_region_fragments < 10000 &
    pct_reads_in_peaks > 60 &
    blacklist_ratio < 0.00005 &
    nucleosome_signal < 5 &
    TSS.enrichment > 3
)

# Export
saveRDS(mtASAP3, file = "mtASAP3_qc")
mtASAP3 <- readRDS(paste0(path, "mtASAP3_qc"))

##### Addition of Antibody counts #####
import_kite_counts <- function(path){
  mtx <- fread(paste0(path,"featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(path,"featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- (fread(paste0(path,"featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

# Import ADT/HTO data
HTO_Lib1 <- import_kite_counts(paste0(path, "data/mtASAP3/CITE/HTO_Lib1/featurecounts/"))
HTO_Lib2 <- import_kite_counts(paste0(path, "data/mtASAP3/CITE/HTO_Lib2/featurecounts/"))
ADT_Lib1 <- import_kite_counts(paste0(path, "data/mtASAP3/CITE/ADT_Lib1/featurecounts/"))
ADT_Lib2 <- import_kite_counts(paste0(path, "data/mtASAP3/CITE/ADT_Lib2/featurecounts/"))

# Correct cellname suffix and exclude cells which are not present in ATAC assay
colnames(HTO_Lib1) <- paste0(colnames(HTO_Lib1), "-1")
HTO_Lib1 <- HTO_Lib1[,colnames(HTO_Lib1)%in%colnames(mtASAP3)]
colnames(ADT_Lib1) <- paste0(colnames(ADT_Lib1), "-1")
ADT_Lib1 <- ADT_Lib1[,colnames(ADT_Lib1)%in%colnames(mtASAP3)]

colnames(HTO_Lib2) <- paste0(colnames(HTO_Lib2), "-2")
HTO_Lib2 <- HTO_Lib2[,colnames(HTO_Lib2)%in%colnames(mtASAP3)]
colnames(ADT_Lib2) <- paste0(colnames(ADT_Lib2), "-2")
ADT_Lib2 <- ADT_Lib2[,colnames(ADT_Lib2)%in%colnames(mtASAP3)]



# Create Seurat objects for ADT/HTO
Lib1_HTO <- CreateSeuratObject(counts = HTO_Lib1, assay = "HTO")
Lib2_HTO <- CreateSeuratObject(counts = HTO_Lib2, assay = "HTO")

Lib1_ADT <- CreateAssayObject(counts = ADT_Lib1)
Lib2_ADT <- CreateAssayObject(counts = ADT_Lib2)


Lib1_HTO <- NormalizeData(Lib1_HTO, assay = "HTO", normalization.method = "CLR")
Lib2_HTO <- NormalizeData(Lib2_HTO, assay = "HTO", normalization.method = "CLR")

# Demultiplex based on HTO Counts
Lib1_HTO <- HTODemux(Lib1_HTO, assay = "HTO", positive.quantile = 0.95)
Lib2_HTO <- HTODemux(Lib2_HTO, assay = "HTO", positive.quantile = 0.95)

table(Lib1_HTO$HTO_classification.global)
table(Lib2_HTO$HTO_classification.global)

# Inspect separation by HTO
Idents(Lib1_HTO) <- "HTO_maxID"
Idents(Lib2_HTO) <- "HTO_maxID"

table(Lib1_HTO$HTO_maxID)
table(Lib2_HTO$HTO_maxID)


RidgePlot(Lib1_HTO, assay = "HTO",
          features = rownames(Lib1_HTO[["HTO"]]), log = T)
RidgePlot(Lib2_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]), log = T)


# Exclude doublets and negative cells
Lib1_HTO <- subset(Lib1_HTO, HTO_classification.global=="Singlet")
Lib2_HTO <- subset(Lib2_HTO, HTO_classification.global=="Singlet")


# Reinspect if doublets were properly excluded
RidgePlot(Lib1_HTO, assay = "HTO",
          features = rownames(Lib1_HTO[["HTO"]]))
RidgePlot(Lib2_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]))
FeatureScatter(Lib1_HTO, feature1 = "Hashtag7", feature2 = "Hashtag8")
FeatureScatter(Lib2_HTO, feature1 = "Hashtag7", feature2 = "Hashtag8")

table(Lib1_HTO$HTO_maxID)
table(Lib2_HTO$HTO_maxID)


# Exclude cells which are not present in HTO assay
mtASAP3 <- subset(mtASAP3, cells = c(colnames(Lib1_HTO), colnames(Lib2_HTO)))

# Add HTO IDs to ATAC Seurat objects
HTO_id <- c(as.character(Lib1_HTO$HTO_maxID), as.character(Lib2_HTO$HTO_maxID))

HTO_id <- factor(HTO_id)

mtASAP3$HTO_maxID <- HTO_id

# Normalize ADT data
Lib1_ADT <- NormalizeData(Lib1_ADT, normalization.method = "CLR")
Lib2_ADT <- NormalizeData(Lib2_ADT, normalization.method = "CLR")


# Scale ADT data
Lib1_ADT <- ScaleData(Lib1_ADT)
Lib2_ADT <- ScaleData(Lib2_ADT)

# Exclude cells which are not present in ADT assay
mtASAP3 <- subset(mtASAP3, cells = c(colnames(Lib1_ADT), colnames(Lib2_ADT)))

# Set key to ADT
Lib1_ADT@key <- "ADT"
Lib2_ADT@key <- "ADT"

ADT_full <- merge(x= Lib1_ADT, y =Lib2_ADT)

ADT_full <- subset(ADT_full, cells = Cells(mtASAP3))

# Add ADT assay
mtASAP3[["ADT"]] <- ADT_full

# Normalize
mtASAP3 <- RunTFIDF(mtASAP3)

# Dim-Reduction
mtASAP3 <- FindTopFeatures(mtASAP3, min.cutoff = 'q1')
mtASAP3 <- RunSVD(mtASAP3)

# Check for correlation with sequencing depth
DepthCor(mtASAP3)

mtASAP3 <- RunUMAP(object = mtASAP3, reduction = 'lsi', dims = 2:30)



# Annotate Hashtags with donor
mtASAP3$donor <- factor(rep("undefined", ncol(mtASAP3)), levels = c("CMVpos2", "CMVpos4"))
mtASAP3$donor[mtASAP3$HTO_maxID=="Hashtag7"] <- "CMVpos2"
mtASAP3$donor[mtASAP3$HTO_maxID=="Hashtag8"] <- "CMVpos4"

# Add experiment information
mtASAP3$experiment <- "mtASAP3"


##### Addition of mitochondrial data from mgatk output #####
# load mgatk output
Lib1_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP3/mgatk_Lib1/final"))
Lib2_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP3/mgatk_Lib2/final"))



# Correct cellname suffixes
Lib2_mito.data$counts@Dimnames[[2]] <- gsub("1", "2", Lib2_mito.data$counts@Dimnames[[2]])
rownames(Lib2_mito.data$depth) <- gsub("1", "2", rownames(Lib2_mito.data$depth))


# Merge
mtASAP3_mitocounts_aggr <- cbind(Lib1_mito.data$counts, Lib2_mito.data$counts)
mtASAP3_mitodepth_aggr <- rbind(Lib1_mito.data$depth, Lib2_mito.data$depth)


# Create Assays, subset to cells present in both objects
mtASAP3_mito <- CreateAssayObject(counts = mtASAP3_mitocounts_aggr)

mtASAP3_mito <- subset(mtASAP3_mito, cells = Cells(mtASAP3)[Cells(mtASAP3)%in%Cells(mtASAP3_mito)])


# add assay and metadata to the seurat object
mtASAP3[["mito"]] <- mtASAP3_mito

mtASAP3 <- AddMetaData(mtASAP3, metadata = (mtASAP3_mitodepth_aggr), col.name = "mtDNA_depth")


# Split donors
CMVpos2_late <- subset(mtASAP3, donor =="CMVpos2")
CMVpos4_late <- subset(mtASAP3, donor =="CMVpos4")

# Use peaks from main analyses of earlier time point to enable merging
# Import early time points
CMVpos2_early <- readRDS("CMVpos2_processed")
CMVpos4_early <- readRDS("CMVpos4_processed")
DefaultAssay(CMVpos2_early) <- "MACS2"
DefaultAssay(CMVpos4_early) <- "MACS2"
DefaultAssay(CMVpos2_late) <- "peaks"
DefaultAssay(CMVpos4_late) <- "peaks"

# Count peaks
peaks <- granges(CMVpos2_early)

New_counts <- FeatureMatrix(
  fragments = Fragments(CMVpos2_late),
  features = peaks,
  cells = colnames(CMVpos2_late)
)

# Add assay to late time point
CMVpos2_assay <- CreateChromatinAssay(New_counts, fragments = Fragments(CMVpos2_late), genome = seq_info)
CMVpos2_late[["MACS2"]] <- CMVpos2_assay

# Count peaks
peaks <- granges(CMVpos4_early)

New_counts1 <- FeatureMatrix(
  fragments = Fragments(CMVpos4_late),
  features = peaks,
  cells = colnames(CMVpos4_late)
)

# Add assay to late time point
CMVpos4_assay <- CreateChromatinAssay(New_counts1, fragments = Fragments(CMVpos4_late), genome = seq_info)
CMVpos4_late[["MACS2"]] <- CMVpos4_assay

# Remove original ATAC assay
DefaultAssay(CMVpos2_early) <- "MACS2"
DefaultAssay(CMVpos4_early) <- "MACS2"
DefaultAssay(CMVpos2_late) <- "MACS2"
DefaultAssay(CMVpos4_late) <- "MACS2"
CMVpos2_late[["peaks"]] <- NULL
CMVpos4_late[["peaks"]] <- NULL
CMVpos2_early[["ATAC"]] <- NULL
CMVpos4_early[["ATAC"]] <- NULL



##### Mitochondrial Genotyping to improve demultiplexing of mtASAP3 (CMVpos2/4) #####
#### CMVpos2 ####
# Filter cells with good mitochondrial coverage
VlnPlot(CMVpos2_late, features = "mtDNA_depth")
CMVpos2_late <- subset(CMVpos2_late, mtDNA_depth > 5)

VlnPlot(CMVpos4_late, features = "mtDNA_depth")
CMVpos4_late <- subset(CMVpos4_late, mtDNA_depth > 5)

# Calling of variable sites and selection of variants with high confidence
mito.refallele <- readRDS(paste0(path, "mito.refallele"))

variable.sites_CMVpos2 <- IdentifyVariants(CMVpos2_late, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos2, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 2)


variable.sites_CMVpos4 <- IdentifyVariants(CMVpos4_late, assay = "mito", refallele = mito.refallele)
VariantPlot(variants = variable.sites_CMVpos4, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 2)


# Select Ubiquituous mutations which can serve to examine/improve hashtag-based demux
ubiq_CMVpos2 <- subset(
  variable.sites_CMVpos2, subset = n_cells_over_20 >= 3900 &
    strand_correlation >= 0.65 &
    vmr < 0.01 & vmr > 0.001
)


ubiq_CMVpos4 <- subset(
  variable.sites_CMVpos4, subset = n_cells_over_20 >= 4000 &
    strand_correlation >= 0.65 &
    vmr < 0.01 & vmr > 0.001
)

# Are there any overlapping mutations? No!
intersect(ubiq_CMVpos2$variant, ubiq_CMVpos4$variant)

# Compute variant frequency per cell for ubiquitous features
CMVpos2_late <- AlleleFreq(
  object = CMVpos2_late,
  variants = c(ubiq_CMVpos2$variant, ubiq_CMVpos4$variant),
  assay = "mito"
)

CMVpos4_late <- AlleleFreq(
  object = CMVpos4_late,
  variants = c(ubiq_CMVpos2$variant, ubiq_CMVpos4$variant),
  assay = "mito"
)

CMVpos2_late[["alleles"]]
CMVpos4_late[["alleles"]]

DefaultAssay(CMVpos2_late) <- "alleles"
DefaultAssay(CMVpos4_late) <- "alleles"


# Find Clonotypes based on ubiq. mutations to identify contaminants and exclude them
CMVpos2_late <- FindClonotypes(CMVpos2_late)
VlnPlot(CMVpos2_late, features = c(ubiq_CMVpos4$variant[1:8],ubiq_CMVpos2$variant[1:8]), ncol = 4) # Cluster 20 contains contaminating cells
table(Idents(CMVpos2_late)) # 47 cells
CMVpos2_late <- subset(CMVpos2_late, idents = "20", invert = TRUE)

CMVpos4_late <- FindClonotypes(CMVpos4_late)
VlnPlot(CMVpos4_late, features = c(ubiq_CMVpos2$variant[1:8],ubiq_CMVpos4$variant[1:8]), ncol = 4) # Cluster 12 contains contaminating cells
table(Idents(CMVpos4_late)) # 25 cells
CMVpos4_late <- subset(CMVpos4_late, idents = "12", invert = TRUE)

# Add serostatus
CMVpos2_late$serostatus <- "CMVpos"
CMVpos4_late$serostatus <- "CMVpos"

# Export processed, cleaned-up datasets from late time point
saveRDS(CMVpos2_late, file = "CMVpos2_late_processed")
saveRDS(CMVpos4_late, file = "CMVpos4_late_processed")