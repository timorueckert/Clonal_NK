# mtASAP5 pre-processing

path <- "~/"

# Import data
counts <- Read10X_h5(filename = paste0(path, "data/mtASAP5/mtASAP5_aggr/outs/filtered_peak_bc_matrix.h5"))


metadata <- read.csv(
  file = paste0(path, "data/mtASAP5/mtASAP5_aggr/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

seq_info <- readRDS(paste0(path, "seq_info"))

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seq_info,
  fragments = paste0(path, "data/mtASAP5/mtASAP5_aggr/outs/fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200,
  verbose = TRUE
)


mtASAP5 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


mtASAP5

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(mtASAP5) <- annotations

# compute nucleosome signal score per cell
mtASAP5 <- NucleosomeSignal(object = mtASAP5)

# compute TSS enrichment score per cell
mtASAP5 <- TSSEnrichment(object = mtASAP5, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
mtASAP5$pct_reads_in_peaks <- mtASAP5$peak_region_fragments / mtASAP5$passed_filters * 100
mtASAP5$blacklist_ratio <- mtASAP5$blacklist_region_fragments / mtASAP5$peak_region_fragments

mtASAP5$high.tss <- ifelse(mtASAP5$TSS.enrichment > 1.2, 'High', 'Low')
TSSPlot(mtASAP5, group.by = 'high.tss') + NoLegend()

# Check nucleosome banding pattern
mtASAP5$nucleosome_group <- ifelse(mtASAP5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = mtASAP5, group.by = 'nucleosome_group')

# Plot QC metrics
VlnPlot(
  object = mtASAP5,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# Remove outliers
mtASAP5 <- subset(
  x = mtASAP5,
  subset = peak_region_fragments < 20000 &
    pct_reads_in_peaks > 52 &
    blacklist_ratio < 0.00005 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2 & 
    TSS.enrichment < 10
)

# Export
saveRDS(mtASAP5, file = "mtASAP5_qc")
mtASAP5 <- readRDS(paste0(path, "mtASAP5_qc"))

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
HTO_Lib1 <- import_kite_counts(paste0(path, "data/mtASAP5/HTO/LibA/"))
HTO_Lib2 <- import_kite_counts(paste0(path, "data/mtASAP5/HTO/LibB/"))
HTO_Lib3 <- import_kite_counts(paste0(path, "data/mtASAP5/HTO/LibC/"))
HTO_Lib4 <- import_kite_counts(paste0(path, "data/mtASAP5/HTO/LibD/"))

ADT_Lib1 <- import_kite_counts(paste0(path, "data/mtASAP5/ADT/LibA/"))
ADT_Lib2 <- import_kite_counts(paste0(path, "data/mtASAP5/ADT/LibB/"))
ADT_Lib3 <- import_kite_counts(paste0(path, "data/mtASAP5/ADT/LibC/"))
ADT_Lib4 <- import_kite_counts(paste0(path, "data/mtASAP5/ADT/LibD/"))

# Correct cellname suffix and exclude cells which are not present in ATAC assay
colnames(HTO_Lib1) <- paste0(colnames(HTO_Lib1), "-1")
HTO_Lib1 <- HTO_Lib1[,colnames(HTO_Lib1)%in%colnames(mtASAP5)]
colnames(ADT_Lib1) <- paste0(colnames(ADT_Lib1), "-1")
ADT_Lib1 <- ADT_Lib1[,colnames(ADT_Lib1)%in%colnames(mtASAP5)]

colnames(HTO_Lib2) <- paste0(colnames(HTO_Lib2), "-2")
HTO_Lib2 <- HTO_Lib2[,colnames(HTO_Lib2)%in%colnames(mtASAP5)]
colnames(ADT_Lib2) <- paste0(colnames(ADT_Lib2), "-2")
ADT_Lib2 <- ADT_Lib2[,colnames(ADT_Lib2)%in%colnames(mtASAP5)]

colnames(HTO_Lib3) <- paste0(colnames(HTO_Lib3), "-3")
HTO_Lib3 <- HTO_Lib3[,colnames(HTO_Lib3)%in%colnames(mtASAP5)]
colnames(ADT_Lib3) <- paste0(colnames(ADT_Lib3), "-3")
ADT_Lib3 <- ADT_Lib3[,colnames(ADT_Lib3)%in%colnames(mtASAP5)]

colnames(HTO_Lib4) <- paste0(colnames(HTO_Lib4), "-4")
HTO_Lib4 <- HTO_Lib4[,colnames(HTO_Lib4)%in%colnames(mtASAP5)]
colnames(ADT_Lib4) <- paste0(colnames(ADT_Lib4), "-4")
ADT_Lib4 <- ADT_Lib4[,colnames(ADT_Lib4)%in%colnames(mtASAP5)]

# Create Seurat objects for ADT/HTO
Lib1_HTO <- CreateSeuratObject(counts = HTO_Lib1, assay = "HTO")
Lib2_HTO <- CreateSeuratObject(counts = HTO_Lib2, assay = "HTO")
Lib3_HTO <- CreateSeuratObject(counts = HTO_Lib3, assay = "HTO")
Lib4_HTO <- CreateSeuratObject(counts = HTO_Lib4, assay = "HTO")

Lib1_ADT <- CreateAssayObject(counts = ADT_Lib1)
Lib2_ADT <- CreateAssayObject(counts = ADT_Lib2)
Lib3_ADT <- CreateAssayObject(counts = ADT_Lib3)
Lib4_ADT <- CreateAssayObject(counts = ADT_Lib4)


Lib1_HTO <- NormalizeData(Lib1_HTO, assay = "HTO", normalization.method = "CLR")
Lib2_HTO <- NormalizeData(Lib2_HTO, assay = "HTO", normalization.method = "CLR")
Lib3_HTO <- NormalizeData(Lib3_HTO, assay = "HTO", normalization.method = "CLR")
Lib4_HTO <- NormalizeData(Lib4_HTO, assay = "HTO", normalization.method = "CLR")

# Demultiplex based on HTO Counts
Lib1_HTO <- HTODemux(Lib1_HTO, assay = "HTO", positive.quantile = 0.95)
Lib2_HTO <- HTODemux(Lib2_HTO, assay = "HTO", positive.quantile = 0.95)
Lib3_HTO <- HTODemux(Lib3_HTO, assay = "HTO", positive.quantile = 0.95)
Lib4_HTO <- HTODemux(Lib4_HTO, assay = "HTO", positive.quantile = 0.95)

table(Lib1_HTO$HTO_classification.global)
table(Lib2_HTO$HTO_classification.global)
table(Lib3_HTO$HTO_classification.global)
table(Lib4_HTO$HTO_classification.global)

# Inspect separation by HTO
Idents(Lib1_HTO) <- "HTO_maxID"
Idents(Lib2_HTO) <- "HTO_maxID"
Idents(Lib3_HTO) <- "HTO_maxID"
Idents(Lib4_HTO) <- "HTO_maxID"

table(Lib1_HTO$HTO_maxID)
table(Lib2_HTO$HTO_maxID)
table(Lib3_HTO$HTO_maxID)
table(Lib4_HTO$HTO_maxID)

# Signal for Hashtag 12 is really low => A lot of cells would be lost here

# Exclude doublets and negative cells
Lib1_HTO <- subset(Lib1_HTO, HTO_classification.global=="Singlet")
Lib2_HTO <- subset(Lib2_HTO, HTO_classification.global=="Singlet")
Lib3_HTO <- subset(Lib3_HTO, HTO_classification.global=="Singlet")
Lib4_HTO <- subset(Lib4_HTO, HTO_classification.global=="Singlet")


# But separation for the cells that are tagged looks good!
RidgePlot(Lib1_HTO, assay = "HTO",
          features = rownames(Lib1_HTO[["HTO"]]), log = T)
RidgePlot(Lib2_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]), log = T)
RidgePlot(Lib3_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]), log = T)
RidgePlot(Lib4_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]), log = T)

# Demultiplex tagged cells and use for identification of homoplasmic mitochondrial
# mutations for each donor, use this for demultiplexing to not lose
# so many cells




##### Addition of mitochondrial data from mgatk output #####
# load mgatk output
Lib1_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP5/mito/libA/final"))
Lib2_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP5/mito/libB/final"))
Lib3_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP5/mito/libC/final"))
Lib4_mito.data <- ReadMGATK(dir = paste0(path, "data/mtASAP5/mito/libD/final"))



# Correct cellname suffixes
Lib2_mito.data$counts@Dimnames[[2]] <- gsub("1", "2", Lib2_mito.data$counts@Dimnames[[2]])
rownames(Lib2_mito.data$depth) <- gsub("1", "2", rownames(Lib2_mito.data$depth))

Lib3_mito.data$counts@Dimnames[[2]] <- gsub("1", "3", Lib3_mito.data$counts@Dimnames[[2]])
rownames(Lib3_mito.data$depth) <- gsub("1", "3", rownames(Lib3_mito.data$depth))

Lib4_mito.data$counts@Dimnames[[2]] <- gsub("1", "4", Lib4_mito.data$counts@Dimnames[[2]])
rownames(Lib4_mito.data$depth) <- gsub("", "4", rownames(Lib4_mito.data$depth))


# Merge
mtASAP5_mitocounts_aggr <- cbind(Lib1_mito.data$counts, Lib2_mito.data$counts,
                                 Lib3_mito.data$counts, Lib4_mito.data$counts)
mtASAP5_mitodepth_aggr <- rbind(Lib1_mito.data$depth, Lib2_mito.data$depth,
                                Lib3_mito.data$depth, Lib4_mito.data$depth)


# Create Assays, subset to cells present in both objects
mtASAP5_mito <- CreateAssayObject(counts = mtASAP5_mitocounts_aggr)

mtASAP5_mito <- subset(mtASAP5_mito, cells = Cells(mtASAP5)[Cells(mtASAP5)%in%Cells(mtASAP5_mito)])


# add assay and metadata to the seurat object
mtASAP5[["mito"]] <- mtASAP5_mito

mtASAP5 <- AddMetaData(mtASAP5, metadata = (mtASAP5_mitodepth_aggr), col.name = "mtDNA_depth")

# Filter based on mtDNA depth
mtASAP5 <- subset(mtASAP5, subset = mtDNA_depth >5)

# Subset cells with clear hashtag identities 
mtASAP5_demux <- subset(mtASAP5, cells = c(colnames(Lib1_HTO), colnames(Lib2_HTO),
                                           colnames(Lib3_HTO), colnames(Lib4_HTO)))

# Subset HTO objects 
Lib1_HTO <- subset(Lib1_HTO, cells = Cells(mtASAP5_demux))
Lib2_HTO <- subset(Lib2_HTO, cells = Cells(mtASAP5_demux))
Lib3_HTO <- subset(Lib3_HTO, cells = Cells(mtASAP5_demux))
Lib4_HTO <- subset(Lib4_HTO, cells = Cells(mtASAP5_demux))

# Add HTO IDs to ATAC Seurat object
HTO_id <- c(as.character(Lib1_HTO$HTO_maxID), as.character(Lib2_HTO$HTO_maxID),
            as.character(Lib3_HTO$HTO_maxID),as.character(Lib4_HTO$HTO_maxID))

HTO_id <- factor(HTO_id)

mtASAP5_demux$HTO_maxID <- HTO_id

# All donors have similar mt-coverage
VlnPlot(mtASAP5_demux, features = "mtDNA_depth", group.by = "HTO_maxID")

# Subset donors based on hashtags
CMVneg3 <- subset(mtASAP5_demux, subset = HTO_maxID == "Hashtag10")
CMVneg4 <- subset(mtASAP5_demux, subset = HTO_maxID == "Hashtag12")
CMVpos1 <- subset(mtASAP5_demux, subset = HTO_maxID == "Hashtag13")
CMVpos3 <- subset(mtASAP5_demux, subset = HTO_maxID == "Hashtag14")

# Find homoplasmic mutations
mito.refallele <- Lib1_mito.data$refallele
variable.sites_CMVneg3 <- IdentifyVariants(CMVneg3, assay = "mito", refallele = mito.refallele)

variable.sites_CMVneg4 <- IdentifyVariants(CMVneg4, assay = "mito", refallele = mito.refallele)

variable.sites_CMVpos1 <- IdentifyVariants(CMVpos1, assay = "mito", refallele = mito.refallele)

variable.sites_CMVpos3 <- IdentifyVariants(CMVpos3, assay = "mito", refallele = mito.refallele)

# Filter for donor-specific homoplasmic mutations

ubiq_CMVneg3 <- subset(
  variable.sites_CMVneg3, subset = n_cells_over_20 >= 2400 &
    strand_correlation >= 0.80
)

ubiq_CMVneg4 <- subset(
  variable.sites_CMVneg4, subset = n_cells_over_20 >= 850 &
    strand_correlation >= 0.80
)

ubiq_CMVpos1 <- subset(
  variable.sites_CMVpos1, subset = n_cells_over_20 >= 1900 &
    strand_correlation >= 0.80
)

ubiq_CMVpos3 <- subset(
  variable.sites_CMVpos3, subset = n_cells_over_20 >= 1800 &
    strand_correlation >= 0.80
)

homo_CMVneg3 <- ubiq_CMVneg3$variant[!ubiq_CMVneg3$variant%in%c(ubiq_CMVneg4$variant, ubiq_CMVpos1$variant, ubiq_CMVpos3$variant)]
homo_CMVneg4 <- ubiq_CMVneg4$variant[!ubiq_CMVneg4$variant%in%c(ubiq_CMVneg3$variant, ubiq_CMVpos1$variant, ubiq_CMVpos3$variant)]
homo_CMVpos1 <- ubiq_CMVpos1$variant[!ubiq_CMVpos1$variant%in%c(ubiq_CMVneg3$variant, ubiq_CMVneg4$variant, ubiq_CMVpos3$variant)]
homo_CMVpos3 <- ubiq_CMVpos3$variant[!ubiq_CMVpos3$variant%in%c(ubiq_CMVneg3$variant, ubiq_CMVneg4$variant, ubiq_CMVpos1$variant)]

homoplasmic_muts <- c(homo_CMVneg3, homo_CMVneg4, homo_CMVpos1, homo_CMVpos3)

# Export for later QC
homoplasmic_muts_df <- data.frame(variant = homoplasmic_muts,
                                  donor = c(rep("CMVneg3", length(homo_CMVneg3)),
                                            rep("CMVneg4", length(homo_CMVneg4)),
                                            rep("CMVpos1", length(homo_CMVpos1)),
                                            rep("CMVpos3", length(homo_CMVpos3))
                                            )) 
write.csv(homoplasmic_muts_df, file = "homoplasmic_muts_mtASAP5")

# Genotype homoplasmic muts in full dataset
mtASAP5 <- AlleleFreq(
  object = mtASAP5,
  variants = homoplasmic_muts,
  assay = "mito"
)

# Dimred
mtASAP5 <- RunTFIDF(mtASAP5)
mtASAP5 <- FindTopFeatures(mtASAP5, min.cutoff = "5")
mtASAP5 <- RunSVD(mtASAP5)
mtASAP5 <- RunUMAP(mtASAP5, dims = 2:30, reduction = "lsi")

# Genotypes nicely match with ATAC UMAP!
# => can confidently be used for demux
FeaturePlot(mtASAP5, features = homo_CMVneg3[1:3])
FeaturePlot(mtASAP5, features = homo_CMVneg4[1:3])
FeaturePlot(mtASAP5, features = homo_CMVpos1[1:3])
FeaturePlot(mtASAP5, features = homo_CMVpos3[1:3])

# Cluster for demux
DefaultAssay(mtASAP5) <- "alleles"

# Expected number of clusters: 4
mtASAP5 <- FindClonotypes(mtASAP5, k = 200)
mtASAP5 <- FindClusters(mtASAP5, resolution = 0.4, group.singletons = F)
DimPlot(mtASAP5, pt.size = 0.4)
# Cluster 4 are doublets: Spread all over and higher number of fragments
# => Exclude
DimPlot(mtASAP5, cells.highlight = WhichCells(mtASAP5, idents = "4"))
VlnPlot(mtASAP5, features = "peak_region_fragments")
table(Idents(mtASAP5))

# Assign donors
# by donors
mtASAP5$donor <- factor(rep("undefined", ncol(mtASAP5)), levels = c("CMVneg3", "CMVneg4","CMVpos1","CMVpos3", "Doublets", "undefined"))
mtASAP5$donor[mtASAP5$seurat_clusters=="0"] <- "CMVneg3"
mtASAP5$donor[mtASAP5$seurat_clusters=="1"] <- "CMVpos3"
mtASAP5$donor[mtASAP5$seurat_clusters=="2"] <- "CMVpos1"
mtASAP5$donor[mtASAP5$seurat_clusters=="3"] <- "CMVneg4"
mtASAP5$donor[mtASAP5$seurat_clusters=="4"] <- "Doublets"

DimPlot(mtASAP5, group.by = "donor", label = T)&NoLegend()
table(mtASAP5$donor)

# Export
saveRDS(mtASAP5, file = "mtASAP5_qc_demux")

# Remove doublets
mtASAP5 <- subset(mtASAP5, idents = "4", invert = TRUE)

# Normalize ADT data
Lib1_ADT <- NormalizeData(Lib1_ADT, normalization.method = "CLR")
Lib2_ADT <- NormalizeData(Lib2_ADT, normalization.method = "CLR")
Lib3_ADT <- NormalizeData(Lib3_ADT, normalization.method = "CLR")
Lib4_ADT <- NormalizeData(Lib4_ADT, normalization.method = "CLR")


# Scale ADT data
Lib1_ADT <- ScaleData(Lib1_ADT)
Lib2_ADT <- ScaleData(Lib2_ADT)
Lib3_ADT <- ScaleData(Lib3_ADT)
Lib4_ADT <- ScaleData(Lib4_ADT)

# Set key to ADT
Lib1_ADT@key <- "ADT"
Lib2_ADT@key <- "ADT"
Lib3_ADT@key <- "ADT"
Lib4_ADT@key <- "ADT"

ADT_full <- merge(x= Lib1_ADT, y =c(Lib2_ADT, Lib3_ADT, Lib4_ADT))

ADT_full <- subset(ADT_full, cells = Cells(mtASAP5))

# Add ADT assay
mtASAP5[["ADT"]] <- ADT_full

# Signal looks good
FeaturePlot(mtASAP5, features = c( "KIR2DL1-S1-S3-S5", "KIR2DL2-L3",
                                   "KIR3DL1","NKp30",
                                   "Anti-PE", "anti-Biotin"), cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE)&NoLegend()&NoAxes()

# Add experiment information
mtASAP5$experiment <- "mtASAP5"

# Count peaks from integrated analysis
integrated <- readRDS("integrated_processed")
peaks <- granges(integrated)
DefaultAssay(mtASAP5) <- "peaks"

New_counts <- FeatureMatrix(
  fragments = Fragments(mtASAP5),
  features = peaks,
  cells = colnames(mtASAP5)
)

# Add assay to late time point
mtASAP5_assay <- CreateChromatinAssay(New_counts, fragments = Fragments(mtASAP5), genome = seq_info)
mtASAP5[["MACS2"]] <- mtASAP5_assay

# Remove original ATAC assay
DefaultAssay(mtASAP5) <- "MACS2"
mtASAP5[["peaks"]] <- NULL

# Add serostatus
mtASAP5$serostatus <- factor(rep("undefined", ncol(mtASAP5)), levels = c("CMVpos", "CMVneg", "undefined"))
mtASAP5$serostatus[mtASAP5$donor%in%c("CMVneg3", "CMVneg4")] <- "CMVneg"
mtASAP5$serostatus[mtASAP5$donor%in%c("CMVpos1","CMVpos3")] <- "CMVpos"

# Add Gene activities
gene.activities <- GeneActivity(mtASAP5)


# add the gene activity matrix to the Seurat object as a new assay and normalize it
mtASAP5[['RNA_atac']] <- CreateAssayObject(counts = gene.activities)
mtASAP5 <- NormalizeData(
  object = mtASAP5,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(mtASAP5$nCount_RNA_atac)
)


# Export processed object
saveRDS(mtASAP5, file = paste0(path, "mtASAP5_processed"))

##### Pre-processing of individual donors #####

#### CMVpos3 ####
CMVpos3_late <- subset(mtASAP5, subset = donor == "CMVpos3")

#### Check demultiplexing of late CMVpos3
FeaturePlot(CMVpos3_late, features = c(homoplasmic_muts_df$variant[1:3],
                                       homoplasmic_muts_df$variant[7:9],
                                       homoplasmic_muts_df$variant[20:22],
                                       homoplasmic_muts_df$variant[41:43]),
            cols = colorscale, ncol = 3, order = TRUE)&NoLegend()&NoAxes()

# There is a contamination with cells from CMVneg3 => remove
DefaultAssay(CMVpos3_late) <- "alleles"
CMVpos3_late <- FindClonotypes(CMVpos3_late, k = 10)
CMVpos3_late <- FindClusters(CMVpos3_late, resolution = 1.4)
VlnPlot(CMVpos3_late, features = c(homoplasmic_muts_df$variant[1:3],
                                   homoplasmic_muts_df$variant[7:9],
                                   homoplasmic_muts_df$variant[20:22],
                                   homoplasmic_muts_df$variant[41:43]), ncol = 3)
# Cells from CMVneg3 fall into cluster 3 => remove
CMVpos3_late <- subset(CMVpos3_late, idents = "3", invert = TRUE)

##### Dimred
DefaultAssay(CMVpos3_late) <- "MACS2"
CMVpos3_late <- RunTFIDF(CMVpos3_late)
CMVpos3_late <- FindTopFeatures(CMVpos3_late, min.cutoff = 10)
CMVpos3_late <- RunSVD(CMVpos3_late)
CMVpos3_late <- RunUMAP(CMVpos3_late, reduction = "lsi", dims = 2:30)
DimPlot(CMVpos3_late)

saveRDS(CMVpos3_late, file = paste0(path, "CMVpos3_late"))

#### CMVpos1 ####
CMVpos1_late <- subset(mtASAP5, subset = donor == "CMVpos1")

# Check demultiplexing from mtASAP5
FeaturePlot(CMVpos1_late, features = c(homoplasmic_muts_df$variant[1:3],
                                       homoplasmic_muts_df$variant[7:9],
                                       homoplasmic_muts_df$variant[20:22],
                                       homoplasmic_muts_df$variant[41:43]),
            cols = colorscale, ncol = 3, order = TRUE)&NoLegend()&NoAxes()

# There is a contamination with cells from CMVneg4 => remove
DefaultAssay(CMVpos1_late) <- "alleles"
CMVpos1_late <- FindClonotypes(CMVpos1_late, k = 10)
CMVpos1_late <- FindClusters(CMVpos1_late, resolution = 0.8)
VlnPlot(CMVpos1_late, features = c(homoplasmic_muts_df$variant[1:3],
                                   homoplasmic_muts_df$variant[7:9],
                                   homoplasmic_muts_df$variant[20:22],
                                   homoplasmic_muts_df$variant[41:43]), ncol = 3)

# Cells from CMVneg4 fall into cluster 5 => remove
CMVpos1_late <- subset(CMVpos1_late, idents = "5", invert = TRUE)

DefaultAssay(CMVpos1_late) <- "MACS2"

# Dimred
CMVpos1_late <- RunTFIDF(CMVpos1_late)
CMVpos1_late <- FindTopFeatures(CMVpos1_late, min.cutoff = 10)
CMVpos1_late <- RunSVD(CMVpos1_late)
CMVpos1_late <- RunUMAP(CMVpos1_late, reduction = "lsi", dims = 2:30, min.dist = 0.2)
DimPlot(CMVpos1_late)

# Clustering and annotation 
DefaultAssay(CMVpos1_late) <- "ADT"
FeaturePlot(CMVpos1_late, features = c("CD56", "CD62L", "CD57", "anti-Biotin",
                                       "Anti-PE", "KIR2DL1-S1-S3-S5", "KIR2DL2-L3",
                                       "KIR3DL1", "NKp30", "CD2", "CD161", "ILT2"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99",
            order = TRUE)&NoLegend()&NoAxes()


DefaultAssay(CMVpos1_late) <- "MACS2"
CMVpos1_late <- FindNeighbors(CMVpos1_late, reduction = "lsi",
                              dims = 2:30)
CMVpos1_late <- FindClusters(CMVpos1_late, resolution = 0.6)
DimPlot(CMVpos1_late, label = T)
# Cluster 3 seems less well defined, has mixed phenotype
# and lower number of fragments
VlnPlot(CMVpos1_late, features = "nCount_MACS2")
#=> Likely low quality cells => remove
CMVpos1_late <- subset(CMVpos1_late, idents = "3", invert = TRUE)

# Re-run UMAP/clustering
CMVpos1_late <- RunUMAP(CMVpos1_late, reduction = "lsi", dims = 2:30)
DimPlot(CMVpos1_late)

CMVpos1_late <- FindNeighbors(CMVpos1_late, reduction = "lsi",
                              dims = 2:30)
CMVpos1_late <- FindClusters(CMVpos1_late, resolution = 0.6)

# Resolve CD56bright
temp <- subset(CMVpos1_late, idents = 0)
temp <- FindNeighbors(temp, reduction = "lsi",
                      dims = 2:30)
temp <- FindClusters(temp, resolution = 0.6)
DimPlot(temp)
Idents(CMVpos1_late, cells = WhichCells(temp, idents = 2)) <- "CD56bright"

CMVpos1_late <- RenameIdents(CMVpos1_late, "0" = "CD56dim",
                             "1" = "Adaptive1", "2" = "Adaptive2", "3" = "Adaptive3",
                             "4" = "Adaptive4", "5" = "Adaptive5")
levels(CMVpos1_late) <- c("CD56bright", "CD56dim",
                          "Adaptive1", "Adaptive2", "Adaptive3", 
                          "Adaptive4", "Adaptive5")
DimPlot(CMVpos1_late, cols = c(clusterpal[c(1,3)], adaptivepal))

# Export
saveRDS(CMVpos1_late, file = paste0(path, "CMVpos1_late"))

#### CMVneg3 ####
CMVneg3 <- subset(mtASAP5, subset = donor == "CMVneg3")
CMVneg3 <- RunTFIDF(CMVneg3)
CMVneg3 <- FindTopFeatures(CMVneg3, min.cutoff = 10)
CMVneg3 <- RunSVD(CMVneg3)
CMVneg3 <- RunUMAP(CMVneg3, reduction = "lsi", dims = 2:30, n.neighbors = 20)

# Check if demultiplexing was complete using homoplasmic mutations
FeaturePlot(CMVneg3, features = c(homoplasmic_muts_df$variant[1:3],
                                  homoplasmic_muts_df$variant[7:9],
                                  homoplasmic_muts_df$variant[20:22],
                                  homoplasmic_muts_df$variant[41:43]),
            cols = colorscale, ncol = 3, order = TRUE)&NoLegend()&NoAxes()
# => Seems complete :)

# Assess phenotype
DefaultAssay(CMVneg3) <- "ADT"
FeaturePlot(CMVneg3, features = c("CD56", "CD62L", "CD57", "Anti-PE",
                                  "NKp30", "KIR2DL1-S1-S3-S5", "KIR2DL2-L3",
                                  "KIR3DL1", "CD328", "ILT2"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE)&NoLegend()&NoAxes()

# Clustering
DefaultAssay(CMVneg3) <- "MACS2"
CMVneg3 <- FindNeighbors(CMVneg3, reduction = "lsi", dims = 2:30)
CMVneg3 <- FindClusters(CMVneg3, resolution = 0.3)
DimPlot(CMVneg3)
CMVneg3$annotation <- Idents(CMVneg3)

saveRDS(CMVneg3, file = paste0(path, "CMVneg3_processed"))



#### CMVneg4 ####
CMVneg4 <- subset(mtASAP5, subset = donor == "CMVneg4")
CMVneg4 <- RunTFIDF(CMVneg4)
CMVneg4 <- FindTopFeatures(CMVneg4, min.cutoff = 10)
CMVneg4 <- RunSVD(CMVneg4)
CMVneg4 <- RunUMAP(CMVneg4, reduction = "lsi", dims = 2:30, n.neighbors = 20)



# Check if demultiplexing was complete using homoplasmic mutations
FeaturePlot(CMVneg4, features = c(homoplasmic_muts_df$variant[1:3],
                                  homoplasmic_muts_df$variant[7:9],
                                  homoplasmic_muts_df$variant[20:22],
                                  homoplasmic_muts_df$variant[41:43]),
            cols = colorscale, ncol = 3, order = TRUE)&NoLegend()&NoAxes()
# => There is a small contamination of cells from CMVneg3 and CMVpos1 => remove
DefaultAssay(CMVneg4) <- "alleles"
CMVneg4 <- FindClonotypes(CMVneg4, k = 10)
CMVneg4 <- FindClusters(CMVneg4, resolution = 0.2)
VlnPlot(CMVneg4, features = c(homoplasmic_muts_df$variant[1:3],
                              homoplasmic_muts_df$variant[7:9],
                              homoplasmic_muts_df$variant[20:22],
                              homoplasmic_muts_df$variant[41:43]), ncol = 3)
# CMVneg3 are in cluster 2
CMVneg4 <- subset(CMVneg4, idents = "2", invert = T)

# Restrict analysis to mutations from CMVpos1/neg4 to catch contaminating cells
# from CMVpos1
homo_muts <- homoplasmic_muts_df %>% dplyr::filter(donor%in%c("CMVneg4", "CMVpos1")) %>% pull(variant)

CMVneg4 <- AlleleFreq(
  object = CMVneg4,
  variants = homo_muts,
  assay = "mito"
)

DefaultAssay(CMVneg4) <- "alleles"
CMVneg4 <- FindClonotypes(CMVneg4, k = 10)
CMVneg4 <- FindClusters(CMVneg4, resolution = 0.6, group.singletons = F)
VlnPlot(CMVneg4, features = c(homoplasmic_muts_df$variant[7:13],
                              homoplasmic_muts_df$variant[20:26]), ncol = 7)
# CMVpos1 is in cluster 6
CMVneg4 <- subset(CMVneg4, idents = "6", invert = T)

# Assess phenotype
DefaultAssay(CMVneg4) <- "ADT"
FeaturePlot(CMVneg4, features = c("CD56", "CD62L", "CD57", "Anti-PE",
                                  "NKp30", "KIR2DL1-S1-S3-S5", "KIR2DL2-L3",
                                  "KIR3DL1", "CD328", "ILT2"),
            cols = colorscale, min.cutoff = "q1", max.cutoff = "q99", order = TRUE)&NoLegend()&NoAxes()

# Clustering
DefaultAssay(CMVneg4) <- "MACS2"
CMVneg4 <- FindNeighbors(CMVneg4, reduction = "lsi", dims = 2:30)
CMVneg4 <- FindClusters(CMVneg4, resolution = 0.3)
DimPlot(CMVneg4)
CMVneg4$annotation <- Idents(CMVneg4)

saveRDS(CMVneg4, file = paste0(path, "CMVneg4_processed"))

