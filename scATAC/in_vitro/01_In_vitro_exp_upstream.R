# Upstream analysis of ASAP seq of NK cells from CMV- individuals stimulated 
# with different culture conditions in vitro

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


# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
viriscale <- viridis(n = 9)
clusterpal <- brewer.pal(8, "Dark2")
pairedpal <- brewer.pal(12, "Paired")
NKG2Cpal <- c("#40BAD5","#120136")
clusterpal2 <- brewer.pal(8, "Set1")


# Import data
counts <- Read10X_h5(filename = paste0(path, "data_NK_act/ASAP_act_aggr/outs/filtered_peak_bc_matrix.h5"))


metadata <- read.csv(
  file = paste0(path, "data_NK_act/ASAP_act_aggr/outs/singlecell.csv"),
  header = TRUE,
  row.names = 1
)




seq_info <- readRDS(paste0(path, "seq_info"))

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seq_info,
  fragments = '/scratch/ATAC_adaptive_NK/data_NK_act/ASAP_act_aggr/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 0,
  verbose = TRUE
)


NK_act_full <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


NK_act_full

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(NK_act_full) <- annotations

# compute nucleosome signal score per cell
NK_act_full <- NucleosomeSignal(object = NK_act_full)

# compute TSS enrichment score per cell
NK_act_full <- TSSEnrichment(object = NK_act_full, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
NK_act_full$pct_reads_in_peaks <- NK_act_full$peak_region_fragments / NK_act_full$passed_filters * 100
NK_act_full$blacklist_ratio <- NK_act_full$blacklist_region_fragments / NK_act_full$peak_region_fragments

NK_act_full$high.tss <- ifelse(NK_act_full$TSS.enrichment > 1.2, 'High', 'Low')
TSSPlot(NK_act_full, group.by = 'high.tss') + NoLegend()

# Check nucleosome banding pattern
NK_act_full$nucleosome_group <- ifelse(NK_act_full$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = NK_act_full, group.by = 'nucleosome_group')

# Plot QC metrics
VlnPlot(
  object = NK_act_full,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# Remove outliers
NK_act_full <- subset(
  x = NK_act_full,
  subset = peak_region_fragments < 16000 &
    pct_reads_in_peaks > 60 &
    blacklist_ratio < 0.00005 &
    nucleosome_signal < 1.25 &
    TSS.enrichment > 2
)


##### Addition of CITE-Seq data #####
import_kite_counts <- function(path){
  mtx <- fread(paste0(path,"featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0(path,"featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- (fread(paste0(path,"featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

# Import CITE-Seq data
HTO_Lib1 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/HTO_Lib1/featurecounts/"))
ADT_Lib1 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/ADT_Lib1/featurecounts/"))
HTO_Lib2 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/HTO_Lib2/featurecounts/"))
ADT_Lib2 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/ADT_Lib2/featurecounts/"))
HTO_Lib3 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/HTO_Lib3/featurecounts/"))
ADT_Lib3 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/ADT_Lib3/featurecounts/"))
HTO_Lib4 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/HTO_Lib4/featurecounts/"))
ADT_Lib4 <- import_kite_counts(paste0(path, "data_NK_act/CITE/output/ADT_Lib4/featurecounts/"))

# Correct cellname suffix and exclude cells which are not present in ATAC assay
colnames(HTO_Lib1) <- paste0(colnames(HTO_Lib1), "-1")
HTO_Lib1 <- HTO_Lib1[,colnames(HTO_Lib1)%in%colnames(NK_act_full)]
colnames(ADT_Lib1) <- paste0(colnames(ADT_Lib1), "-1")
ADT_Lib1 <- ADT_Lib1[,colnames(ADT_Lib1)%in%colnames(NK_act_full)]

colnames(HTO_Lib2) <- paste0(colnames(HTO_Lib2), "-2")
HTO_Lib2 <- HTO_Lib2[,colnames(HTO_Lib2)%in%colnames(NK_act_full)]
colnames(ADT_Lib2) <- paste0(colnames(ADT_Lib2), "-2")
ADT_Lib2 <- ADT_Lib2[,colnames(ADT_Lib2)%in%colnames(NK_act_full)]

colnames(HTO_Lib3) <- paste0(colnames(HTO_Lib3), "-3")
HTO_Lib3 <- HTO_Lib3[,colnames(HTO_Lib3)%in%colnames(NK_act_full)]
colnames(ADT_Lib3) <- paste0(colnames(ADT_Lib3), "-3")
ADT_Lib3 <- ADT_Lib3[,colnames(ADT_Lib3)%in%colnames(NK_act_full)]

colnames(HTO_Lib4) <- paste0(colnames(HTO_Lib4), "-4")
HTO_Lib4 <- HTO_Lib4[,colnames(HTO_Lib4)%in%colnames(NK_act_full)]
colnames(ADT_Lib4) <- paste0(colnames(ADT_Lib4), "-4")
ADT_Lib4 <- ADT_Lib4[,colnames(ADT_Lib4)%in%colnames(NK_act_full)]

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
Lib1_HTO <- HTODemux(Lib1_HTO, assay = "HTO", positive.quantile = 0.99)
Lib2_HTO <- HTODemux(Lib2_HTO, assay = "HTO", positive.quantile = 0.99)
Lib3_HTO <- HTODemux(Lib3_HTO, assay = "HTO", positive.quantile = 0.99)
Lib4_HTO <- HTODemux(Lib4_HTO, assay = "HTO", positive.quantile = 0.99)
table(Lib1_HTO$HTO_classification.global)
table(Lib2_HTO$HTO_classification.global)
table(Lib3_HTO$HTO_classification.global)
table(Lib4_HTO$HTO_classification.global)

# Inspect separation by HTO
table(Lib1_HTO$HTO_maxID)
table(Lib2_HTO$HTO_maxID)
table(Lib3_HTO$HTO_maxID)
table(Lib4_HTO$HTO_maxID)


RidgePlot(Lib1_HTO, assay = "HTO",
          features = rownames(Lib1_HTO[["HTO"]]), log = T)
RidgePlot(Lib2_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]), log = T)
RidgePlot(Lib3_HTO, assay = "HTO",
          features = rownames(Lib3_HTO[["HTO"]]), log = T)
RidgePlot(Lib4_HTO, assay = "HTO",
          features = rownames(Lib4_HTO[["HTO"]]), log = T)

# Exclude doublets and negative cells
Lib1_HTO <- subset(Lib1_HTO, HTO_classification.global=="Singlet")
Lib2_HTO <- subset(Lib2_HTO, HTO_classification.global=="Singlet")
Lib3_HTO <- subset(Lib3_HTO, HTO_classification.global=="Singlet")
Lib4_HTO <- subset(Lib4_HTO, HTO_classification.global=="Singlet")




# Reinspect if doublets were properly excluded
RidgePlot(Lib1_HTO, assay = "HTO",
          features = rownames(Lib1_HTO[["HTO"]]))
RidgePlot(Lib2_HTO, assay = "HTO",
          features = rownames(Lib2_HTO[["HTO"]]))
RidgePlot(Lib3_HTO, assay = "HTO",
          features = rownames(Lib3_HTO[["HTO"]]))
RidgePlot(Lib4_HTO, assay = "HTO",
          features = rownames(Lib4_HTO[["HTO"]]))

table(Lib1_HTO$HTO_maxID)
table(Lib2_HTO$HTO_maxID)
table(Lib3_HTO$HTO_maxID)
table(Lib4_HTO$HTO_maxID)

# Exclude cells which are not present in HTO assay
NK_act_full <- subset(NK_act_full, cells = c(colnames(Lib1_HTO), colnames(Lib2_HTO), colnames(Lib3_HTO), colnames(Lib4_HTO)))

# Add HTO IDs to ATAC Seurat objects
HTO_id <- c(as.character(Lib1_HTO$HTO_maxID), as.character(Lib2_HTO$HTO_maxID),
            as.character(Lib3_HTO$HTO_maxID), as.character(Lib4_HTO$HTO_maxID))

HTO_id <- factor(HTO_id)

NK_act_full$HTO_maxID <- HTO_id

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


# Exclude cells which are not present in ADT assay
NK_act_full <- subset(NK_act_full, cells = c(colnames(Lib1_ADT), colnames(Lib2_ADT),
                                             colnames(Lib3_ADT), colnames(Lib4_ADT)))

# Set key to ADT
Lib1_ADT@key <- "ADT"
Lib2_ADT@key <- "ADT"
Lib3_ADT@key <- "ADT"
Lib4_ADT@key <- "ADT"

ADT_full <- merge(x= Lib1_ADT, y =c(Lib2_ADT, Lib3_ADT, Lib4_ADT))

ADT_full <- subset(ADT_full, cells = Cells(NK_act_full))

# Add ADT assay
NK_act_full[["ADT"]] <- ADT_full

# Annotate Hashtags
# by conditions
NK_act_full$condition <- factor(rep("undefined", ncol(NK_act_full)), levels = c("PQS","PQS_12-18", "LFL", "LFL_12-18", "undefined"))
NK_act_full$condition[NK_act_full$HTO_maxID%in%c("Hashtag1", "Hashtag8")] <- "PQS"
NK_act_full$condition[NK_act_full$HTO_maxID%in%c("Hashtag2", "Hashtag9")] <- "PQS_12-18"
NK_act_full$condition[NK_act_full$HTO_maxID%in%c("Hashtag3", "Hashtag10")] <- "LFL"
NK_act_full$condition[NK_act_full$HTO_maxID%in%c("Hashtag4", "Hashtag12")] <- "LFL_12-18"

# by donors
NK_act_full$donor<- factor(rep("undefined", ncol(NK_act_full)), levels = c("CMVneg1", "CMVneg2", "undefined"))
NK_act_full$donor[NK_act_full$HTO_maxID%in%c("Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4")] <- "CMVneg1"
NK_act_full$donor[NK_act_full$HTO_maxID%in%c("Hashtag8", "Hashtag9", "Hashtag10", "Hashtag12")] <- "CMVneg2"

# Call peaks with MACS2
DefaultAssay(NK_act_full) <- 'peaks'
peaks <- CallPeaks(NK_act_full, macs2.path = "/opt/anaconda3/bin/macs2")

# Count peaks
New_counts <- FeatureMatrix(
  fragments = Fragments(NK_act_full),
  features = peaks,
  cells = colnames(NK_act_full)
)

# Generate Assay object and add to NK_act_full
MACS2_peaks <- CreateChromatinAssay(New_counts, fragments = Fragments(NK_act_full), genome = seq_info)

NK_act_full[["MACS2"]] <- MACS2_peaks 
DefaultAssay(NK_act_full) <- "MACS2"
rm(MACS2_peaks)

# Normalize
NK_act_full <- RunTFIDF(NK_act_full)

# Dim-Reduction
NK_act_full <- FindTopFeatures(NK_act_full, min.cutoff = 'q1')
NK_act_full <- RunSVD(NK_act_full)

# Check for correlation with sequencing depth
DepthCor(NK_act_full)

NK_act_full <- RunUMAP(object = NK_act_full, reduction = 'lsi', dims = 2:30)

NK_act_full <- FindNeighbors(object = NK_act_full, reduction = 'lsi', dims = 2:30, k.param = 30)
NK_act_full <- FindClusters(object = NK_act_full, verbose = FALSE, algorithm = 3, resolution = 0.2)




# Inspect UMAP
pdf(file = "Plots/SupplFig4_dimplot.pdf", paper = "a4")
DimPlot(NK_act_full, group.by = "condition", cols = clusterpal, pt.size = 0.3)
DimPlot(NK_act_full, group.by = "donor", cols = clusterpal, pt.size = 0.3)
DimPlot(NK_act_full,  pt.size = 0.3, cols = c("black", "grey", "grey53"))
dev.off()


# Harmony integration of the two donors
library(harmony)

NK_act_full <- RunHarmony(
  object = NK_act_full,
  group.by.vars = 'donor',
  reduction = 'lsi',
  assay.use = 'MACS2',
  project.dim = FALSE
)


NK_act_full <- RunUMAP(object = NK_act_full, reduction = 'harmony', dims = 2:30)

# Inspect UMAP
DimPlot(NK_act_full, group.by = "condition", cols = pairedpal[c(1,2,5,6)], pt.size = 1)
DimPlot(NK_act_full, group.by = "donor", cols = clusterpal, pt.size = 1)
DimPlot(NK_act_full,  pt.size = 1)


#### Gene Activity Matrix: Proximity Model, to get a first understanding of clusters #####
gene.activities <- GeneActivity(NK_act_full, assay = "MACS2")

# add the gene activity matrix to the Seurat object as a new assay and normalize it
NK_act_full[['RNA_atac']] <- CreateAssayObject(counts = gene.activities)
NK_act_full <- NormalizeData(
  object = NK_act_full,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(NK_act_full$nCount_RNA_atac)
)

# Visualize gene activity
DefaultAssay(NK_act_full) <- 'RNA_atac'

# FindMarkers
DEG_activity <- FindAllMarkers(NK_act_full, only.pos = TRUE, test.use = 'LR',
                               latent.vars = 'peak_region_fragments')

top_DEG <- DEG_activtiy %>% group.by(cluster) %>% top_n(wt = , n = 15)

##### TF motif Analysis #####
DefaultAssay(NK_act_full) <- 'MACS2'
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
  subject = granges(NK_act_full),
  out = 'positions',
  genome = 'hg38'
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(NK_act_full),
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
NK_act_full <- SetAssayData(
  object = NK_act_full,
  assay = 'MACS2',
  slot = 'motifs',
  new.data = motif
)
NK_act_full[["MACS2"]]

# Calculate regional stats
NK_act_full <- RegionStats(object = NK_act_full, genome = BSgenome.Hsapiens.UCSC.hg38)

##### Calculate motif acitivies #####
NK_act_full <- RunChromVAR(
  object = NK_act_full,
  genome = BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(NK_act_full) <- 'chromvar'

# Add assay with motif names rather than IDs for more convenient plotting##
chromvar_names <- GetAssayData(NK_act_full, assay = "chromvar", slot = "data")
rownames(chromvar_names) <- c(ConvertMotifID(NK_act_full, id = rownames(NK_act_full[["chromvar"]])[1:633], assay = "MACS2"),
                              names(Stat_ids))
NK_act_full[["chromvar_names"]] <- CreateAssayObject(data = chromvar_names)
DefaultAssay(NK_act_full) <- 'chromvar_names'

# Export seurat object
saveRDS(NK_act_full, file = paste0(path, "NK_act_full"))
NK_act_full <- readRDS(paste0(path, "NK_act_full"))


# Final QC plots
# Add annotation of peaks as promoter/exon/intron/intergenic
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annots <- c('hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)
DefaultAssay(NK_act_full) <- "MACS2"

# Global overview of peak distribution
all_features <- ClosestFeature(NK_act_full, regions = rownames(NK_act_full[["MACS2"]]), annotation = annotations)

# Counts types
all_features %>% count(type) -> Peak_type_numbers
# Show as pie chart
pdf(file = "Plots/SupplFig4_ATAC_peak_types.pdf")
Peak_type_numbers %>% 
  ggplot(aes(x = "", y = n, fill = type))+
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette = "Accent")+
  theme_void()
dev.off()





##### Remove cells that do not cluster by conditions and come mainly from one donor
DefaultAssay(NK_act_full) <- "MACS2"



pdf(file = "Plots/SupplFig4_Dimplots.pdf", paper = "a4r")
DimPlot(NK_act_full, group.by = "condition", pt.size = 0.4, cols = pairedpal[c(3,4,7,8)])&NoLegend()
DimPlot(NK_act_full,  group.by = "donor", shuffle = T, pt.size = 0.4, cols = clusterpal[c(5,3)])&NoLegend()
DimPlot(NK_act_full, group.by = "seurat_clusters", cols = c("black", "grey", "grey56"), pt.size = 0.4)&NoLegend()
dev.off()


FeaturePlot(NK_act_full, features = c("CD56", "CD16", "anti-Biotin", "CD57",
                                 "KIR2DL1-S1-S3-S5", "KIR2DL2-L3"), cols = colorscale, min.cutoff = "q5", max.cutoff = "q99", order = TRUE, ncol =3)&NoAxes()&NoLegend()


# Subset cells represented from both donors & clustering by condition: Cluster 0
NK_act <- subset(NK_act_full, idents = c("0"))

NK_act <- RunUMAP(object = NK_act, reduction = 'harmony', dims = 2:30, n.neighbors = 30)
DimPlot(NK_act, group.by = "condition", cols = pairedpal[c(1,2,5,6)])

NK_act <- FindNeighbors(object = NK_act, reduction = 'harmony', dims = 2:30, k.param = 30)
NK_act <- FindClusters(object = NK_act, verbose = FALSE, algorithm = 3, resolution = 0.2)

DimPlot(NK_act)

# Subcluster cluster 0 to resolve peptide effect
temp <- subset(NK_act, idents = 0)
temp <- FindNeighbors(object = temp, reduction = 'harmony', dims = 2:30)
temp <- FindClusters(object = temp, verbose = FALSE, algorithm = 3, resolution = 0.2)
DimPlot(temp)
Idents(NK_act, cells = WhichCells(temp, idents = 1)) <- "Peptide+Cytokines"

# Isolate PQS/LFL conditions to resolve peptide effect also here (masked by strong cytokine effect otherwise)
NK_peptide <- subset(NK_act, subset = condition%in%c("PQS", "LFL"))


# Cluster at high resolution to catch LFL-enriched, NKG2A- CD137+ cluster
DimPlot(NK_peptide, group.by = "condition", cols = pairedpal[c(1,2,5,6)])
FeaturePlot(NK_peptide, features = c("anti-Biotin", "CD137"), cols = colorscale, min.cutoff = "q1", max.cutoff = "q99")
NK_peptide <- FindNeighbors(NK_peptide, reduction = "harmony", dims = 2:30)
NK_peptide <- FindClusters(NK_peptide, verbose = FALSE, algorithm = 3, resolution = 1)
DimPlot(NK_peptide, label = TRUE)


NK_peptide <- RenameIdents(NK_peptide, "2" = "Peptide")
VlnPlot(NK_peptide, features = c("anti-Biotin", "CD137"), log = TRUE)
Idents(NK_act, cells = WhichCells(NK_peptide, idents = "Peptide")) <- "Peptide"
NK_act <- RenameIdents(NK_act, "0" = "Cytokines", "1" = "Control")

# Reorder levels
levels(NK_act) <- c("Control", "Peptide", "Cytokines", "Peptide+Cytokines")
NK_act$annotation <- Idents(NK_act)

# Check composition
condition_df <- FetchData(NK_act, vars = c("condition", "annotation"))
condition_freq <- condition_df %>% group_by(annotation) %>% dplyr::count(condition)
condition_freq <- bind_rows(condition_freq, 
                            data_frame(annotation = rep("Total", 4),condition_df %>% dplyr::count(condition)))
condition_freq <- condition_freq %>% group_by(annotation) %>%  mutate(freq = n/sum(n))

condition_freq %>% 
  ggplot(aes(x=annotation, y=freq, fill=condition)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = pairedpal[c(1,2,5,6)]) +
  ylab("Frequency (%)")+
  theme_classic()


donor_df <- FetchData(NK_act, vars = c("donor", "annotation"))
donor_freq <- donor_df %>% group_by(annotation) %>% dplyr::count(donor)
donor_freq <- bind_rows(donor_freq, 
                            data_frame(annotation = rep("Total", 2),donor_df %>% dplyr::count(donor)))
donor_freq <- donor_freq %>% group_by(annotation) %>%  mutate(freq = n/sum(n))

donor_freq %>% 
  ggplot(aes(x=annotation, y=freq, fill=donor)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = pairedpal[c(1,2,5,6)]) +
  ylab("Frequency (%)")+
  theme_classic()
# Export
saveRDS(NK_act, file = paste0(path, "NK_act_dim"))



