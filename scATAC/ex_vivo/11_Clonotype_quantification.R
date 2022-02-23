# Quantification of clonotype to cluster association, identification of
# clonotype-specific peaks


library(BuenColors)
library(ggrastr)
library(SummarizedExperiment)
library(dplyr)
library(data.table)

# Load data
CMVpos2_full <- readRDS(paste0(path, "CMVpos2_full_alleles"))
CMVpos3_full <- readRDS(paste0(path, "CMVpos3_full_alleles"))
CMVpos4_full <- readRDS(paste0(path, "CMVpos4_full_alleles"))
CMVpos1_mito <- readRDS(paste0(path, "CMVpos1_mito"))

# Plot some representative mutations
DefaultAssay(CMVpos2_full) <- "alleles"

FeaturePlot(
  object = CMVpos2_full,
  features = c("5590G>A", "11157T>C", "4933T>C", "1415G>A"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2
)&NoAxes()

DefaultAssay(CMVpos3_full) <- "alleles"
FeaturePlot(
  object = CMVpos3_full,
  features = c("14422T>C","14165A>G", "13976A>G", "9943A>G"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2
)&NoAxes()


DefaultAssay(CMVpos4_full) <- "alleles"
FeaturePlot(
  object = CMVpos4_full,
  features = c("2983G>A", "9343G>A", "12396T>C", "12093T>C"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2
)&NoAxes()

DefaultAssay(CMVpos1_mito) <- "alleles"
FeaturePlot(
  object = CMVpos1_mito,
  features = c("13710A>G", "7395T>G", "4048G>A", "7621T>C"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2
)&NoAxes()




# Look at chi-square associations between clones and the cell type cluster

# Import data
All_NK_clone_CMVpos2 <- read.csv(file = paste0(path, "All_NK_clone_CMVpos2"), row.names = 1)
Adaptive_clone_CMVpos2 <- read.csv(file = paste0(path, "Adaptive_clone_CMVpos2"), row.names = 1)
Conventional_clone_CMVpos2 <- read.csv(file = paste0(path, "Conventional_clone_CMVpos2"), row.names = 1)

All_NK_clone_CMVpos3 <- read.csv(file = paste0(path, "All_NK_clone_CMVpos3"), row.names = 1)
Adaptive_clone_CMVpos3 <- read.csv(file = paste0(path, "Adaptive_clone_CMVpos3"), row.names = 1)
Conventional_clone_CMVpos3 <- read.csv(file = paste0(path, "Conventional_clone_CMVpos3"), row.names = 1)

All_NK_clone_CMVpos4 <- read.csv(file = paste0(path, "All_NK_clone_CMVpos4"), row.names = 1)
Adaptive_clone_CMVpos4 <- read.csv(file = paste0(path, "Adaptive_clone_CMVpos4"), row.names = 1)
Conventional_clone_CMVpos4 <- read.csv(file = paste0(path, "Conventional_clone_CMVpos4"), row.names = 1)

All_NK_clone_CMVpos1 <- read.csv(file = paste0(path, "All_NK_clone_CMVpos1"), row.names = 1)
Adaptive_clone_CMVpos1 <- read.csv(file = paste0(path, "Adaptive_clone_CMVpos1"), row.names = 1)
Conventional_clone_CMVpos1 <- read.csv(file = paste0(path, "Conventional_clone_CMVpos1"), row.names = 1)



# Function from Caleb Lareau's paper
get_chisquare_stats <- function(clone_df_input, seed = 0){
  
  if(seed == 0){
    clone_df <- clone_df_input
  } else {
    set.seed(seed)
    clone_df_input$mito_cluster <- sample(clone_df_input$mito_cluster)
    clone_df <- clone_df_input
  }
  
  # Pull the clone abundance
  clones <- clone_df %>% group_by(mito_cluster) %>% summarize(n=n()) %>% dplyr::filter(n >= 5) %>% pull(mito_cluster)
  chromatin_clusters <- as.character(sort(unique(clone_df$Clusters)))
  n_global <- clone_df %>% group_by(Clusters) %>% summarize(prop = n() / dim(clone_df)[1], n = n()) %>% pull(n)
  prop_global <- n_global/sum(n_global)
  
  # Loop over the clones and compute statistics
  lapply(clones, function(clone){
    
    ndf <- clone_df %>% dplyr::filter(mito_cluster == clone) %>%
      group_by(Clusters) %>% summarize(n = n())
    
    # Make a named vector of observations per cluster
    n_vec <- ndf$n; names(n_vec) <- ndf$Clusters
    n_clone <- unname(n_vec[chromatin_clusters])
    n_clone <- ifelse(is.na(n_clone), 0, n_clone)
    
    cs <- chisq.test(n_clone, p = prop_global)
    
    data.frame(
      clone_id = clone,
      CS_stat = unname(cs$statistic),
      CS_pvalue = unname(cs$p.value),
      n_cells = sum(ndf$n)
    )
  }) %>% rbindlist() %>% data.frame() -> clone_stat_df
  clone_stat_df$FDR <- p.adjust(clone_stat_df$CS_pvalue, "fdr")
  clone_stat_df$log10FDR <- -1*log10(clone_stat_df$FDR)
  clone_stat_df$log10pvalue <- -1*log10(clone_stat_df$CS_pvalue)
  
  print(length(clones))
  clone_stat_df %>% arrange(desc(log10pvalue)) %>% mutate(rank = 1:n(), seed = seed)
}

# Do it for both observed and permuted
All_NK_obs_CMVpos1 <- get_chisquare_stats(All_NK_clone_CMVpos1); All_NK_obs_CMVpos1$what <- "Observed"
All_NK_perm_CMVpos1 <- get_chisquare_stats(All_NK_clone_CMVpos1, 1); All_NK_perm_CMVpos1$what <- "Permuted"
Conventional_NK_obs_CMVpos1 <- get_chisquare_stats(Conventional_clone_CMVpos1); Conventional_NK_obs_CMVpos1$what <- "Observed"
Conventional_NK_perm_CMVpos1 <- get_chisquare_stats(Conventional_clone_CMVpos1, 1); Conventional_NK_perm_CMVpos1$what <- "Permuted"
Adaptive_NK_obs_CMVpos1 <- get_chisquare_stats(Adaptive_clone_CMVpos1); Adaptive_NK_obs_CMVpos1$what <- "Observed"
Adaptive_NK_perm_CMVpos1 <- get_chisquare_stats(Adaptive_clone_CMVpos1, 1); Adaptive_NK_perm_CMVpos1$what <- "Permuted"

All_NK_obs_CMVpos2 <- get_chisquare_stats(All_NK_clone_CMVpos2); All_NK_obs_CMVpos2$what <- "Observed"
All_NK_perm_CMVpos2 <- get_chisquare_stats(All_NK_clone_CMVpos2, 1); All_NK_perm_CMVpos2$what <- "Permuted"
Conventional_NK_obs_CMVpos2 <- get_chisquare_stats(Conventional_clone_CMVpos2); Conventional_NK_obs_CMVpos2$what <- "Observed"
Conventional_NK_perm_CMVpos2 <- get_chisquare_stats(Conventional_clone_CMVpos2, 1); Conventional_NK_perm_CMVpos2$what <- "Permuted"
Adaptive_NK_obs_CMVpos2 <- get_chisquare_stats(Adaptive_clone_CMVpos2); Adaptive_NK_obs_CMVpos2$what <- "Observed"
Adaptive_NK_perm_CMVpos2 <- get_chisquare_stats(Adaptive_clone_CMVpos2, 1); Adaptive_NK_perm_CMVpos2$what <- "Permuted"

All_NK_obs_CMVpos3 <- get_chisquare_stats(All_NK_clone_CMVpos3); All_NK_obs_CMVpos3$what <- "Observed"
All_NK_perm_CMVpos3 <- get_chisquare_stats(All_NK_clone_CMVpos3, 1); All_NK_perm_CMVpos3$what <- "Permuted"
Conventional_NK_obs_CMVpos3 <- get_chisquare_stats(Conventional_clone_CMVpos3); Conventional_NK_obs_CMVpos3$what <- "Observed"
Conventional_NK_perm_CMVpos3 <- get_chisquare_stats(Conventional_clone_CMVpos3, 1); Conventional_NK_perm_CMVpos3$what <- "Permuted"
Adaptive_NK_obs_CMVpos3 <- get_chisquare_stats(Adaptive_clone_CMVpos3); Adaptive_NK_obs_CMVpos3$what <- "Observed"
Adaptive_NK_perm_CMVpos3 <- get_chisquare_stats(Adaptive_clone_CMVpos3, 1); Adaptive_NK_perm_CMVpos3$what <- "Permuted"

All_NK_obs_CMVpos4 <- get_chisquare_stats(All_NK_clone_CMVpos4); All_NK_obs_CMVpos4$what <- "Observed"
All_NK_perm_CMVpos4 <- get_chisquare_stats(All_NK_clone_CMVpos4, 1); All_NK_perm_CMVpos4$what <- "Permuted"
Conventional_NK_obs_CMVpos4 <- get_chisquare_stats(Conventional_clone_CMVpos4); Conventional_NK_obs_CMVpos4$what <- "Observed"
Conventional_NK_perm_CMVpos4 <- get_chisquare_stats(Conventional_clone_CMVpos4, 1); Conventional_NK_perm_CMVpos4$what <- "Permuted"
Adaptive_NK_obs_CMVpos4 <- get_chisquare_stats(Adaptive_clone_CMVpos4); Adaptive_NK_obs_CMVpos4$what <- "Observed"
Adaptive_NK_perm_CMVpos4 <- get_chisquare_stats(Adaptive_clone_CMVpos4, 1); Adaptive_NK_perm_CMVpos4$what <- "Permuted"


# Bind donors together
All_NK_obs <- rbind(All_NK_obs_CMVpos2, All_NK_obs_CMVpos3, All_NK_obs_CMVpos4)
All_NK_perm <- rbind(All_NK_perm_CMVpos2, All_NK_perm_CMVpos3, All_NK_perm_CMVpos4)
Conventional_NK_obs <- rbind(Conventional_NK_obs_CMVpos2, Conventional_NK_obs_CMVpos3, Conventional_NK_obs_CMVpos4)
Conventional_NK_perm <- rbind(Conventional_NK_perm_CMVpos2, Conventional_NK_perm_CMVpos3, Conventional_NK_perm_CMVpos4)
Adaptive_NK_obs <- rbind(Adaptive_NK_obs_CMVpos2, Adaptive_NK_obs_CMVpos3, Adaptive_NK_obs_CMVpos4)
Adaptive_NK_perm <- rbind(Adaptive_NK_perm_CMVpos2, Adaptive_NK_perm_CMVpos3, Adaptive_NK_perm_CMVpos4)

# Re-rank clones so that donors can be plotted together
All_NK_obs %>% arrange(FDR) %>% mutate(rank = 1:n(), seed = seed) -> All_NK_obs
All_NK_perm %>% arrange(FDR) %>% mutate(rank = 1:n(), seed = seed) -> All_NK_perm
Conventional_NK_obs %>% arrange(FDR) %>% mutate(rank = 1:n(), seed = seed) -> Conventional_NK_obs
Conventional_NK_perm %>% arrange(FDR) %>% mutate(rank = 1:n(), seed = seed) -> Conventional_NK_perm
Adaptive_NK_obs %>% arrange(FDR) %>% mutate(rank = 1:n(), seed = seed) -> Adaptive_NK_obs
Adaptive_NK_perm %>% arrange(FDR) %>% mutate(rank = 1:n(), seed = seed) -> Adaptive_NK_perm

p1 <- ggplot(rbind(Conventional_NK_obs, Conventional_NK_perm) %>% arrange(desc(what)), aes(x = rank, y = log10FDR, color = what)) + 
  geom_point(size = 1.2)+
  scale_color_manual(values = c("black", "lightgrey"))  + 
  labs(x = "Rank sorted clones", y = "-log10 FDR") + 
  pretty_plot(fontsize = 14) + L_border() +  
  scale_y_continuous(limits = c(0,80)) +
  theme(legend.position = "none") + ggtitle("Conventional NK")

p2 <- ggplot(rbind(All_NK_obs, All_NK_perm) %>% arrange(desc(what)), aes(x = rank, y = log10FDR, color = what)) + 
  geom_point(size = 1.2)+
  scale_color_manual(values = c("black", "lightgrey"))  + 
  labs(x = "Rank sorted clones", y = "-log10 FDR")  +
  pretty_plot(fontsize = 14) + L_border() +  
  scale_y_continuous(limits = c(0,80)) +
  theme(legend.position = "none") + ggtitle("All NK")
p3 <- ggplot(rbind(Adaptive_NK_obs, Adaptive_NK_perm) %>% arrange(desc(what)), aes(x = rank, y = log10FDR, color = what)) + 
  geom_point(size = 1.2)+
  scale_color_manual(values = c("black", "lightgrey"))  + 
  labs(x = "Rank sorted clones", y = "-log10 FDR")  +
  pretty_plot(fontsize = 14) + L_border() +
  scale_y_continuous(limits = c(0,80)) +
  theme(legend.position = "none") + ggtitle("Adaptive NK")
p2+p1+p3

# Export to plot in graphpad (duh!)
write.csv(All_NK_obs, file = "Plots/All_NK_obs.csv")
write.csv(All_NK_perm, file = "Plots/All_NK_perm.csv")
write.csv(Conventional_NK_obs, file = "Plots/Conventional_NK_obs.csv")
write.csv(Conventional_NK_perm, file = "Plots/Conventional_NK_perm.csv")
write.csv(Adaptive_NK_obs, file = "Plots/Adaptive_NK_obs.csv")
write.csv(Adaptive_NK_perm, file = "Plots/Adaptive_NK_perm.csv")



# Manually curate clones for downstream analysis

# Get significantly associated clonotypes
sig_clones_CMVpos2 <- All_NK_obs_CMVpos2 %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)
sig_clones_CMVpos3 <- All_NK_obs_CMVpos3 %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)
sig_clones_CMVpos4 <- All_NK_obs_CMVpos4 %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)
sig_clones_CMVpos1 <- All_NK_obs_CMVpos1 %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)

# Inspect which clones are associated to adaptive population
#CMVpos2
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "3"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "8"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "56"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "20"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "13"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "27"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "19"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "30"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "33"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "23"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "44"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "12"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "98"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "14"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "26"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "39"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "24"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "11"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "29"))#no
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "72"))#yes
DimPlot(CMVpos2_full, cells.highlight = WhichCells(CMVpos2_full, expression = clonotype == "36"))#no


sig_clones_CMVpos2_curated <- paste0(c("3", "8", "56", "20", "27",
                                       "33", "44", "98", "11", "72", "44"), "CMVpos2")

#CMVpos3
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "6"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "11"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "7"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "19"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "8"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "18"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "15"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "34"))#yes
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "35"))#no
DimPlot(CMVpos3_full, cells.highlight = WhichCells(CMVpos3_full, expression = clonotype == "25"))#yes


sig_clones_CMVpos3_curated <- paste0(c("6", "11", "7", "19", "8", "18",
                                       "15", "34", "25"), "CMVpos3")

#CMVpos4 
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "66"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "20"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "7"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "8"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "13"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "11"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "47"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "53"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "44"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "14"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "69"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "24"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "23"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "40"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "21"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "37"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "46"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "22"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "56"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "62"))#yes
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "25"))#no
DimPlot(CMVpos4_full, cells.highlight = WhichCells(CMVpos4_full, expression = clonotype == "9"))#yes



sig_clones_CMVpos4_curated <- paste0(c("66", "20", "7", "8", "47",
                                       "53", "14", "69", "40", "21",
                                       "37", "22", "62", "9"), "CMVpos4")

#CMVpos1  
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "13"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "2"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "5"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "7"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "9"))#no
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "18"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "15"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "10"))#no
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "21"))#no
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "12"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "14"))#no
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "17"))#yes
DimPlot(CMVpos1_mito, cells.highlight = WhichCells(CMVpos1_mito, expression = clonotype == "16"))#no


sig_clones_CMVpos1_curated <- paste0(c("13", "2", "5", "7", "18",
                                       "15", "12", "17"), "CMVpos1")

# Plot KIR3DL1 and "5590G>A" for CMVpos2
CMVpos2_mtASAP3 <- subset(CMVpos2_full, subset = experiment == "mtASAP3")
VlnPlot(CMVpos2_mtASAP3, features = "KIR3DL1")&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[1:3]))
VlnPlot(CMVpos2_full, features = "5590G>A")&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[1:3]))

# Plot KIR2DL1 and "13710A>G" for CMVpos1
VlnPlot(CMVpos1_mito, features = "KIR2DL1-S1-S3-S5", y.max = 3)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))
VlnPlot(CMVpos1_mito, features = "13710A>G")&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal))

#### Peak-clonotype association analysis #####
# Functions modified from Caleb Lareau's Nature Biotech
getX2stats <- function(y, cid){
  obs <- chisq.test(cid, y)
  perm <- chisq.test(cid, sample(y))
  data.frame(
    obs_p = obs$p.value,
    obs_x2 = obs$statistic,
    perm_p = perm$p.value,
    perm_x2 = perm$statistic
  )
}


# Perform Chisquare test for Peak clonotype association
peak_clonotype_association <- function(seu, key){
  # Pull data from Seurat container and format as for Caleb's script
  peaks <- as.data.frame(granges(seu))[,1:3]
  colnames(peaks) <- c("chr", "start", "end")
  peaks_gr <- makeGRangesFromDataFrame(peaks)
  
  barcodes <- colnames(seu)
  data <- GetAssayData(seu, slot = "counts")
  
  # Assemble binary matrix
  mat <- Matrix::sparseMatrix(i = c(summary(data)[[1]],length(peaks)), 
                              j = c(summary(data)[[2]],length(barcodes)),
                              x = c(as.numeric(summary(data)[[3]]) > 0, 0))
  
  colnames(mat) <- barcodes
  
  # Import cluster annotation from mtDNA
  boo <- Matrix::rowSums(mat) >= 10
  mat <- mat[boo,]
  clusters <- as.character(seu$clonotype2)
  peaks <- peaks[boo, ]
  
  # Reformat for faster subsetting
  sm <- data.table(summary(mat))
  
  # Look at all peaks
  set.seed(1)
  x2df <- lapply(1:dim(mat)[1], function(idx){ 
    if((idx %% 1000) == 0){
      print(idx)
    }
    Y = 1:length(clusters) %in% sm[i == idx][["j"]]
    data.frame(peaks[idx,], getX2stats(Y, clusters))
    
  }) %>% rbindlist() %>% data.frame()
  
  saveRDS(x2df, file = paste0(key,"_X2.rds"))
}


# Merge adaptive subclusters with clonotype information
seu_adaptive <- merge(subset(CMVpos2_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
                      c(subset(CMVpos3_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
                        subset(CMVpos4_full, idents = c("Adaptive1", "Adaptive2", "Adaptive3")),
                        subset(CMVpos1_mito, idents = c("Adaptive1", "Adaptive2", "Adaptive3",
                                                        "Adaptive4", "Adaptive5")))
)

seu_adaptive$clonotype2 <- paste0(seu_adaptive$clonotype, seu_adaptive$donor)

# Subset adaptive clonotypes for analysis of clonotype-peak association
seu_adaptive_clones <- subset(seu_adaptive, 
                              subset= clonotype2%in%c(sig_clones_CMVpos2_curated,
                                                           sig_clones_CMVpos4_curated,
                                                           sig_clones_CMVpos3_curated,
                                                      sig_clones_CMVpos1_curated))

peak_clonotype_association(seu_adaptive_clones, "merged_clones")

# Plot
dt <- readRDS("merged_clones_X2.rds")
dt <- dt %>% arrange(desc(obs_x2)) %>%
  mutate(rank = 1:n(), padj = p.adjust(obs_p), obs_log10p = -1*log10(obs_p))
dt$perm_x2 <- sort(dt$perm_x2, decreasing = TRUE)
dt$perm_log10p <- -1*log10(sort(dt$perm_p, decreasing = FALSE))

sum(dt$padj < 0.05)

p1 <- ggplot(dt, aes(x = rank, y = obs_x2, color = padj < 0.05)) + 
  geom_point(size = 1) + scale_color_manual(values = c("black", "firebrick"))  + 
  geom_point(inherit.aes = FALSE, data = dt, aes(x = rank, y = perm_x2), color = "lightgrey", size = 0.1) +
  labs(x = "Rank sorted peaks", y = "X2 Statistic") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")
p1


# Plot some of the significantly associated peaks for respective clones
# IL7R
ranges.show <- StringToGRanges("chr5-35853528-35854727")
ranges.show$color <- "orange"
CoveragePlot(
  object = seu_adaptive_clones,
  region = "chr5-35853528-35854727",
  extend.upstream = 10000,
  extend.downstream = 10000,
  scale.factor = 1e7,
  region.highlight = ranges.show,
  group.by = "clonotype2"
)

# GLI3
ranges.show <- StringToGRanges("chr7-42227427-42228268")
ranges.show$color <- "orange"
CoveragePlot(
  object = seu_adaptive_clones,
  region = "chr7-42227427-42228268",
  extend.upstream = 10000,
  extend.downstream = 10000,
  scale.factor = 1e7,
  region.highlight = ranges.show,
  group.by = "clonotype2"
)

# GLI3: Clone 7 fr. CMVpos4; IL7R: Clone 8 CMVpos2
seu_adaptive_clones$clone <- "rest"
seu_adaptive_clones$clone[seu_adaptive_clones$clonotype2=="7CMVpos4"] <- "7CMVpos4"
seu_adaptive_clones$clone[seu_adaptive_clones$clonotype2=="8CMVpos2"] <- "8CMVpos2"

CoveragePlot(
  object = seu_adaptive_clones,
  region = "chr7-42227427-42228268",
  extend.upstream = 10000,
  extend.downstream = 10000,
  scale.factor = 1e7,
  group.by = "clone"
)


CoveragePlot(
  object = seu_adaptive_clones,
  region = "chr5-35853528-35854727",
  extend.upstream = 10000,
  extend.downstream = 10000,
  scale.factor = 1e7,
  group.by = "clone"
)


# Some unique peak examples

#CMVpos 4
seu_adaptive_clones$clone <- "rest"
seu_adaptive_clones$clone[seu_adaptive_clones$clonotype2=="7CMVpos4"] <- "7CMVpos4"


CoveragePlot(
  object = seu_adaptive_clones,
  region = "chr17-82550724-82551104",
  extend.upstream = 10000,
  extend.downstream = 10000,
  scale.factor = 1e7,
  group.by = "clone"
)

#CMVpos2
seu_adaptive_clones$clone <- "rest"
seu_adaptive_clones$clone[seu_adaptive_clones$clonotype2=="3CMVpos2"] <- "3CMVpos2"

CoveragePlot(
  object = seu_adaptive_clones,
  region = "chr2-208358929-208359964",
  extend.upstream = 10000,
  extend.downstream = 10000,
  scale.factor = 1e7,
  group.by = "clone"
)


