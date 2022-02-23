# Longitudinal analysis of P2,P3,P4

# Dimplots for both time points
CMVpos2_list <- SplitObject(CMVpos2_full, split.by = "experiment") 
CMVpos4_list <- SplitObject(CMVpos4_full, split.by = "experiment")
CMVpos3_list <- SplitObject(CMVpos3_full, split.by = "experiment") 

p1 <- lapply(X = CMVpos2_list, FUN = function(x) {
  DimPlot(x, reduction = "umap", cols = c(clusterpal[1:3], adaptivepal[1:3]))&NoLegend()&scale_y_continuous(limits = c(-4.5, 6.5), breaks = c(-4,-2,0,2,4,6))
})

p2 <- lapply(X = CMVpos4_list, FUN = function(x) {
  DimPlot(x, reduction = "umap", cols = c(clusterpal[1:3], adaptivepal[4:6]))&NoLegend()
})

p3 <- lapply(X = CMVpos3_list, FUN = function(x) {
  DimPlot(x, reduction = "umap", cols = c(clusterpal[1:3], adaptivepal[7:9]))&NoLegend()
})

pdf("Plots/new/Fig7_dimplot.pdf", width = 8, height = 8)
p1
p2
p3
dev.off()


#### Heatmaps of subcluster-specific peaks over time ####
DefaultAssay(CMVpos2_full) <- "MACS2"
DefaultAssay(CMVpos4_full) <- "MACS2"
DefaultAssay(CMVpos3_full) <- "MACS2"

avg_exp_tidy_df <- function(seu, features, assay){
  avg_exp <- AverageExpression(seu, features = features, assays = assay)
  avg_exp <- as.data.frame(avg_exp[[as.character(assay)]])
  avg_exp$region <- rownames(avg_exp)
  avg_exp_tidy <- gather(avg_exp, key = "cluster", value = "avg_expression", -region)
  return(avg_exp_tidy)
}

# CMVpos2
average_acc_CMVpos2_t0 <- avg_exp_tidy_df(subset(CMVpos2_full, subset = experiment == "CMVpos2"),
                                          features = sub_adaptive_DAP_CMVpos2$gene, assay = "MACS2")
average_acc_CMVpos2_t1 <- avg_exp_tidy_df(subset(CMVpos2_full, subset = experiment == "mtASAP3"),
                                          features = sub_adaptive_DAP_CMVpos2$gene, assay = "MACS2")

average_acc_CMVpos2_t0$cluster <- paste0(average_acc_CMVpos2_t0$cluster, "t0")
average_acc_CMVpos2_t1$cluster <- paste0(average_acc_CMVpos2_t1$cluster, "t1")

average_acc_CMVpos2 <- rbind(average_acc_CMVpos2_t0, average_acc_CMVpos2_t1)

# CMVpos4
average_acc_CMVpos4_t0 <- avg_exp_tidy_df(subset(CMVpos4_full, subset = experiment == "mtASAP2"),
                                          features = sub_adaptive_DAP_CMVpos4$gene, assay = "MACS2")
average_acc_CMVpos4_t1 <- avg_exp_tidy_df(subset(CMVpos4_full, subset = experiment == "mtASAP3"),
                                          features = sub_adaptive_DAP_CMVpos4$gene, assay = "MACS2")


average_acc_CMVpos4_t0$cluster <- paste0(average_acc_CMVpos4_t0$cluster, "t0")
average_acc_CMVpos4_t1$cluster <- paste0(average_acc_CMVpos4_t1$cluster, "t1")

average_acc_CMVpos4 <- rbind(average_acc_CMVpos4_t0, average_acc_CMVpos4_t1)

# CMVpos3
average_acc_CMVpos3_t0 <- avg_exp_tidy_df(subset(CMVpos3_full, subset = experiment == "mtASAP2"),
                                          features = sub_adaptive_DAP_CMVpos3$gene, assay = "MACS2")
average_acc_CMVpos3_t1 <- avg_exp_tidy_df(subset(CMVpos3_full, subset = experiment == "mtASAP5"),
                                          features = sub_adaptive_DAP_CMVpos3$gene, assay = "MACS2")


average_acc_CMVpos3_t0$cluster <- paste0(average_acc_CMVpos3_t0$cluster, "t0")
average_acc_CMVpos3_t1$cluster <- paste0(average_acc_CMVpos3_t1$cluster, "t1")

average_acc_CMVpos3 <- rbind(average_acc_CMVpos3_t0, average_acc_CMVpos3_t1)


# Plots
average_acc_CMVpos4 %>% filter(cluster%in%c("Adaptive1t0", "Adaptive2t0", "Adaptive3t0",
                                            "Adaptive1t1", "Adaptive2t1", "Adaptive3t1")) %>%
  tidy_heatmap(
  rows = cluster,
  columns = region,
  values = avg_expression,
  cluster_rows = TRUE,
  #cluster_cols = TRUE,
  scale = "column",
  colors =  rev(mapal)
)

average_acc_CMVpos4_t1 %>% filter(cluster%in%c("Adaptive1", "Adaptive2", "Adaptive3")) %>%  tidy_heatmap(
  rows = cluster,
  columns = region,
  values = avg_expression,
  #cluster_rows = TRUE,
  #cluster_cols = TRUE,
  scale = "column",
  colors =  rev(mapal)
)

average_acc_CMVpos2_t0 %>% filter(cluster%in%c("Adaptive1", "Adaptive2", "Adaptive3")) %>%  tidy_heatmap(
  rows = cluster,
  columns = region,
  values = avg_expression,
  #cluster_rows = TRUE,
  #cluster_cols = TRUE,
  scale = "column",
  colors =  rev(mapal)
)
average_acc_CMVpos2_t1 %>% filter(cluster%in%c("Adaptive1", "Adaptive2", "Adaptive3")) %>%  tidy_heatmap(
  rows = cluster,
  columns = region,
  values = avg_expression,
  #cluster_rows = TRUE,
  #cluster_cols = TRUE,
  scale = "column",
  colors =  rev(mapal)
)

#
average_acc_CMVpos3 %>% dplyr::filter(cluster%in%c("Adaptive1t0", "Adaptive2t0", "Adaptive3t0",
                                                   "Adaptive1t1", "Adaptive2t1", "Adaptive3t1")) %>%
  tidy_heatmap(
    rows = cluster,
    columns = region,
    values = avg_expression,
    cluster_rows = TRUE,
    #cluster_cols = TRUE,
    scale = "column",
    colors =  rev(mapal)
  )

# Subcluster-specific peaks
# CMVpos2
CoveragePlot(
  object = subset(CMVpos2_full, subset = experiment == "CMVpos2"),
  region = "chr21-31087852-31088428",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[1:3]))

CoveragePlot(
  object = subset(CMVpos2_full, subset = experiment == "mtASAP3"),
  region = "chr21-31087852-31088428",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[1:3]))


# CMVpos4
CoveragePlot(
  object = subset(CMVpos4_full, subset = experiment == "mtASAP2"),
  region = "chr17-67398748-67399002",
  extend.upstream = 4000,
  extend.downstream = 4000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[4:6]))

CoveragePlot(
  object = subset(CMVpos4_full, subset = experiment == "mtASAP3"),
  region = "chr17-67398748-67399002",
  extend.upstream = 4000,
  extend.downstream = 4000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[4:6]))

# CMVpos3
CoveragePlot(
  object = subset(CMVpos3_full, subset = experiment == "mtASAP2"),
  region = "chr16-81678354-81679431",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[7:9]))

CoveragePlot(
  object = subset(CMVpos3_full, subset = experiment == "mtASAP5"),
  region = "chr16-81678354-81679431",
  extend.upstream = 5000,
  extend.downstream = 5000,
  scale.factor = 1e7,
)&scale_fill_manual(values = c(clusterpal[1:3], adaptivepal[7:9]))




### Clonotype stability

# Plot some representative mutations
DefaultAssay(CMVpos2_full) <- "alleles"
FeaturePlot(
  object = CMVpos2_full,
  features = c("5590G>A", "11157T>C", "4933T>C", "1415G>A"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2,
  split.by = "experiment"
)&NoAxes()

DefaultAssay(CMVpos3_full) <- "alleles"
FeaturePlot(
  object = CMVpos3_full,
  features = c("14422T>C","14165A>G", "13976A>G", "9943A>G"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2,
  split.by = "experiment"
)&NoAxes()

DefaultAssay(CMVpos4_full) <- "alleles"
FeaturePlot(
  object = CMVpos4_full,
  features = c("2983G>A", "9343G>A", "12396T>C", "15876T>C"),
  order = TRUE,
  cols = c("grey", "darkred"),
  coord.fixed = T,
  ncol = 2,
  split.by = "experiment"
)&NoAxes()


##### Quantifiy clonotype frequencies #####

get_clonotype_counts <- function(seu, clonotypes, donor){
  foo <- data.frame(cluster = Idents(seu), clonotype = paste0(seu$clonotype, donor), experiment = seu$experiment)
  # Filter for adaptive cells
  foo %>% dplyr::filter(cluster%in%c("Adaptive1", "Adaptive2", "Adaptive3")) -> foo
  # Filter for defined clonotypes and count their occurence
  foo %>% group_by(experiment) %>%  dplyr::filter(clonotype%in%clonotypes) %>% dplyr::count(clonotype) -> out
  # Add number of cells not in clonotypes
  clone_number <- out %>% group_by(experiment) %>% summarise(sum(n))
  cell_number <- foo %>% group_by(experiment) %>% summarise(n())
  non_clonotypes <- data.frame(experiment = clone_number[,1], clonotype = "non_clonotypes", n = cell_number[,2]-clone_number[,2])
  colnames(non_clonotypes) <- colnames(out)
  out <- rbind(out, non_clonotypes)
  return(out)
}



# Quantification
clone_counts_CMVpos2 <- get_clonotype_counts(CMVpos2_full, sig_clones_CMVpos2_curated, "CMVpos2")
clone_counts_CMVpos3 <- get_clonotype_counts(CMVpos3_full, sig_clones_CMVpos3_curated, "CMVpos3")
clone_counts_CMVpos4 <- get_clonotype_counts(CMVpos4_full, sig_clones_CMVpos4_curated, "CMVpos4")

# Calculate frequencies
clone_counts_CMVpos2 <- clone_counts_CMVpos2 %>% group_by(experiment) %>%  mutate(freq = n/sum(n))
clone_counts_CMVpos3 <- clone_counts_CMVpos3 %>%  group_by(experiment) %>%  mutate(freq = n/sum(n))
clone_counts_CMVpos4 <- clone_counts_CMVpos4 %>% group_by(experiment) %>% mutate(freq = n/sum(n))




# Reorder by clone size and plot as bargraph, include only clonotypes with
# at least 

clone_counts_CMVpos2 %>% 
  mutate(clonotype = forcats::fct_reorder(.f = clonotype, .x = freq, .fun = mean)) %>% 
  dplyr::filter(clonotype != "non_clonotypes") %>% 
  ggplot(aes(x=experiment, y=freq*100, fill=clonotype)) +
  geom_bar(stat="identity", width=0.7, color="white") +
  scale_y_continuous(limits = c(0,12))+
  scale_fill_manual(values = large_pal)+
  ylab("Clonotype Frequency (%)")+
  theme_classic()&NoLegend()

clone_counts_CMVpos3 %>% 
  mutate(clonotype = forcats::fct_reorder(.f = clonotype, .x = freq, .fun = mean)) %>% 
  dplyr::filter(clonotype != "non_clonotypes") %>% 
  ggplot(aes(x=experiment, y=freq*100, fill=clonotype)) +
  geom_bar(stat="identity", width=0.7, color="white") +
  scale_y_continuous(limits = c(0,12))+
  scale_fill_manual(values = large_pal)+
  ylab("Clonotype Frequency (%)")+
  theme_classic()&NoLegend()


clone_counts_CMVpos4 %>% 
  mutate(clonotype = forcats::fct_reorder(.f = clonotype, .x = freq, .fun = mean)) %>% 
  dplyr::filter(clonotype != "non_clonotypes") %>% 
  ggplot(aes(x=experiment, y=freq*100, fill=clonotype)) +
  geom_bar(stat="identity", width=0.7, color="white") +
  scale_y_continuous(limits = c(0,12))+
  scale_fill_manual(values = large_pal)+
  ylab("Clonotype Frequency (%)")+
  theme_classic()&NoLegend()



# Test for significance
library(tidyr)
# CMVpos2
clone_df_CMVpos2 <- data.frame(CMVpos2= ungroup(clone_counts_CMVpos2) %>% 
                                 dplyr::filter(clonotype!= "non_clonotypes") %>%  
                                 dplyr::filter(experiment=="CMVpos2") %>% dplyr::select(n),
                               
                               mtASAP3= ungroup(clone_counts_CMVpos2) %>% 
                                 dplyr::filter(clonotype!= "non_clonotypes") %>% 
                                 dplyr::filter(experiment=="mtASAP3") %>% dplyr::select(n)
)
rownames(clone_df_CMVpos2) <- ungroup(clone_counts_CMVpos2) %>% 
  dplyr::filter(clonotype!= "non_clonotypes") %>% 
  dplyr::filter(experiment=="mtASAP3") %>% dplyr::pull(clonotype)

colnames(clone_df_CMVpos2) <- c("CMVpos2", "mtASAP3")

test <- fisher.test(clone_df_CMVpos2, simulate.p.value = T)
# CMVpos3
clone_df_CMVpos3 <- data.frame(mtASAP2= ungroup(clone_counts_CMVpos3) %>% 
                                 dplyr::filter(clonotype!= "non_clonotypes") %>%  
                                 dplyr::filter(experiment=="mtASAP2") %>% dplyr::select(n),
                               
                               mtASAP5= ungroup(clone_counts_CMVpos3) %>% 
                                 dplyr::filter(clonotype!= "non_clonotypes") %>% 
                                 dplyr::filter(experiment=="mtASAP5") %>% dplyr::select(n)
)
rownames(clone_df_CMVpos3) <- ungroup(clone_counts_CMVpos3) %>% 
  dplyr::filter(clonotype!= "non_clonotypes") %>% 
  dplyr::filter(experiment=="mtASAP5") %>% dplyr::pull(clonotype)

colnames(clone_df_CMVpos3) <- c("mtASAP2", "mtASAP5")

test1 <- fisher.test(clone_df_CMVpos3, simulate.p.value = T)

# CMvpos4

clone_df_CMVpos4 <- data.frame(CMVpos4= ungroup(clone_counts_CMVpos4) %>% 
                                 dplyr::filter(clonotype!= "non_clonotypes") %>%  
                                 dplyr::filter(experiment=="mtASAP2") %>% dplyr::select(n),
                               
                               mtASAP3= ungroup(clone_counts_CMVpos4) %>% 
                                 dplyr::filter(clonotype!= "non_clonotypes") %>% 
                                 dplyr::filter(experiment=="mtASAP3") %>% dplyr::select(n)
)
rownames(clone_df_CMVpos4) <- ungroup(clone_counts_CMVpos4) %>% 
  dplyr::filter(clonotype!= "non_clonotypes") %>% 
  dplyr::filter(experiment=="mtASAP2") %>% dplyr::pull(clonotype)

colnames(clone_df_CMVpos4) <- c("mtASAP2", "mtASAP3")

test2 <- fisher.test(clone_df_CMVpos4, simulate.p.value = T)

test
test1
test2

# Export counts
write.csv(clone_df_CMVpos2, file = "Plots/new/clone_counts_CMVpos2")
write.csv(clone_df_CMVpos3, file = "Plots/new/clone_counts_CMVpos3")
write.csv(clone_df_CMVpos4, file = "Plots/new/clone_counts_CMVpos4")

# ALternative from Caleb: KS test and density distribution

# Create clonotype time point df to match input to Caleb's script
get_clonotype_df <- function(seu, donor, clonotypes){
  foo <- data.frame(cluster = seu$annotation, clonotype = paste0(seu$clonotype, donor), experiment = seu$experiment)
  # Filter for adaptive cells
  foo %>% dplyr::filter(cluster%in%c("Adaptive1", "Adaptive2", "Adaptive3")) -> foo
  # Filter for defined clonotypes
  foo %>% group_by(experiment) %>% dplyr::filter(clonotype%in%paste0(clonotypes, donor)) %>% dplyr::select(clonotype, experiment)  -> out
  return(out)
}


# CMVpos2
set.seed(1234)
clone_df_in <- get_clonotype_df(CMVpos2_full, "CMVpos2", filtered_clonotypes_CMVpos2)


time_df <- clone_df_in %>%
  group_by(clonotype,experiment) %>% summarize(count = n())

dd <- reshape2::dcast(time_df, clonotype ~ experiment, value.var = "count", fill = 1)
dd$t0p <- dd$CMVpos2/sum(dd$CMVpos2)
dd$t1p <- dd$mtASAP3/sum(dd$mtASAP3)

observed_foldchange <- log2((dd$t0p )/(dd$t1p))


perm_df <- data.frame(clonotype = clone_df_in$clonotype,
                      experiment = sample(clone_df_in$experiment)) %>%
  group_by(clonotype,experiment) %>% summarize(count = n())
dd2 <- reshape2::dcast(perm_df, clonotype ~ experiment, value.var = "count", fill = 1)
dd2$t0p <- dd2$CMVpos2/sum(dd2$CMVpos2)
dd2$t1p <- dd2$mtASAP3/sum(dd2$mtASAP3)

permuted_foldchange <- log2((dd2$t0p )/(dd2$t1p))

density_df <- data.frame(
  logFC = c(observed_foldchange, permuted_foldchange),
  what = c(rep("observed", length(observed_foldchange)),
           rep("permuted", length(permuted_foldchange))
  ))
pDens <- ggplot(density_df, aes(x = logFC, fill = what)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values =c ("firebrick", "grey")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-4.5,4.5))+
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 FC", y = "empirical density")+
  theme(legend.position = "none")
pDens

ks.test(permuted_foldchange, observed_foldchange)

# CMVpos3
set.seed(1234)
clone_df_in <- get_clonotype_df(CMVpos3_full, "CMVpos3", filtered_clonotypes_CMVpos3)


time_df <- clone_df_in %>%
  group_by(clonotype,experiment) %>% summarize(count = n())

dd <- reshape2::dcast(time_df, clonotype ~ experiment, value.var = "count", fill = 1)
dd$t0p <- dd$mtASAP2/sum(dd$mtASAP2)
dd$t1p <- dd$mtASAP5/sum(dd$mtASAP5)

observed_foldchange <- log2((dd$t0p )/(dd$t1p))


perm_df <- data.frame(clonotype = clone_df_in$clonotype,
                      experiment = sample(clone_df_in$experiment)) %>%
  group_by(clonotype,experiment) %>% summarize(count = n())
dd2 <- reshape2::dcast(perm_df, clonotype ~ experiment, value.var = "count", fill = 1)
dd2$t0p <- dd2$mtASAP2/sum(dd2$mtASAP2)
dd2$t1p <- dd2$mtASAP5/sum(dd2$mtASAP5)

permuted_foldchange <- log2((dd2$t0p )/(dd2$t1p))

density_df <- data.frame(
  logFC = c(observed_foldchange, permuted_foldchange),
  what = c(rep("observed", length(observed_foldchange)),
           rep("permuted", length(permuted_foldchange))
  ))
pDens <- ggplot(density_df, aes(x = logFC, fill = what)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values =c ("firebrick", "grey")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-4.5,4.5))+
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 FC", y = "empirical density")+
  theme(legend.position = "none")
pDens

ks.test(permuted_foldchange, observed_foldchange)


# CMVpos4
set.seed(1234)
clone_df_in <- get_clonotype_df(CMVpos4_full, "CMVpos4", filtered_clonotypes_CMVpos4)


time_df <- clone_df_in %>%
  group_by(clonotype,experiment) %>% summarize(count = n())

dd <- reshape2::dcast(time_df, clonotype ~ experiment, value.var = "count", fill = 1)
dd$t0p <- dd$mtASAP2/sum(dd$mtASAP2)
dd$t1p <- dd$mtASAP3/sum(dd$mtASAP3)

observed_foldchange <- log2((dd$t0p )/(dd$t1p))


perm_df <- data.frame(clonotype = clone_df_in$clonotype,
                      experiment = sample(clone_df_in$experiment)) %>%
  group_by(clonotype,experiment) %>% summarize(count = n())

dd2 <- reshape2::dcast(perm_df, clonotype ~ experiment, value.var = "count", fill = 1)
dd2$t0p <- dd2$mtASAP2/sum(dd2$mtASAP2)
dd2$t1p <- dd2$mtASAP3/sum(dd2$mtASAP3)

permuted_foldchange <- log2((dd2$t0p )/(dd2$t1p))

density_df <- data.frame(
  logFC = c(observed_foldchange, permuted_foldchange),
  what = c(rep("observed", length(observed_foldchange)),
           rep("permuted", length(permuted_foldchange))
  ))
pDens <- ggplot(density_df, aes(x = logFC, fill = what)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values =c ("firebrick", "grey")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-4.5,4.5))+
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 FC", y = "empirical density")+
  theme(legend.position = "none")
pDens

ks.test(permuted_foldchange, observed_foldchange)

