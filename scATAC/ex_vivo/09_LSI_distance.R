# Calculation of median LSI distance as a measure of cluster
# heterogeneity (conceived by Caleb Lareau)

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

library(RANN)
library(data.table)
library(dplyr)
library(BuenColors)

set.seed(1234)
path <- "~/"

estimate_distances_nn <- function(seurat_object, k = 10, n_cells_test = 200, n_sims = 100, populations){
  test_me <- populations
  lapply(1:n_sims, function(sim_n){
    lapply(test_me, function(celltype){ 
      
      # Pull data
      set.seed(sim_n)
      cells_analyze <- sample(rownames(seurat_object@meta.data)[seurat_object@meta.data$annotation == celltype], n_cells_test)
      lsi_mat <- seurat_object@reductions$lsi@cell.embeddings[cells_analyze,2:30] # this should already be normalized but consider increasing the 30
      
      # Do LSI distance
      dist_lsi <- RANN::nn2(lsi_mat, k = k+1) # returns both the index and the dinstances of the nearest neighbors in LSI space
      per_cell_lsi_dist <- rowMeans(dist_lsi[["nn.dists"]][,2:(k+1)])
      
      data.frame(
        cell = cells_analyze,
        per_cell_lsi_dist,
        sim_n = sim_n,
        celltype,
        n = length(cells_analyze)
      )
    }) %>% rbindlist() %>% data.frame() -> all_distances_df_sim1
    all_distances_df_sim1
  }) %>% rbindlist() -> all_distances_df
  
  all_distances_df_sim <- all_distances_df %>%
    group_by(celltype, sim_n) %>%
    summarize(mean_dist = mean(per_cell_lsi_dist), median_dist = median(per_cell_lsi_dist))
  
  all_distances_df_sim
}

df_10_CMVpos4 <- estimate_distances_nn(CMVpos4_full, k = 10,
                                       populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_10_CMVpos2 <- estimate_distances_nn(CMVpos2_full, k = 10,
                                       populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_10_CMVpos3 <- estimate_distances_nn(CMVpos3_full, k = 10,
                                       populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_10_CMVpos1 <- estimate_distances_nn(CMVpos1_full, k = 10,
                                      populations = c("Adaptive1","Adaptive2",
                                                      "Adaptive3", "Adaptive4", "Adaptive5", "CD56dim"))



pdf(file = "Plots/new/Fig5_Cluster_heterogeneity.pdf", paper= "a4r")
df_10_CMVpos2 %>%
  ggplot(aes(x = celltype, y = median_dist, color = celltype)) +
  scale_color_manual(values = c(adaptivepal[1:3], clusterpal[3]))+
  geom_boxplot() + ggtitle("CMVpos2") +
  labs(x = "", y = "Mean within celltype LSI distance", color = "") +
  pretty_plot() + L_border()&NoLegend()

df_10_CMVpos4 %>%
  ggplot(aes(x = celltype, y = median_dist, color = celltype)) +
  scale_color_manual(values = c(adaptivepal[4:6], clusterpal[3]))+
  geom_boxplot() + ggtitle("CMVpos4") +
  labs(x = "", y = "Mean within celltype LSI distance", color = "") +
  pretty_plot() + L_border()&NoLegend()


df_10_CMVpos3 %>%
  ggplot(aes(x = celltype, y = median_dist, color = celltype)) +
  scale_color_manual(values = c(adaptivepal[7:9], clusterpal[3]))+
  geom_boxplot() + ggtitle("CMVpos3") +
  labs(x = "", y = "Mean within celltype LSI distance", color = "") +
  pretty_plot() + L_border()&NoLegend()

df_10_CMVpos1 %>%
  ggplot(aes(x = celltype, y = median_dist, color = celltype)) +
  scale_color_manual(values = c(adaptivepal[1:5], clusterpal[3]))+
  geom_boxplot() + ggtitle("CMVpos3") +
  labs(x = "", y = "Mean within celltype LSI distance", color = "") +
  pretty_plot() + L_border()&NoLegend()

dev.off()
# Do a few other distances for robustness
df_5 <- estimate_distances_nn(CMVpos4_full, k = 5, populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_10 <- estimate_distances_nn(CMVpos4_full, k = 10, populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_15 <- estimate_distances_nn(CMVpos4_full, k = 15, populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_20 <- estimate_distances_nn(CMVpos4_full, k = 20, populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_25 <- estimate_distances_nn(CMVpos4_full, k = 25, populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))
df_30 <- estimate_distances_nn(CMVpos4_full, k = 30, populations = c("Adaptive1","Adaptive2","Adaptive3","CD56dim"))

df_5$k <-  "05"
df_10$k <- "10"
df_15$k  <- "15"
df_20$k  <- "20"
df_25$k  <- "25"
df_30$k  <- "30"

# Combine for robustness plot
df_allk <- rbind(df_5, df_10, df_15, df_20, df_25, df_30)
df_allk %>%
  ggplot(aes(x = k, y = median_dist, color = celltype)) +
  scale_color_manual(values = c(adaptivepal[1:3], clusterpal[3]))+
  geom_boxplot(width = 0.3,   position = position_dodge(0.3) ) + ggtitle("CMVpos4") +
  labs(x = "", y = "Mean within celltype LSI distance", color = "") +
  pretty_plot() + L_border()
