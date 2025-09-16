load('rdata/extinction.integrated_data_with_metadata.Rdata')

Idents(extinction.integrated_data) <- "celltype_number"
Exc_cluster <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30","Exc_32")

p_values_list <- list()  # 每个cluster保存一个1000长度的p值向量
sampled_ext_cells <- list()
sampled_noExt_cells <- list()

for (cluster in Exc_cluster) {
  message("Processing: ", cluster)
  
  ext_cells <- WhichCells(
    extinction.integrated_data, 
    idents = cluster, 
    expression = stim == "extinction"
  )
  noExt_cells <- WhichCells(
    extinction.integrated_data, 
    idents = cluster, 
    expression = stim == "noExtinction"
  )
  
  if (length(ext_cells) < 2 || length(noExt_cells) < 2) {
    p_values_list[[cluster]] <- rep(NA, 1000)
    next
  }
  
  min_cells <- min(length(ext_cells), length(noExt_cells))
  p_values <- numeric(1000)
  
  for (seed in 1:1000) {
    set.seed(seed)
    sampled_ext <- sample(ext_cells, min_cells)
    sampled_noExt <- sample(noExt_cells, min_cells)
    
    expr_ext <- FetchData(
      extinction.integrated_data, 
      vars = "Fos", 
      cells = sampled_ext
    )$Fos
    expr_noExt <- FetchData(
      extinction.integrated_data, 
      vars = "Fos", 
      cells = sampled_noExt
    )$Fos
    
    # 两组都至少 2 个样本才做 t.test
    if (length(expr_ext) >= 2 && length(expr_noExt) >= 2) {
      p_values[seed] <- t.test(expr_ext, expr_noExt)$p.value
    } else {
      p_values[seed] <- NA
    }
    
    # 记录抽样的细胞 ID
    sampled_ext_cells[[cluster]][[seed]] <- sampled_ext
    sampled_noExt_cells[[cluster]][[seed]] <- sampled_noExt
  }
  
  p_values_list[[cluster]] <- p_values
}

save(p_values_list, file = "p_values_list_Exc_cluster_downsample_1000times_fos_expression.rdata")
save(sampled_ext_cells, file = "sampled_ext_cells_Exc_cluster_downsample_1000times_fos_expression.rdata")
save(sampled_noExt_cells, file = "sampled_noExt_cells_Exc_cluster_downsample_1000times_fos_expression.rdata")

# 示例：查看Exc_3的p值90th百分位数及对应seed
for (cluster in Exc_cluster) {
  p_90th <- quantile(p_values_list[[cluster]], 0.9, na.rm = TRUE)
  seed_90th <- which.min(abs(p_values_list[[cluster]] - p_90th))
  message(cluster, " 90th percentile p-value: ", p_90th)
  message("Corresponding seed: ", seed_90th)
}

# Exc_3 90th percentile p-value: 0.204045480170654
# Corresponding seed: 461
# Exc_4 90th percentile p-value: 0.21858592694769
# Corresponding seed: 418
# Exc_10 90th percentile p-value: 0.460293960813242
# Corresponding seed: 882
# Exc_16 90th percentile p-value: 0.839391866005146
# Corresponding seed: 56
# Exc_27 90th percentile p-value: 0.944410902332056
# Corresponding seed: 20
# Exc_28 90th percentile p-value: 0.227968386915345
# Corresponding seed: 3
# Exc_30 90th percentile p-value: 0.32107649342061
# Corresponding seed: 2
# Exc_32 90th percentile p-value: 0.747299364858321
# Corresponding seed: 45

load('rdata/extinction.integrated_data_with_metadata.Rdata')

Idents(extinction.integrated_data) <- "celltype_number"
Exc_cluster <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30","Exc_32")
Exc_cluster_seed <- c(461, 418, 882, 56, 20, 3, 2, 45)
sampled_merge_subset <- list()
for (number in seq_along(Exc_cluster)) {
  ext_cells <- WhichCells(
    extinction.integrated_data, 
    idents = Exc_cluster[number], 
    expression = stim == "extinction"
  )
  noExt_cells <- WhichCells(
    extinction.integrated_data, 
    idents = Exc_cluster[number], 
    expression = stim == "noExtinction"
  )
  
  min_cells <- min(length(ext_cells), length(noExt_cells))
  message("Processing: ", Exc_cluster[number], "_min_cells_", min_cells)
  set.seed(Exc_cluster_seed[number])
  sampled_ext <- sample(ext_cells, min_cells)
  sampled_noExt <- sample(noExt_cells, min_cells)
  sampled_merge_subset[[Exc_cluster[number]]] <- subset(extinction.integrated_data, cells = c(sampled_ext, sampled_noExt))
}
# Processing: Exc_3_min_cells_873
# Processing: Exc_4_min_cells_375
# Processing: Exc_10_min_cells_74
# Processing: Exc_16_min_cells_182
# Processing: Exc_27_min_cells_90
# Processing: Exc_28_min_cells_78
# Processing: Exc_30_min_cells_65
# Processing: Exc_32_min_cells_32

sampled_Exc_cluster_subset <- purrr::reduce(sampled_merge_subset, merge)
sampled_Exc_cluster_subset$celltype_number <-  factor(sampled_Exc_cluster_subset$celltype_number, levels = Exc_cluster)
sampled_Exc_cluster_subset$stim <-  factor(sampled_Exc_cluster_subset$stim, levels =  c("noExtinction", "extinction"))

Exc_cluster_subset <- subset(extinction.integrated_data, idents = Exc_cluster)
Exc_cluster_subset$celltype_number <-  factor(Exc_cluster_subset$celltype_number, levels = Exc_cluster)
Exc_cluster_subset$stim <-  factor(Exc_cluster_subset$stim, levels =  c("noExtinction", "extinction"))

pdf("plot/Fos_expression_level_Exc_cluster_downsample_1000times_p_val_90th_percentile.pdf",width = 6, height = 6)
DotPlot(Exc_cluster_subset,
        features = "Fos",
        group.by = "celltype_number",
        split.by = "stim") +
  ggplot2::ggtitle("raw_data")
DotPlot(sampled_Exc_cluster_subset,
        features = "Fos",
        group.by = "celltype_number",
        split.by = "stim") +
  ggplot2::ggtitle("downsample")
dev.off()

sampled_Exc_cluster_subset_Fos_expressed <- subset(sampled_Exc_cluster_subset, subset = Fos > 0)
sampled_Exc_cluster_subset_Fos_expressed$celltype_number <- factor(sampled_Exc_cluster_subset_Fos_expressed$celltype_number, levels = Exc_cluster)
sampled_Exc_cluster_subset_Fos_expressed$stim <- factor(sampled_Exc_cluster_subset_Fos_expressed$stim, levels = c("noExtinction", "extinction"))
sampled_Exc_cluster_subset_Fos_expressed$celltype_number.stim <-  factor(sampled_Exc_cluster_subset_Fos_expressed$celltype_number.stim, levels = paste0(rep(Exc_cluster,each=2), "_", rep(c("noExtinction", "extinction"),4)))

