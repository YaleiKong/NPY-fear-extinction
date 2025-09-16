##### celltype (MAST + downsample + filter) (seed 1-1000) ####
# logfc.threshold = 0.25,
# min.pct = 0.1,
library(future)

# Set to multisession with appropriate workers (adjust based on your system)
plan("multisession", workers = 39)  # Example: 4 workers; change to your available cores -1

# Load data
load('rdata/extinction.integrated_data_with_metadata.Rdata')

# Start timing
start_time <- proc.time()

celltype <- c("Excitatory neuron", "Inhibitory neuron", "Astrocyte", "Oligodendrocyte", 
              "OPC", "Microglia", "Pericyte", "Endothelium, smooth muscle", "Macrophage")

# 初始化列表
celltype_subset <- vector("list", length(celltype))
avg_celltype <- vector("list", length(celltype))
DEGs_MAST_filter_runs <- vector("list", length(celltype))  # List of lists for 1000 runs per celltype

for (i in seq_along(celltype)) {
  # 子集提取
  Idents(extinction.integrated_data) <- "celltype"
  celltype_subset[[i]] <- subset(extinction.integrated_data, idents = celltype[i])
  
  # 平均表达 (不变，无需多次运行)
  Idents(celltype_subset[[i]]) <- "stim"
  avg_expr <- log1p(AverageExpression(celltype_subset[[i]], verbose = FALSE)$RNA) # log(1 + x)
  avg_celltype[[i]] <- avg_expr <- as.data.frame(avg_expr)
  avg_expr$gene <- rownames(avg_expr)
  
  # MAST分析：运行1000次，种子1-1000
  Idents(extinction.integrated_data) <- "celltype.stim"
  celltype_noExtinction <- paste0(celltype[i], "_noExtinction")
  celltype_extinction <- paste0(celltype[i], "_extinction")
  cells1 <- WhichCells(extinction.integrated_data, idents = celltype_extinction)
  cells2 <- WhichCells(extinction.integrated_data, idents = celltype_noExtinction)
  min_cells <- min(length(cells1), length(cells2))
  
  # 初始化当前 celltype 的 1000 次运行列表
  DEGs_MAST_filter_runs[[i]] <- vector("list", 1000)
  
  for (seed in 1:1000) {
    # MAST without filter, with seed 1-1000
    DEGs_filter <- FindMarkers(extinction.integrated_data, 
                               ident.1 = celltype_extinction,
                               ident.2 = celltype_noExtinction,
                               test.use = "MAST",
                               max.cells.per.ident = min_cells,
                               random.seed = seed,  # Seed from 1 to 1000
                               verbose = FALSE)
    DEGs_MAST_filter_runs[[i]][[seed]] <- DEGs_filter  # Save to sublist
    
    # 匹配并添加到 avg_expr (可选: 如果需要平均多个运行的结果，这里只处理一个；否则平均 p_val 等)
    match_idx_filter <- match(rownames(avg_expr), rownames(DEGs_filter))
    avg_expr[[paste0("p_val_MAST_filter_seed_", seed)]] <- DEGs_filter$p_val[match_idx_filter]
    avg_expr[[paste0("p_adj_MAST_filter_seed_", seed)]] <- DEGs_filter$p_val_adj[match_idx_filter]
    avg_expr[[paste0("avg_logFC_MAST_filter_seed_", seed)]] <- DEGs_filter$avg_log2FC[match_idx_filter]
    avg_expr[[paste0("pct.1_seed_", seed)]] <- DEGs_filter$pct.1[match_idx_filter]
    avg_expr[[paste0("pct.2_seed_", seed)]] <- DEGs_filter$pct.2[match_idx_filter]
  }
  
  # 保存 avg_expr (现在包含所有种子的列)
  avg_celltype[[i]] <- avg_expr
}

# End timing
end_time <- proc.time()

# Calculate elapsed time
elapsed_time <- end_time - start_time
print(paste("Elapsed time:", elapsed_time[3], "seconds"))

# Reset to sequential
plan("sequential")

save(avg_celltype,file = "rdata/DE_results_avg_celltype_downsample_1000times.rdata")
save(DEGs_MAST_filter_runs,file = "rdata/DEGs_MAST_filter_celltype_downsample_1000times.rdata")

#### BH adjust ####
# 加载所需包
library(stats)  # p.adjust 来自 base R 的 stats 包

load("rdata/DE_results_avg_celltype_downsample_1000times.rdata")

# 创建一个新列表来存储处理后的 data.frames（或直接修改原列表）
processed_celltypes <- avg_celltype_all  # 可以直接修改 avg_celltype_all

# 循环遍历每个 celltype
for (i in seq_along(processed_celltypes)) {
  df <- processed_celltypes[[i]]
  
  # 识别 p_val 列（假设有 1000 个 seed 的 p_val 列）
  pval_cols <- grep("^p_val_MAST_filter_seed_", colnames(df), value = TRUE)
  
  # 确保有正好 1000 个这样的列
  if (length(pval_cols) != 1000) {
    warning(paste("celltype", i, "has", length(pval_cols), "p_val columns instead of 1000."))
  }
  
  # 对于每个 p_val 列，进行 BH 校正并创建对应的 p_adj 列
  for (pval_col in pval_cols) {
    # 提取当前 seed 的 p 值向量
    p_vals <- df[[pval_col]]
    
    # 识别非 NA 的索引
    non_na_idx <- !is.na(p_vals)
    
    # 提取非 NA 的 p 值
    p_vals_non_na <- p_vals[non_na_idx]
    
    # 如果向量为空，创建全 NA 的 p_adj 列并跳过
    if (length(p_vals_non_na) == 0) {
      padj_col_name <- sub("^p_val_", "p_val_adj_", pval_col)  # 生成对应的 p_adj 列名，如 p_adj_MAST_filter_seed_1
      df[[padj_col_name]] <- NA
      next
    }
    
    # 使用 BH 方法校正 p 值
    # 默认：n = length(p_vals_non_na)，基于当前测试数（过滤后的基因数）
    # p_vals_adj <- p.adjust(p_vals_non_na, method = "BH")
    
    # 可选：如果想用固定总基因数作为 n（例如 Seurat 风格，更保守）
    total_genes <- nrow(extinction.integrated_data) 
    p_vals_adj <- p.adjust(p_vals_non_na, method = "BH", n = total_genes)
    
    # 创建对应的 p_adj 列名
    padj_col_name <- sub("^p_val_", "p_val_adj_", pval_col)
    
    # 初始化 p_adj 列为 NA，然后填充非 NA 位置
    df[[padj_col_name]] <- NA
    df[[padj_col_name]][non_na_idx] <- p_vals_adj
  }
  
  # 更新列表
  processed_celltypes[[i]] <- df
}

# 查看结果（例如第一个 celltype 的前几行和列）
head(processed_celltypes[[1]][, c(pval_cols[1], sub("^p_val_", "p_val_adj_", pval_cols[1]))])

sum(processed_celltypes[[2]]["Npy",grep("p_val_adj_",colnames(processed_celltypes[[2]]))]<0.05)
# 947
save(processed_celltypes, file = "rdata/DE_results_avg_celltype_downsample_1000times_BH_adjust.rdata")

for (i in seq_along(processed_celltypes)) {
  df <- processed_celltypes[[i]]
  
  # Identify columns for adjusted p-values
  padj_cols <- grep("^p_val_adj_MAST_filter_seed_", colnames(df), value = TRUE)
  
  # Ensure there are exactly 1000 such columns
  if (length(padj_cols) != 1000) {
    warning(paste("celltype", i, "has", length(padj_cols), "p_adj columns instead of 1000."))
  }
  
  # Calculate proportion: number of times p_adj < 0.05 divided by 1000
  # Handle NA values (e.g., if gene not detected in a run)
  df$prop_sig_padj_0.05 <- rowSums(df[, padj_cols] < 0.05, na.rm = TRUE) / 1000
  
  # Update the list
  processed_celltypes[[i]] <- df
}
processed_celltypes[[2]]$gene[which(processed_celltypes[[2]]$prop_sig_padj_0.05 > 0.9)]
save(processed_celltypes, file = "rdata/DE_results_avg_celltype_downsample_1000times_BH_adjust_with_prop.rdata")

sum(processed_celltypes[[2]]$prop_sig_padj_0.05 > 0.9)
sum(processed_celltypes[[2]]$prop_sig_padj_0.05 > 0.8)
# 144

##### cluster (MAST + downsample + filter) (seed 1-1000) ####
# logfc.threshold = 0.25,
# min.pct = 0.1,

library(future)

# Set to multisession with appropriate workers (adjust based on your system)
plan("multisession", workers = 39)  # Example: 4 workers; change to your available cores -1

# Load data
load('rdata/extinction.integrated_data_with_metadata.Rdata')

# Start timing
start_time <- proc.time()

Idents(extinction.integrated_data) <- "celltype_number"
cluster <- unique(extinction.integrated_data@meta.data$celltype_number) %>% 
  as.character()
# [1] "Endo_9"        "Exc_27"        "Exc_3"         "Astro_5"       "Endo_2"        "OPC_8"        
# [7] "Micro_13"      "Pericyte_17"   "Olig_20"       "Pericyte_7"    "Pericyte_18"   "Olig_0"       
# [13] "Micro_29"      "Astro_1"       "Exc_16"        "Endo_21"       "Macrophage_15" "Micro_25"     
# [19] "Inh_6"         "Micro_12"      "Pericyte_23"   "Astro_26"      "Exc_4"         "Olig_14"      
# [25] "Endo_11"       "OPC_19"        "Exc_30"        "Exc_10"        "Olig_22"       "Exc_28"       
# [31] "Micro_24"      "Olig_33"       "Astro_31"      "Astro_34"      "Exc_32"
# 初始化列表
cluster_subset <- vector("list", length(cluster))
avg_cluster <- vector("list", length(cluster))
DEGs_MAST_filter_runs <- vector("list", length(cluster))  # List of lists for 1000 runs per cluster

for (i in seq_along(cluster)) {
  # 子集提取
  Idents(extinction.integrated_data) <- "celltype_number"
  cluster_subset[[i]] <- subset(extinction.integrated_data, idents = cluster[i])
  
  # 平均表达 (不变，无需多次运行)
  Idents(cluster_subset[[i]]) <- "stim"
  avg_expr <- log1p(AverageExpression(cluster_subset[[i]], verbose = FALSE)$RNA) # log(1 + x)
  avg_cluster[[i]] <- avg_expr <- as.data.frame(avg_expr)
  avg_expr$gene <- rownames(avg_expr)
  
  # MAST分析：运行1000次，种子1-1000
  Idents(extinction.integrated_data) <- "celltype_number.stim"
  cluster_noExtinction <- paste0(cluster[i], "_noExtinction")
  cluster_extinction <- paste0(cluster[i], "_extinction")
  cells1 <- WhichCells(extinction.integrated_data, idents = cluster_extinction)
  cells2 <- WhichCells(extinction.integrated_data, idents = cluster_noExtinction)
  min_cells <- min(length(cells1), length(cells2))
  
  # 初始化当前 cluster 的 1000 次运行列表
  DEGs_MAST_filter_runs[[i]] <- vector("list", 1000)
  
  for (seed in 1:1000) {
    # MAST without filter, with seed 1-1000
    DEGs_filter <- FindMarkers(extinction.integrated_data, 
                               ident.1 = cluster_extinction,
                               ident.2 = cluster_noExtinction,
                               test.use = "MAST",
                               max.cells.per.ident = min_cells,
                               random.seed = seed,  # Seed from 1 to 1000
                               verbose = FALSE)
    DEGs_MAST_filter_runs[[i]][[seed]] <- DEGs_filter  # Save to sublist
    
    # 匹配并添加到 avg_expr (可选: 如果需要平均多个运行的结果，这里只处理一个；否则平均 p_val 等)
    match_idx_filter <- match(rownames(avg_expr), rownames(DEGs_filter))
    avg_expr[[paste0("p_val_MAST_filter_seed_", seed)]] <- DEGs_filter$p_val[match_idx_filter]
    avg_expr[[paste0("p_adj_MAST_filter_seed_", seed)]] <- DEGs_filter$p_val_adj[match_idx_filter]
    avg_expr[[paste0("avg_logFC_MAST_filter_seed_", seed)]] <- DEGs_filter$avg_log2FC[match_idx_filter]
    avg_expr[[paste0("pct.1_seed_", seed)]] <- DEGs_filter$pct.1[match_idx_filter]
    avg_expr[[paste0("pct.2_seed_", seed)]] <- DEGs_filter$pct.2[match_idx_filter]
  }
  
  # 保存 avg_expr (现在包含所有种子的列)
  avg_cluster[[i]] <- avg_expr
}

# End timing
end_time <- proc.time()

# Calculate elapsed time
elapsed_time <- end_time - start_time
print(paste("Elapsed time:", elapsed_time[3], "seconds"))

# Reset to sequential
plan("sequential")

save(avg_cluster,file = "rdata/DE_results_avg_cluster_downsample_1000times.rdata")
save(DEGs_MAST_filter_runs,file = "rdata/DEGs_MAST_filter_cluster_downsample_1000times.rdata")

#### BH adjust ####
# 加载所需包
library(stats)  # p.adjust 来自 base R 的 stats 包

load("rdata/DE_results_avg_cluster_downsample_1000times.rdata")

# 创建一个新列表来存储处理后的 data.frames（或直接修改原列表）
processed_clusters <- avg_cluster_all  # 可以直接修改 avg_cluster_all

# 循环遍历每个 cluster
for (i in seq_along(processed_clusters)) {
  df <- processed_clusters[[i]]
  
  # 识别 p_val 列（假设有 1000 个 seed 的 p_val 列）
  pval_cols <- grep("^p_val_MAST_filter_seed_", colnames(df), value = TRUE)
  
  # 确保有正好 1000 个这样的列
  if (length(pval_cols) != 1000) {
    warning(paste("cluster", i, "has", length(pval_cols), "p_val columns instead of 1000."))
  }
  
  # 对于每个 p_val 列，进行 BH 校正并创建对应的 p_adj 列
  for (pval_col in pval_cols) {
    # 提取当前 seed 的 p 值向量
    p_vals <- df[[pval_col]]
    
    # 识别非 NA 的索引
    non_na_idx <- !is.na(p_vals)
    
    # 提取非 NA 的 p 值
    p_vals_non_na <- p_vals[non_na_idx]
    
    # 如果向量为空，创建全 NA 的 p_adj 列并跳过
    if (length(p_vals_non_na) == 0) {
      padj_col_name <- sub("^p_val_", "p_val_adj_", pval_col)  # 生成对应的 p_adj 列名，如 p_adj_MAST_filter_seed_1
      df[[padj_col_name]] <- NA
      next
    }
    
    # 使用 BH 方法校正 p 值
    # 默认：n = length(p_vals_non_na)，基于当前测试数（过滤后的基因数）
    # p_vals_adj <- p.adjust(p_vals_non_na, method = "BH")
    
    # 可选：如果想用固定总基因数作为 n（例如 Seurat 风格，更保守）
    total_genes <- nrow(extinction.integrated_data)  # 或 dim(seurat_object@assays$RNA)[1]，取决于上下文
    p_vals_adj <- p.adjust(p_vals_non_na, method = "BH", n = total_genes)
    
    # 创建对应的 p_adj 列名
    padj_col_name <- sub("^p_val_", "p_val_adj_", pval_col)
    
    # 初始化 p_adj 列为 NA，然后填充非 NA 位置
    df[[padj_col_name]] <- NA
    df[[padj_col_name]][non_na_idx] <- p_vals_adj
  }
  
  # 更新列表
  processed_clusters[[i]] <- df
}

sum(processed_clusters[[19]]["Npy",grep("p_val_adj_",colnames(processed_clusters[[19]]))]<0.05)
save(processed_clusters, file = "rdata/DE_results_avg_cluster_downsample_1000times_BH_adjust.rdata")

# For each cluster, calculate the proportion of significant adjusted p-values (< 0.05) across 1000 seeds
for (i in seq_along(processed_clusters)) {
  df <- processed_clusters[[i]]
  
  # Identify columns for adjusted p-values
  padj_cols <- grep("^p_val_adj_MAST_filter_seed_", colnames(df), value = TRUE)
  
  # Calculate proportion: number of times p_adj < 0.05 divided by 1000
  # Handle NA values (e.g., if gene not detected in a run)
  df$prop_sig_padj_0.05 <- rowSums(df[, padj_cols] < 0.05, na.rm = TRUE) / 1000
  
  # Update the list
  processed_clusters[[i]] <- df
  
}
save(processed_clusters, file = "rdata/DE_results_avg_cluster_downsample_1000times_BH_adjust_with_prop.rdata")

#### DE volcano plot ####
###### celltype ####
load("rdata/avg_celltype.rdata")
avg_celltype
load("rdata/DE_results_avg_celltype_downsample_1000times_BH_adjust_with_prop.rdata")
processed_celltypes
all(rownames(avg_celltype[[2]]) == rownames(processed_celltypes[[2]]))

# Loop through indices 1 to 9
for (i in seq_along(processed_celltypes)) {
  # Extract data frames
  proc_df <- processed_celltypes[[i]]
  avg_df <- avg_celltype[[i]]
  
  # Check if prop_sig_padj_0.05 exists in processed_celltypes
  if (!"prop_sig_padj_0.05" %in% colnames(proc_df)) {
    warning(paste("prop_sig_padj_0.05 not found in processed_celltypes[[", i, "]]"))
    next
  }
  
  # Verify matching genes (assuming rownames are genes)
  if (!identical(rownames(proc_df), rownames(avg_df))) {
    # If gene column exists instead of rownames
    if ("gene" %in% colnames(proc_df) && "gene" %in% colnames(avg_df)) {
      if (!identical(proc_df$gene, avg_df$gene)) {
        warning(paste("Gene mismatch in celltype", i))
        next
      }
    } else {
      warning(paste("Rownames mismatch and no gene column in celltype", i))
      next
    }
  }
  
  # Add prop_sig_padj_0.05 to avg_celltype
  avg_df$prop_sig_padj_0.05 <- proc_df$prop_sig_padj_0.05
  
  # Update the list
  avg_celltype[[i]] <- avg_df
}

avg_celltype[[2]]["Npy",]

celltype_volcanoPlot_list <- list()
for (i in seq_along(celltypes)) {
  
  if (!all(required_cols %in% colnames(avg_celltype[[i]]))) {
    warning(sprintf("Missing required columns in avg_celltype[[%d]]", i))
    next
  }
  
  plot_data <- avg_celltype[[i]] %>%
    filter(!is.na(avg_logFC_MAST_no_filter),
           !is.na(p_adj_MAST_no_filter),
           !is.na(prop_sig_padj_0.05),
           is.finite(avg_logFC_MAST_no_filter),
           p_adj_MAST_no_filter > 0) %>%
    mutate(neg_log10_p = -log10(p_adj_MAST_no_filter))
  
  if (nrow(plot_data) == 0) {
    warning(sprintf("No valid data for celltype %s", celltypes[i]))
    next
  }
  
  sig_subset <- filter(plot_data, prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST) > 0.585)
  cat(sprintf("Celltype %s has %d significant genes\n", celltypes[i], nrow(sig_subset)))
  
  p <- ggplot() +
    # 非显著点
    geom_point(
      data = filter(plot_data, !(prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST) > 0.585)),
      aes(x = avg_logFC_MAST_no_filter, y = neg_log10_p),
      colour = "grey50", size = 1
    ) +
    # 显著点（带 size legend）
    geom_point(
      data = sig_subset,
      aes(
        x = avg_logFC_MAST_no_filter, 
        y = neg_log10_p,
        size = prop_sig_padj_0.05
      ),
      colour = "brown1"
    ) +
    scale_size_continuous(
      range = c(1, 2),
      limits = c(0.9, 1),
      name = "Prop. Significant\n(p_adj < 0.05)"
    ) +
    xlab("log2FC") + ylab("-log10(FDR)") +
    ggtitle(celltypes[i]) +
    xlim(-3, 3) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    ) +
    guides(size = guide_legend(
      override.aes = list(colour = "red")  # legend bar 改为红色
    ))
  
  if (nrow(sig_subset) > 0) {
    p <- p +
      geom_text_repel(
        data = sig_subset,
        aes(
          x = avg_logFC_MAST_no_filter,
          y = neg_log10_p,
          label = gene
        ),
        colour = "red",
        size = 3
      )
  }
  
  celltype_volcanoPlot_list[[i]] <- p
}

dir.create("volcano_plots", showWarnings = FALSE)
for (i in seq_along(celltypes)) {
  if (!is.null(celltype_volcanoPlot_list[[i]])) {
    ggsave(
      filename = sprintf("volcano_plots/volcano_plot_%s.pdf", gsub("[ ,]", "_", celltypes[i])),
      plot = celltype_volcanoPlot_list[[i]],
      width = 7, height = 6, dpi = 300
    )
  }
}

###### cluster ####
load("rdata/avg_celltype_number.rdata")
celltype_number <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                     "Inh_6", 
                     "Astro_1","Astro_5","Astro_26","Astro_31", 'Astro_34',
                     "Olig_0", "Olig_14","Olig_20","Olig_22","Olig_33",
                     "OPC_8", "OPC_19", 
                     "Micro_12","Micro_13","Micro_24","Micro_25","Micro_29", 
                     "Pericyte_7","Pericyte_17","Pericyte_18","Pericyte_23",
                     "Endo_2","Endo_9","Endo_11","Endo_21", 
                     "Macrophage_15")
names(avg_celltype_number) <-celltype_number
clusters <- unique(extinction.integrated_data@meta.data$celltype_number) %>% 
  as.character()
load("rdata/DE_results_avg_cluster_downsample_1000times_BH_adjust_with_prop.rdata")
names(processed_clusters) <- clusters
processed_clusters_reordered <- processed_clusters[celltype_number]
all(rownames(avg_celltype_number[[2]]) == rownames(processed_clusters_reordered[[2]]))

# Loop through indices 1 to 9
for (i in seq_along(processed_clusters_reordered)) {
  # Extract data frames
  proc_df <- processed_clusters_reordered[[i]]
  avg_df <- avg_celltype_number[[i]]
  
  # Check if prop_sig_padj_0.05 exists in processed_clusters_reordered
  if (!"prop_sig_padj_0.05" %in% colnames(proc_df)) {
    warning(paste("prop_sig_padj_0.05 not found in processed_clusters_reordered[[", i, "]]"))
    next
  }
  
  # Verify matching genes (assuming rownames are genes)
  if (!identical(rownames(proc_df), rownames(avg_df))) {
    # If gene column exists instead of rownames
    if ("gene" %in% colnames(proc_df) && "gene" %in% colnames(avg_df)) {
      if (!identical(proc_df$gene, avg_df$gene)) {
        warning(paste("Gene mismatch in cluster", i))
        next
      }
    } else {
      warning(paste("Rownames mismatch and no gene column in cluster", i))
      next
    }
  }
  
  # Add prop_sig_padj_0.05 to avg_celltype_number
  avg_df$prop_sig_padj_0.05 <- proc_df$prop_sig_padj_0.05
  
  # Update the list
  avg_celltype_number[[i]] <- avg_df
}

avg_celltype_number[[9]]["Npy",]

cluster_volcanoPlot_list <- list()
for (i in seq_along(celltype_number)) {
  
  if (!all(required_cols %in% colnames(avg_celltype_number[[i]]))) {
    warning(sprintf("Missing required columns in avg_celltype_number[[%d]]", i))
    next
  }
  
  plot_data <- avg_celltype_number[[i]] %>%
    filter(!is.na(avg_logFC_MAST_no_filter),
           !is.na(p_adj_MAST_no_filter),
           !is.na(prop_sig_padj_0.05),
           is.finite(avg_logFC_MAST_no_filter),
           p_adj_MAST_no_filter > 0) %>%
    mutate(neg_log10_p = -log10(p_adj_MAST_no_filter))
  
  if (nrow(plot_data) == 0) {
    warning(sprintf("No valid data for cluster %s", celltype_number[i]))
    next
  }
  
  sig_subset <- filter(plot_data, prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST) > 0.585)
  cat(sprintf("cluster %s has %d significant genes\n", celltype_number[i], nrow(sig_subset)))
  
  p <- ggplot() +
    # 非显著点
    geom_point(
      data = filter(plot_data, !(prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST) > 0.585)),
      aes(x = avg_logFC_MAST_no_filter, y = neg_log10_p),
      colour = "grey50", size = 1
    ) +
    # 显著点（带 size legend）
    geom_point(
      data = sig_subset,
      aes(
        x = avg_logFC_MAST_no_filter, 
        y = neg_log10_p,
        size = prop_sig_padj_0.05
      ),
      colour = "brown1"
    ) +
    scale_size_continuous(
      range = c(1, 2),
      limits = c(0.9, 1),
      name = "Prop. Significant\n(p_adj < 0.05)"
    ) +
    xlab("log2FC") + ylab("-log10(FDR)") +
    ggtitle(celltype_number[i]) +
    xlim(-3, 3) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    ) +
    guides(size = guide_legend(
      override.aes = list(colour = "red")  # legend bar 改为红色
    ))
  
  if (nrow(sig_subset) > 0) {
    p <- p +
      geom_text_repel(
        data = sig_subset,
        aes(
          x = avg_logFC_MAST_no_filter,
          y = neg_log10_p,
          label = gene
        ),
        colour = "red",
        size = 3
      )
  }
  
  cluster_volcanoPlot_list[[i]] <- p
}

dir.create("volcano_plots", showWarnings = FALSE)
for (i in seq_along(celltype_number)) {
  if (!is.null(cluster_volcanoPlot_list[[i]])) {
    ggsave(
      filename = sprintf("volcano_plots/volcano_plot_%s.pdf", gsub("[ ,]", "_", celltype_number[i])),
      plot = cluster_volcanoPlot_list[[i]],
      width = 7, height = 6, dpi = 300
    )
  }
}