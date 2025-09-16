#### Npy1r and Npy2r DEGs plot (vocalno plot and scatter plot) (Exc + Inh) ####
load("rdata/extinction.integrated_data_with_metadata.Rdata")
Idents(extinction.integrated_data) <- "stim"

Ext <- subset(extinction.integrated_data, idents = "extinction")
noExt <- subset(extinction.integrated_data, idents = "noExtinction")

Idents(Ext) <- "celltype"
Exc_Inh_neuron_ext <- subset(Ext, idents = c("Excitatory neuron", "Inhibitory neuron"))

Idents(noExt) <- "celltype"
Exc_Inh_neuron_noExt <- subset(noExt, idents = c("Excitatory neuron", "Inhibitory neuron"))

neuron_list <- list(Exc_Inh_neuron_ext, Exc_Inh_neuron_noExt)
neuron_name_list <- c("Exc_Inh_neuron_ext", "Exc_Inh_neuron_noExt")

Npy1r_only_neuron_list <- list()
Npy2r_only_neuron_list <- list()
diff_Npy1r_Npy2r_filter <- list()
avg_Npy1r_Npy2r <- list()

for (i in 1:2){
  Npy1r_only_neuron_list[[i]] <- subset(subset(neuron_list[[i]], subset = Npy1r > 0), subset = Npy2r <= 0)
  Npy2r_only_neuron_list[[i]] <- subset(subset(neuron_list[[i]], subset = Npy2r > 0), subset = Npy1r <= 0)
  
  avg_Npy1r <- log1p(AverageExpression(Npy1r_only_neuron_list[[i]], verbose = FALSE)$RNA)
  avg_Npy2r <- log1p(AverageExpression(Npy2r_only_neuron_list[[i]], verbose = FALSE)$RNA)
  
  avg_Npy1r_Npy2r[[i]] <- data.frame(
    gene = rownames(avg_Npy1r),
    Npy1r_all = avg_Npy1r[,1],
    Npy2r_all = avg_Npy2r[,1]
  )
  
  cells1 <- WhichCells(Npy1r_only_neuron_list[[i]])
  cells2 <- WhichCells(Npy2r_only_neuron_list[[i]])
  min_cells <- min(length(cells1), length(cells2))
  
  # 初始化存储
  diff_Npy1r_Npy2r_filter[[i]] <- vector("list", 1000)
  
  for (seed in 1:1000) {
    DEGs_filter <- FindMarkers(
      extinction.integrated_data,
      ident.1 = cells1,
      ident.2 = cells2,
      test.use = "MAST",
      max.cells.per.ident = min_cells,
      random.seed = seed,
      verbose = FALSE
    )
    diff_Npy1r_Npy2r_filter[[i]][[seed]] <- DEGs_filter  
    
    match_idx_filter <- match(rownames(avg_Npy1r_Npy2r[[i]]), rownames(DEGs_filter))
    avg_Npy1r_Npy2r[[i]][[paste0("p_val_MAST_filter_seed_", seed)]] <- DEGs_filter$p_val[match_idx_filter]
    avg_Npy1r_Npy2r[[i]][[paste0("p_adj_MAST_filter_seed_", seed)]] <- DEGs_filter$p_val_adj[match_idx_filter]
    avg_Npy1r_Npy2r[[i]][[paste0("avg_logFC_MAST_filter_seed_", seed)]] <- DEGs_filter$avg_log2FC[match_idx_filter]
    avg_Npy1r_Npy2r[[i]][[paste0("pct.1_seed_", seed)]] <- DEGs_filter$pct.1[match_idx_filter]
    avg_Npy1r_Npy2r[[i]][[paste0("pct.2_seed_", seed)]] <- DEGs_filter$pct.2[match_idx_filter]
  }
}
save(avg_Npy1r_Npy2r, file = "rdata/avg_Npy1r_Npy2r.rdata")
save(diff_Npy1r_Npy2r_filter, file = "rdata/diff_Npy1r_Npy2r_filter.rdata")

load("rdata/avg_Npy1r_Npy2r.rdata")
load("rdata/diff_Npy1r_Npy2r_filter.rdata")

# For each cluster, calculate the proportion of significant adjusted p-values (< 0.05) across 1000 seeds
for (i in seq_along(avg_Npy1r_Npy2r)) {
  df <- avg_Npy1r_Npy2r[[i]]
  
  # Identify columns for adjusted p-values
  padj_cols <- grep("^p_adj_MAST_filter_seed_", colnames(df), value = TRUE)
  
  # Calculate proportion: number of times p_adj < 0.05 divided by 1000
  # Handle NA values (e.g., if gene not detected in a run)
  df$prop_sig_padj_0.05 <- rowSums(df[, padj_cols] < 0.05, na.rm = TRUE) / 1000
  
  pval_cols <- grep("^p_val_MAST_filter_seed_", colnames(df), value = TRUE)
  df$prop_sig_pval_0.01 <- rowSums(df[, pval_cols] < 0.01, na.rm = TRUE) / 1000
  
  # Update the list
  avg_Npy1r_Npy2r[[i]] <- df
  
}
avg_Npy1r_Npy2r[[1]]$gene[avg_Npy1r_Npy2r[[1]]$prop_sig_padj_0.05 > 0.9 & abs(avg_Npy1r_Npy2r[[1]]$avg_logFC_MAST_filter_seed_1) >= 0.25]
# "Olfm1"  "Npy2r"  "Tmsb10" "Npy1r" 

avg_Npy1r_Npy2r[[2]]$gene[avg_Npy1r_Npy2r[[2]]$prop_sig_padj_0.05 > 0.9 & abs(avg_Npy1r_Npy2r[[2]]$avg_logFC_MAST_filter_seed_1) >= 0.25]
# "Npy2r" "Npy1r"

# Plot
expression_plot <- list()
for (i in 1:2){
  expression_plot[[i]] <- ggplot(data = avg_Npy1r_Npy2r[[i]], 
                                 aes(x = Npy1r_all, y = Npy2r_all, 
                                     label = gene, 
                                     colour = ((prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.25 )==T))) +
    scale_color_manual(values = c("grey50", "brown1")) +
    geom_point() +
    geom_text_repel(data=subset(avg_Npy1r_Npy2r[[i]], prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.25 ), colour='red') + 
    xlab("Npy1r +") +
    ylab("Npy2r +") + 
    xlim(0,7) +
    ylim(0,7) +
    ggtitle(neuron_name_list[i]) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1)),
          axis.text = element_text(size = rel(0.75))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) 
}

pdf("plot/Npy1r_Npy2r_expression_scatter_plot_downsample.pdf", width = 5, height = 5)
for (i in 1:2){
  print(expression_plot[[i]])
}
dev.off()

#### BH adjust ####
# 加载所需包
library(stats)  # p.adjust 来自 base R 的 stats 包

load("rdata/avg_Npy1r_Npy2r.rdata")

# 创建一个新列表来存储处理后的 data.frames（或直接修改原列表）
processed_avg_Npy1r_Npy2r <- avg_Npy1r_Npy2r  # 可以直接修改 avg_Npy1r_Npy2r

# 循环遍历每个 cluster
for (i in seq_along(processed_avg_Npy1r_Npy2r)) {
  df <- processed_avg_Npy1r_Npy2r[[i]]
  
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
  processed_avg_Npy1r_Npy2r[[i]] <- df
}

for (i in seq_along(processed_avg_Npy1r_Npy2r)) {
  df <- processed_avg_Npy1r_Npy2r[[i]]
  
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
  processed_avg_Npy1r_Npy2r[[i]] <- df
}

processed_avg_Npy1r_Npy2r[[2]]$gene[which(processed_avg_Npy1r_Npy2r[[2]]$prop_sig_padj_0.05 > 0.9 & abs(processed_avg_Npy1r_Npy2r[[2]]$avg_logFC_MAST_filter_seed_1) >= 0.585)]
# "Npy2r" "Npy1r"
processed_avg_Npy1r_Npy2r[[1]]$gene[which(processed_avg_Npy1r_Npy2r[[1]]$prop_sig_padj_0.05 > 0.9 & abs(processed_avg_Npy1r_Npy2r[[1]]$avg_logFC_MAST_filter_seed_1) >= 0.585)]
# "Rasgrp1" "Tshz2"   "Npy2r"   "Tmsb10"  "Cend1"   "Npy1r"   "Nptx1"  

# Plot
neuron_name_list <- c("Exc_Inh_neuron_ext", "Exc_Inh_neuron_noExt")

expression_plot <- list()
for (i in 1:2){
  expression_plot[[i]] <- ggplot(data = subset(processed_avg_Npy1r_Npy2r[[i]], (prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.585 )==F),
                                 aes(x = Npy1r_all, y = Npy2r_all, 
                                     label = gene)) +
    geom_point(colour ="grey50")+
    geom_point(data = subset(processed_avg_Npy1r_Npy2r[[i]], prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.585), aes(size = prop_sig_padj_0.05),
               colour = "brown1") +
    scale_size_continuous(
      range = c(1, 2),
      limits = c(0.9, 1),
      name = "Prop. Significant\n(p_adj < 0.05)"
    ) +
    geom_text_repel(data=subset(processed_avg_Npy1r_Npy2r[[i]], prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.585 ), colour='red') + 
    xlab("Npy1r +") +
    ylab("Npy2r +") + 
    xlim(0,7) +
    ylim(0,7) +
    ggtitle(neuron_name_list[i]) +
    scale_size_continuous(
      range = c(1, 2),
      limits = c(0.9, 1),
      name = "Prop. Significant\n(p_adj < 0.05)"
    )+
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1)),
          axis.text = element_text(size = rel(0.75))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) 
}

pdf("plot/Npy1r_Npy2r_expression_scatter_plot_downsample_BH_padj_log2FC_0.585.pdf", width = 5, height = 5)
for (i in 1:2){
  print(expression_plot[[i]])
}
dev.off()

expression_plot <- list()
for (i in 1:2){
  expression_plot[[i]] <- ggplot(data = subset(processed_avg_Npy1r_Npy2r[[i]], (prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.25 )==F),
                                 aes(x = Npy1r_all, y = Npy2r_all, 
                                     label = gene)) +
    geom_point(colour ="grey50")+
    geom_point(data = subset(processed_avg_Npy1r_Npy2r[[i]], prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.25), aes(size = prop_sig_padj_0.05),
               colour = "brown1") +
    scale_size_continuous(
      range = c(1, 2),
      limits = c(0.9, 1),
      name = "Prop. Significant\n(p_adj < 0.05)"
    ) +
    geom_text_repel(data=subset(processed_avg_Npy1r_Npy2r[[i]], prop_sig_padj_0.05 > 0.9 & abs(avg_logFC_MAST_filter_seed_1)>=0.25 ), colour='red') + 
    xlab("Npy1r + (capped at 5)") +
    ylab("Npy2r + (capped at 5)") + 
     scale_x_continuous(
    limits = c(0, 5),
    oob = squish   # 超出范围的点压到边界，而不是被去掉
  )  +
    scale_y_continuous(
      limits = c(0, 5),
      oob = squish   # 超出范围的点压到边界，而不是被去掉
    )  +
    ggtitle(neuron_name_list[i]) +
    scale_size_continuous(
      range = c(1, 2),
      limits = c(0.9, 1),
      name = "Prop. Significant\n(p_adj < 0.05)"
    )+
    theme(plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1)),
          axis.text = element_text(size = rel(0.75))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) 
}

pdf("plot/Npy1r_Npy2r_expression_scatter_plot_downsample_BH_padj_log2FC_0.25.pdf", width = 5, height = 5)
for (i in 1:2){
  print(expression_plot[[i]])
}
dev.off()



