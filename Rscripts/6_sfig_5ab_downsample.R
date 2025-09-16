load("rdata/extinction.integrated_data_with_metadata.Rdata")
Idents(extinction.integrated_data) <- "celltype"
celltypes <- c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage")

IEG_list <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb",
              "Jund", "Egr1", "Egr2", "Egr3", "Egr4", "Nr4a1", 
              "Nr4a2", "Nr4a3", "Arc", "Homer1", "Rheb", "Rgs2", 
              "Plk2", "Ptgs2", "Bdnf", "Inhba", "Nptx2", "Plat", 
              "Nrn1", "Myc", "Dusp1", "Dusp5", "Dusp6", "Pcdh8", 
              "Cyr61", "Gadd45b", "Trib1", "Gem", "Btg2", "Ier2", 
              "Npas4", "Rasd1", "Crem", "Mbnl2", "Arf4", "Gadd45g", 
              "Arih1", "Nup98", "Ppp1r15a", "Fbxo33", "Per1", "Per2", 
              "Maff", "Zfp36", "Srf", "Mcl1", "Ctgf", "Il6", 
              "Atf3", "Rcan1", "Ncoa7", "Cxcl2", "Bhlhe40", "Slc2a3", 
              "Nfkbia", "Ier3", "Sgk1", "Klf6", "Klf10", "Nfkbiz", 
              "Tnfaip3", "Cebpd", "Hbegf", "Ldlr", "Tsc22d1", "F3", 
              "Ccl2", "Csrnp1", "Pmaip1", "Zfp36l2", "Plau", "Ccl5", 
              "Saa3", "Tnf" , "Irf1", "Cd83", "Map3k8", "Socs3", 
              "Il1a", "Cxcl1", "Il1b", "Sod2", "Pim1", "Peli1", 
              "Tlr2", "Ccl3", "Noct", "Bcl3", "Ifit2", "Icam1", 
              "Ifit1", "Tnfsf9", "Ccrl2", "Cxcl10", "Gbp2", "Mmp13", 
              "Il23a", "Arhgef3", "Serpine1", "Traf1", "Vcam1", "Marcksl1",
              "Nfkbid", "Ikbke", "Ccl12", "Ifit3", "Cebpb", "Zfp36l1", 
              "Txnip", "Nfib", "Hes1", "Pias1", "Klf2", "Dusp2", 
              "Wee1", "Thbs1", "Sik1", "Gdf15", "Ier5", "Rgs1", 
              "Id2", "Apold1")

IEG_expr_ext <- list()
IEG_expr_noExt <- list()
for (celltype in celltypes) {
  ext_cells <- WhichCells(
    extinction.integrated_data, 
    idents = celltype, 
    expression = stim == "extinction"
  )
  
  noExt_cells <- WhichCells(
    extinction.integrated_data, 
    idents = celltype, 
    expression = stim == "noExtinction"
  )
  
  min_cells <- min(length(ext_cells), length(noExt_cells))
  
  for (seed in 1:1000) {
    set.seed(seed)
    sampled_ext <- sample(ext_cells, min_cells)
    sampled_noExt <- sample(noExt_cells, min_cells)
    
    expr_ext <- FetchData(
      extinction.integrated_data, 
      vars = IEG_list, 
      cells = sampled_ext
    )
    expr_noExt <- FetchData(
      extinction.integrated_data, 
      vars = IEG_list, 
      cells = sampled_noExt
    )
    
    IEG_expr_ext[[celltype]][[seed]] <- expr_ext
    IEG_expr_noExt[[celltype]][[seed]] <- expr_noExt
  }
}
  
save(IEG_expr_ext, file = "rdata/IEG_expr_ext_celltype_downsample_1000times.rdata")
save(IEG_expr_noExt, file = "rdata/IEG_expr_noExt_celltype_downsample_1000times.rdata")

load("rdata/extinction.integrated_data_with_metadata.Rdata")
Idents(extinction.integrated_data) <- "celltype_number"
Exc_clusters <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30","Exc_32")

IEG_list <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb",
              "Jund", "Egr1", "Egr2", "Egr3", "Egr4", "Nr4a1", 
              "Nr4a2", "Nr4a3", "Arc", "Homer1", "Rheb", "Rgs2", 
              "Plk2", "Ptgs2", "Bdnf", "Inhba", "Nptx2", "Plat", 
              "Nrn1", "Myc", "Dusp1", "Dusp5", "Dusp6", "Pcdh8", 
              "Cyr61", "Gadd45b", "Trib1", "Gem", "Btg2", "Ier2", 
              "Npas4", "Rasd1", "Crem", "Mbnl2", "Arf4", "Gadd45g", 
              "Arih1", "Nup98", "Ppp1r15a", "Fbxo33", "Per1", "Per2", 
              "Maff", "Zfp36", "Srf", "Mcl1", "Ctgf", "Il6", 
              "Atf3", "Rcan1", "Ncoa7", "Cxcl2", "Bhlhe40", "Slc2a3", 
              "Nfkbia", "Ier3", "Sgk1", "Klf6", "Klf10", "Nfkbiz", 
              "Tnfaip3", "Cebpd", "Hbegf", "Ldlr", "Tsc22d1", "F3", 
              "Ccl2", "Csrnp1", "Pmaip1", "Zfp36l2", "Plau", "Ccl5", 
              "Saa3", "Tnf" , "Irf1", "Cd83", "Map3k8", "Socs3", 
              "Il1a", "Cxcl1", "Il1b", "Sod2", "Pim1", "Peli1", 
              "Tlr2", "Ccl3", "Noct", "Bcl3", "Ifit2", "Icam1", 
              "Ifit1", "Tnfsf9", "Ccrl2", "Cxcl10", "Gbp2", "Mmp13", 
              "Il23a", "Arhgef3", "Serpine1", "Traf1", "Vcam1", "Marcksl1",
              "Nfkbid", "Ikbke", "Ccl12", "Ifit3", "Cebpb", "Zfp36l1", 
              "Txnip", "Nfib", "Hes1", "Pias1", "Klf2", "Dusp2", 
              "Wee1", "Thbs1", "Sik1", "Gdf15", "Ier5", "Rgs1", 
              "Id2", "Apold1")

IEG_expr_ext <- list()
IEG_expr_noExt <- list()
for (cluster in Exc_clusters) {
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
  
  min_cells <- min(length(ext_cells), length(noExt_cells))
  
  for (seed in 1:1000) {
    set.seed(seed)
    sampled_ext <- sample(ext_cells, min_cells)
    sampled_noExt <- sample(noExt_cells, min_cells)
    
    expr_ext <- FetchData(
      extinction.integrated_data, 
      vars = IEG_list, 
      cells = sampled_ext
    )
    expr_noExt <- FetchData(
      extinction.integrated_data, 
      vars = IEG_list, 
      cells = sampled_noExt
    )
    
    IEG_expr_ext[[cluster]][[seed]] <- expr_ext
    IEG_expr_noExt[[cluster]][[seed]] <- expr_noExt
  }
}

save(IEG_expr_ext, file = "rdata/IEG_expr_ext_cluster_downsample_1000times.rdata")
save(IEG_expr_noExt, file = "rdata/IEG_expr_noExt_cluster_downsample_1000times.rdata")

##### celltype mean log2FC ####
# 去除不在两组中均有表达的IEG
load("rdata/IEG_expr_ext_celltype_downsample_1000times.rdata")
load("rdata/IEG_expr_noExt_celltype_downsample_1000times.rdata")
IEG_expr_ext_celltype <- IEG_expr_ext
IEG_expr_noExt_celltype <- IEG_expr_noExt

IEGs_log2fc_results_celltype <- list()

for (celltype in celltypes) {
  message("Processing celltype: ", celltype)
  
  ext_list <- IEG_expr_ext_celltype[[celltype]]
  noExt_list <- IEG_expr_noExt_celltype[[celltype]]
  
  # 初始化矩阵
  ext_means_mat <- matrix(NA, nrow = length(IEG_list), ncol = 1000,
                          dimnames = list(IEG_list, 1:1000))
  noExt_means_mat <- matrix(NA, nrow = length(IEG_list), ncol = 1000,
                            dimnames = list(IEG_list, 1:1000))
  
  for (seed in 1:1000) {
    expr_ext <- ext_list[[seed]]
    expr_noExt <- noExt_list[[seed]]
    
    if (is.null(expr_ext) || is.null(expr_noExt)) next
    
    # 计算平均表达量，应用 expm1 变换，不加伪计数
    ext_means_mat[, seed] <- colMeans(expm1(expr_ext))
    noExt_means_mat[, seed] <- colMeans(expm1(expr_noExt))
  }
  
  # 每个基因 across seeds 的均值 & 标准差
  ext_mean_across_seeds <- rowMeans(ext_means_mat, na.rm = TRUE)
  noExt_mean_across_seeds <- rowMeans(noExt_means_mat, na.rm = TRUE)
  
  # 计算均值的 log2FC，如果任一组均值为 0，则设为 NA
  log2fc_mean_expr <- ifelse(ext_mean_across_seeds == 0 | noExt_mean_across_seeds == 0,
                             NA,
                             log2(ext_mean_across_seeds / noExt_mean_across_seeds))
  
  IEGs_log2fc_results_celltype[[celltype]] <- list(
    ext_means_mat = ext_means_mat,
    noExt_means_mat = noExt_means_mat,
    ext_mean_across_seeds = ext_mean_across_seeds,
    noExt_mean_across_seeds = noExt_mean_across_seeds,
    log2fc_mean_expr = log2fc_mean_expr
  )
}

log2fc_df <- bind_rows(
  lapply(names(IEGs_log2fc_results_celltype), function(ct) {
    data.frame(
      celltype = ct,
      gene = names(IEGs_log2fc_results_celltype[[ct]]$log2fc_mean_expr),
      mean_log2FC = IEGs_log2fc_results_celltype[[ct]]$log2fc_mean_expr
    )
  }),
  .id = NULL
)
log2fc_df_wide <- log2fc_df %>%
  pivot_wider(names_from = celltype, values_from = mean_log2FC) %>% 
  column_to_rownames(var = "gene") %>% na.omit()
dim(log2fc_df_wide)
# 81  9
max(log2fc_df_wide)
# 5.801238
min(log2fc_df_wide)
# -3.991529

pdf("plot/IEG_log2FC_heatmap_celltype_downsample_1000times_average.pdf", height = 8, width = 5)
my_palette <- c(
  colorRampPalette(c("blue4", "white"))(40),
  "white",
  colorRampPalette(c("white", "red"))(58)
)

# 颜色映射范围
breaks <- seq(-2, 3, length.out = length(my_palette))

pheatmap(log2fc_df_wide,
         cluster_col = FALSE, 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_palette,
         breaks = breaks,                         # 固定颜色范围
         border_color = NA,
         legend = TRUE,
         legend_breaks = seq(-2, 3, by = 1),  # 修改 legend 范围
         fontsize_row = 4,
         fontsize_col = 8,
         fontsize = 9)
dev.off()

##### cluster mean log2FC ####
load("rdata/IEG_expr_ext_cluster_downsample_1000times.rdata")
load("rdata/IEG_expr_noExt_cluster_downsample_1000times.rdata")
IEG_expr_ext_cluster <- IEG_expr_ext
IEG_expr_noExt_cluster <- IEG_expr_noExt

IEGs_log2fc_results_cluster <- list()

for (cluster in Exc_clusters) {
  message("Processing cluster: ", cluster)
  
  ext_list <- IEG_expr_ext_cluster[[cluster]]
  noExt_list <- IEG_expr_noExt_cluster[[cluster]]
  
  # 初始化矩阵
  ext_means_mat <- matrix(NA, nrow = length(IEG_list), ncol = 1000,
                          dimnames = list(IEG_list, 1:1000))
  noExt_means_mat <- matrix(NA, nrow = length(IEG_list), ncol = 1000,
                            dimnames = list(IEG_list, 1:1000))
  
  for (seed in 1:1000) {
    expr_ext <- ext_list[[seed]]
    expr_noExt <- noExt_list[[seed]]
    
    if (is.null(expr_ext) || is.null(expr_noExt)) next
    
    # 计算平均表达量，应用 expm1 变换，不加伪计数
    ext_means_mat[, seed] <- colMeans(expm1(expr_ext))
    noExt_means_mat[, seed] <- colMeans(expm1(expr_noExt))
  }
  
  # 每个基因 across seeds 的均值 & 标准差
  ext_mean_across_seeds <- rowMeans(ext_means_mat, na.rm = TRUE)
  noExt_mean_across_seeds <- rowMeans(noExt_means_mat, na.rm = TRUE)
  
  # 计算均值的 log2FC，如果任一组均值为 0，则设为 NA
  log2fc_mean_expr <- ifelse(ext_mean_across_seeds == 0 | noExt_mean_across_seeds == 0,
                             NA,
                             log2(ext_mean_across_seeds / noExt_mean_across_seeds))
  
  IEGs_log2fc_results_cluster[[cluster]] <- list(
    ext_means_mat = ext_means_mat,
    noExt_means_mat = noExt_means_mat,
    ext_mean_across_seeds = ext_mean_across_seeds,
    noExt_mean_across_seeds = noExt_mean_across_seeds,
    log2fc_mean_expr = log2fc_mean_expr
  )
}

log2fc_df <- bind_rows(
  lapply(names(IEGs_log2fc_results_cluster), function(ct) {
    data.frame(
      cluster = ct,
      gene = names(IEGs_log2fc_results_cluster[[ct]]$log2fc_mean_expr),
      mean_log2FC = IEGs_log2fc_results_cluster[[ct]]$log2fc_mean_expr
    )
  }),
  .id = NULL
)
log2fc_df_wide <- log2fc_df %>%
  pivot_wider(names_from = cluster, values_from = mean_log2FC) %>% 
  column_to_rownames(var = "gene") %>% na.omit()
dim(log2fc_df_wide)
# 57  8
max(log2fc_df_wide)
# 2.907987
min(log2fc_df_wide)
# -2.529752

pdf("plot/IEG_log2FC_heatmap_cluster_downsample_1000times_average.pdf", height = 8, width = 5)
my_palette <- c(colorRampPalette(c("blue4", "white"))(50), "white",colorRampPalette(c("white", "red"))(58))
pheatmap(log2fc_df_wide,
         cluster_col = FALSE, 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_palette,
         border_color = NA,
         legend = TRUE,
         legend_breaks = seq(-2,2),
         fontsize_row = 4,      # 行名字号（比如基因名）
         fontsize_col = 8,      # 列名字号（比如细胞类型）
         fontsize = 9)           # 其余整体字体
dev.off()


