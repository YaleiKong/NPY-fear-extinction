##### Npy expression (IEG+ Inh neuron)  (at least one IEG⁺ in the 90th percentile among all inhibitory neurons) #####
#### seed (p-value 90th) ####
for (seed in 1:1000) {
  ext_subset <- Inh_ext_downsample_subset_IEG_expressed_list[[seed]]
  noExt_subset <- Inh_noExt_subset_IEG_expressed_list[[seed]]
  # 获取Npy表达数据
  ext_expr <- FetchData(ext_subset, vars = "Npy")$Npy
  noExt_expr <- FetchData(noExt_subset, vars = "Npy")$Npy
  
  # t.test
  t_res <- t.test(ext_expr, noExt_expr)
  p_values[seed] <- t_res$p.value
}
p_90th <- quantile(p_values, 0.9, na.rm = TRUE)
seed_90th <- which.min(abs(p_values - p_90th))
# 722
p_values[722]
# 0.004703304
#### seed 722 ####
Idents(extinction.integrated_data) <- "celltype.stim"
Inh_noExt_cells <- WhichCells(extinction.integrated_data, idents = "Inhibitory neuron_noExtinction")
Inh_ext_cells <- WhichCells(extinction.integrated_data, idents = "Inhibitory neuron_extinction")
min_Inh_cells <- min(length(Inh_noExt_cells), length(Inh_ext_cells))
Inh_ext_subset <- subset(extinction.integrated_data, idents = "Inhibitory neuron_extinction")
Inh_noExt_subset <- subset(extinction.integrated_data, idents = "Inhibitory neuron_noExtinction")

set.seed(722)
Inh_ext_cells_downsample_seed_722 <- sample(Inh_ext_cells, min_Inh_cells)
# 抽取子集
Inh_ext_downsample_subset_seed_722 <- subset(
  Inh_ext_subset, 
  cells = Inh_ext_cells_downsample_seed_722
)

Inh_downsample_subset_seed_722 <- merge(Inh_ext_downsample_subset_seed_722, Inh_noExt_subset)


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

###### Npy + percent (ext vs. noExt) (fig s7) #####
Inh_downsample_subset_seed_722_data <- GetAssayData(Inh_downsample_subset_seed_722, assay = "RNA", slot = "data")
# slot = c("data", "scale.data", "counts")

# Replace zeros with NA to avoid including them in percentile calculations
IEGs_expression_level_90th_percentile_Inh_downsample_subset_seed_722 <- Inh_downsample_subset_seed_722_data[IEG_list, ]
IEGs_expression_level_90th_percentile_Inh_downsample_subset_seed_722[IEGs_expression_level_90th_percentile_Inh_downsample_subset_seed_722 == 0] <- NA

# Calculate the 90th percentile for each gene (row) in the dataset
IEGs_expression_level_90th_percentile_Inh_downsample_subset_seed_722 <- apply(IEGs_expression_level_90th_percentile_Inh_downsample_subset_seed_722, 
                                                                     1, 
                                                                     function(x) quantile(x, probs = 0.90, na.rm = TRUE))

IEGs_data_Inh_downsample_subset_seed_722 <- Inh_downsample_subset_seed_722_data[IEG_list,]
IEGs_data_Inh_downsample_subset_seed_722[which(IEGs_data_Inh_downsample_subset_seed_722 <= IEGs_expression_level_90th_percentile_Inh_downsample_subset_seed_722)] <- NA
IEG_number_Inh_downsample_subset_seed_722 <- apply(IEGs_data_Inh_downsample_subset_seed_722, 2, function(x) sum(x > 0, na.rm = T))
IEG_positive_cells_Inh_downsample_subset_seed_722 <- names(IEG_number_Inh_downsample_subset_seed_722[IEG_number_Inh_downsample_subset_seed_722>0])
length(IEG_positive_cells_Inh_downsample_subset_seed_722)
# 305
ncol(IEGs_data_Inh_downsample_subset_seed_722)
# 556

IEGs_positive_Inh_downsample_subset_seed_722 <- extinction.integrated_data[,IEG_positive_cells_Inh_downsample_subset_seed_722]
Idents(IEGs_positive_Inh_downsample_subset_seed_722) <- "stim"
IEGs_positive_Inh_downsample_subset_seed_722_ext <- subset(IEGs_positive_Inh_downsample_subset_seed_722, idents = "extinction")
ncol(IEGs_positive_Inh_downsample_subset_seed_722_ext)
# 171
Idents(Inh_downsample_subset_seed_722) <- "stim"
ncol(subset(Inh_downsample_subset_seed_722, idents = "extinction"))
# 278
171/278
# 0.6151079
IEGs_positive_Inh_downsample_subset_seed_722_noExt <- subset(IEGs_positive_Inh_downsample_subset_seed_722, idents = "noExtinction")
ncol(IEGs_positive_Inh_downsample_subset_seed_722_noExt)
# 134
ncol(subset(Inh_downsample_subset_seed_722, idents = "noExtinction"))
# 278
134/278
# 0.4820144
ncol(subset(IEGs_positive_Inh_downsample_subset_seed_722_ext, subset = Npy > 0))
# 41
41/171
# 0.2397661
ncol(subset(IEGs_positive_Inh_downsample_subset_seed_722_noExt, subset = Npy > 0))
# 20
20/134
# 0.1492537

# 定义计数值
npy_counts <- list(
  extinction = c("Npy+" = 41, "Npy-" = 171 - 41),
  noExtinction = c("Npy+" = 20, "Npy-" = 134 - 20)
)

# 输出 PDF 文件
pdf("plot/pie_plot_npy_expressing_percent_downsample_seed_722_IEG_pos_Inh_neuron.pdf", width = 7, height = 4)

# 两张图并排
par(mfrow = c(1, 2))

# 绘图
for (group in names(npy_counts)) {
  pie(npy_counts[[group]],
      labels = paste0(names(npy_counts[[group]]), "\n", 
                      round(100 * npy_counts[[group]] / sum(npy_counts[[group]]), 1), "%"),
      col = c("hotpink1", "lavenderblush"),
      main = group)
}

dev.off()
par(mfrow = c(1, 1))

###### Npy expression level (ext vs. noExt) (fig s7) ####
Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722 <- subset(IEGs_positive_Inh_downsample_subset_seed_722, subset = Npy > 0)
Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722$stim <- factor(Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722$stim, levels = c("noExtinction", "extinction"))

Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722_ext <- subset(Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722, idents = "extinction")
Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722_noExt <- subset(Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722, idents = "noExtinction")
Npy_expression_level_Inh_downsample_subset_seed_722_ext <- Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722_ext["Npy",] %>% GetAssayData(., assay = "RNA", slot = "data") 
Npy_expression_level_Inh_downsample_subset_seed_722_noExt <-Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722_noExt["Npy",] %>% GetAssayData(., assay = "RNA", slot = "data") 
# 判断是否符合正态分布
shapiro.test(Npy_expression_level_Inh_downsample_subset_seed_722_ext["Npy",])
# Shapiro-Wilk normality test
# 
# data:  Npy_expression_level_Inh_downsample_subset_seed_722_ext["Npy", ]
# W = 0.95826, p-value = 0.1366
shapiro.test(Npy_expression_level_Inh_downsample_subset_seed_722_noExt["Npy",])
# Shapiro-Wilk normality test
# 
# data:  Npy_expression_level_Inh_downsample_subset_seed_722_noExt["Npy", ]
# W = 0.90375, p-value = 0.04851

# 判断方差是否相等
var.test(Npy_expression_level_Inh_downsample_subset_seed_722_ext["Npy",], Npy_expression_level_Inh_downsample_subset_seed_722_noExt["Npy",])
# F test to compare two variances
# 
# data:  Npy_expression_level_Inh_downsample_subset_seed_722_ext["Npy", ] and Npy_expression_level_Inh_downsample_subset_seed_722_noExt["Npy", ]
# F = 1.3939, num df = 40, denom df = 19, p-value = 0.4409
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.5974953 2.9082066
# sample estimates:
#   ratio of variances 
# 1.393911 

t.test(Npy_expression_level_Inh_downsample_subset_seed_722_ext, Npy_expression_level_Inh_downsample_subset_seed_722_noExt)
# Welch Two Sample t-test
# 
# data:  Npy_expression_level_Inh_downsample_subset_seed_722_ext and Npy_expression_level_Inh_downsample_subset_seed_722_noExt
# t = 2.3839, df = 43.967, p-value = 0.02151
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1084241 1.2945309
# sample estimates:
#   mean of x mean of y 
# 2.857191  2.155713 

wilcox.test(Npy_expression_level_Inh_downsample_subset_seed_722_ext["Npy",], Npy_expression_level_Inh_downsample_subset_seed_722_noExt["Npy",], alternative = "two.sided")
# Wilcoxon rank sum exact test
# 
# data:  Npy_expression_level_Inh_downsample_subset_seed_722_ext["Npy", ] and Npy_expression_level_Inh_downsample_subset_seed_722_noExt["Npy", ]
# W = 552, p-value = 0.02885
# alternative hypothesis: true location shift is not equal to 0
pdf("plot/vlnplot_Npy_expression_level_IEGs_positive_Inh_downsample_seed_722.pdf")
VlnPlot(Npy_positive_IEGs_positive_Inh_downsample_subset_seed_722, 
        features = "Npy", 
        pt.size = 1, 
        group.by = "stim")
dev.off()
