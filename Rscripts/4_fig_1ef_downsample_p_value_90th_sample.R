#### seed (p-value 90th) ####
noExt_subset <- Inh_noExt_subset
p_values <- numeric(1000)
for (seed in 1:1000) {
  ext_subset <- Inh_ext_downsample_subset_Npy_expressed_list[[seed]]
  # 获取Npy表达数据
  ext_expr <- FetchData(ext_subset, vars = "Npy")$Npy
  noExt_expr <- FetchData(noExt_subset, vars = "Npy")$Npy
  
  # t.test
  t_res <- t.test(ext_expr, noExt_expr)
  p_values[seed] <- t_res$p.value
}
p_90th <- quantile(p_values, 0.9, na.rm = TRUE)
seed_90th <- which.min(abs(p_values - p_90th))
# 736
p_values[736]

#### seed 736 ####
Idents(extinction.integrated_data) <- "celltype.stim"
Inh_noExt_cells <- WhichCells(extinction.integrated_data, idents = "Inhibitory neuron_noExtinction")
Inh_ext_cells <- WhichCells(extinction.integrated_data, idents = "Inhibitory neuron_extinction")
min_Inh_cells <- min(length(Inh_noExt_cells), length(Inh_ext_cells))
Inh_ext_subset <- subset(extinction.integrated_data, idents = "Inhibitory neuron_extinction")
Inh_noExt_subset <- subset(extinction.integrated_data, idents = "Inhibitory neuron_noExtinction")

set.seed(736)
Inh_ext_cells_downsample_seed_736 <- sample(Inh_ext_cells, min_Inh_cells)
# 抽取子集
Inh_ext_downsample_subset_seed_736 <- subset(
  Inh_ext_subset, 
  cells = Inh_ext_cells_downsample_seed_736
)

Inh_downsample_subset_seed_736 <- merge(Inh_ext_downsample_subset_seed_736, Inh_noExt_subset)
# Npy 阳性子集
Inh_downsample_subset_seed_736_Npy_expressed <- subset(
  Inh_downsample_subset_seed_736, 
  subset = Npy > 0
)
##### statistics #####
# 提取 Npy 表达数据（log-normalized）
npy_data <- GetAssayData(Inh_downsample_subset_seed_736_Npy_expressed, slot = "data")["Npy", ]

# 拆分为两组
npy_ext <- npy_data[Inh_downsample_subset_seed_736_Npy_expressed$stim == "extinction"]
npy_noExt <- npy_data[Inh_downsample_subset_seed_736_Npy_expressed$stim == "noExtinction"]

# 判断是否符合正态分布
shapiro.test(npy_ext)
# Shapiro-Wilk normality test
# 
# data:  npy_ext
# W = 0.93364, p-value = 0.001766
# 不符合正态分布
shapiro.test(npy_noExt)
# Shapiro-Wilk normality test
# 
# data:  npy_noExt
# W = 0.90635, p-value = 0.006779
# 不符合正态分布

# 判断方差是否相等
var.test(npy_ext, npy_noExt)
# F test to compare two variances
# 
# data:  npy_ext and npy_noExt
# F = 1.6701, num df = 64, denom df = 33, p-value = 0.1099
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.887972 2.962694
# sample estimates:
#   ratio of variances 
# 1.67007
# 方差相等

# 进行 t 检验
t.test(x = as.numeric(npy_ext),
       y = as.numeric(npy_noExt),
       alternative = "two.sided",
       var.equal = TRUE)

# Two Sample t-test
# 
# data:  as.numeric(npy_ext) and as.numeric(npy_noExt)
# t = 2.9742, df = 97, p-value = 0.003707
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.2635734 1.3209385
# sample estimates:
#   mean of x mean of y 
# 3.011797  2.219541

wilcox.test(npy_ext, npy_noExt, alternative = "two.sided")
# Wilcoxon rank sum test with continuity correction
# 
# data:  npy_ext and npy_noExt
# W = 1495, p-value = 0.004103
# alternative hypothesis: true location shift is not equal to 0

##### vln plot ####
Idents(Inh_downsample_subset_seed_736_Npy_expressed) <- "stim"
Inh_downsample_subset_seed_736_Npy_expressed$stim <- factor(Inh_downsample_subset_seed_736_Npy_expressed$stim, levels = c("noExtinction","extinction"))
pdf("plot/violinplot_npy_expression_downsample_seed_736_inh_neuron_split_by_stim.pdf", width = 5, height = 5)
VlnPlot(object = Inh_downsample_subset_seed_736_Npy_expressed, 
        group.by = "stim", 
        features = "Npy", 
        pt.size = 0.5, 
        cols = alpha(c("brown1", "royalblue3"), 0.8)) &  
  NoLegend() 
dev.off()
##### pie plot ####
Idents(Inh_downsample_subset_seed_736) <- "stim"
npy_expressing_info <- data.frame(expressing_number = as.numeric(table(Idents(Inh_downsample_subset_seed_736_Npy_expressed))),
                                  all_number = as.numeric(table(Idents(Inh_downsample_subset_seed_736))),
                                  row.names = unique(Inh_downsample_subset_seed_736@active.ident)) %>% 
  mutate(expressing_percent = expressing_number/all_number) %>% 
  mutate(no_expressing_percent = 1-expressing_percent)
# expressing_number all_number expressing_percent no_expressing_percent
# extinction                  65        278          0.2338129             0.7661871
# noExtinction                34        278          0.1223022             0.8776978
pdf("plot/pie_plot_npy_expressing_percent_downsample_seed_736_Inh_neuron.pdf", width = 7, height = 4)

# 设置布局：两行一列
par(mfrow = c(1, 2))  # side-by-side

# extinction
values_ext <- t(npy_expressing_info)[3:4, 1]
labels_ext <- paste0(c("Npy+", "Npy-"), "\n", round(100 * values_ext, 1), "%")
pie(values_ext,
    labels = labels_ext,
    col = c("hotpink1", "lavenderblush"),
    main = "extinction")

# noExtinction
values_noExt <- t(npy_expressing_info)[3:4, 2]
labels_noExt <- paste0(c("Npy+", "Npy-"), "\n", round(100 * values_noExt, 1), "%")
pie(values_noExt,
    labels = labels_noExt,
    col = c("hotpink1", "lavenderblush"),
    main = "noExtinction")

dev.off()
par(mfrow = c(1, 1))
