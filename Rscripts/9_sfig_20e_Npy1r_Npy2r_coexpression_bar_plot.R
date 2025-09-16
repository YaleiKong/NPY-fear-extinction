###### Npyr no-overlap bar plot (Exc + Inh) ####
Idents(extinction.integrated_data) <- "celltype"
neuron <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"))
# 8338 samples
no_npy1r_npy2r_neuron <- subset(neuron, subset = Npy1r == 0 & Npy2r == 0)
# 7383
only_npy1r_neuron <- subset(neuron, subset = Npy1r > 0 & Npy2r == 0)
# 432
only_npy2r_neuron <- subset(neuron, subset = Npy2r > 0 & Npy1r == 0)
# 487
npy1r_npy2r_neuron <- subset(neuron, subset = Npy1r > 0 & Npy2r > 0)
# 36 samples

npy1r_npy2r_expressing_cellNumber_neuron_mat <- matrix(c(7383, 432, 487, 36), nrow = 2, byrow = TRUE)
colnames(npy1r_npy2r_expressing_cellNumber_neuron_mat) <- c("Npy1r_Negative", "Npy1r_Positive")
rownames(npy1r_npy2r_expressing_cellNumber_neuron_mat)<- c("Npy2r_Negative", "Npy2r_Positive")
# Npy1r_Negative Npy1r_Positive
# Npy2r_Negative           7383            432
# Npy2r_Positive            487             36

chisq.test(npy1r_npy2r_expressing_cellNumber_neuron_mat)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  npy1r_npy2r_expressing_cellNumber_neuron_mat
# X-squared = 1.4539, df = 1, p-value = 0.2279

# 原假设 Npy1r 和 Npy2r 表达状态是相互独立的（没有显著关联）。
# 备择假设 Npy1r 和 Npy2r 表达状态之间存在显著关联。

npy1r_npy2r_expressing_cellNumber_neuron_mat_long <- as.data.frame(npy1r_npy2r_expressing_cellNumber_neuron_mat) %>%
  rownames_to_column("Npy2r_status") %>%
  pivot_longer(cols = -Npy2r_status, 
               names_to = "Npy1r_status", 
               values_to = "Cell_Count")

pdf("plot/barplot_Npy1r_npy2r_no_overlap_neuron_20250112.pdf")
p <- ggplot(npy1r_npy2r_expressing_cellNumber_neuron_mat_long, aes(x = Npy1r_status, y = Cell_Count, fill = Npy2r_status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.6) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  labs(
    title = "Expression of Npy1r and Npy2r in Neurons \n(p-value = 0.2279) ",
    x = "Npy1r Status", 
    y = "Cell Count", 
    fill = "Npy2r Status"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    
    axis.title = element_text(size = 14, face = "bold"),
    
    plot.title = element_text(size = 16, hjust = 0.5),
    
    legend.position = "bottom",  
    legend.key = element_rect(fill = 'transparent'), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
library(gg.gap)
# 断轴
p_gap <- gg.gap(
  plot = p,
  segments = c(600, 7000), # 定义断轴范围
  tick_width = 200,       # 断轴后两部分之间的刻度间隔
  rel_heights = c(0.2, 0.8), # 定义上下两部分的相对高度
  ylim = c(0, 7500)        # 设置 y 轴范围
)
print(p)
print(p_gap)
dev.off()

