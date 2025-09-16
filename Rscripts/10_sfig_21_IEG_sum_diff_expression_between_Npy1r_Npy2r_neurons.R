#### step 9: vlnplot: IEGs_sum of Npy1r+ vs. Npy2r+ neuron in extinction or noExtinction stim ####
# 将overlap的细胞画两次
Exc_Npy1r_pos <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"), subset = Npy1r > 0)
Exc_Npy1r_pos$Npy1r_Npy2r <- "Npy1r"
Exc_Npy2r_pos <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"), subset = Npy2r > 0)
Exc_Npy2r_pos$Npy1r_Npy2r <- "Npy2r"
Exc_Npy1r_Npy2r_pos <- merge(Exc_Npy1r_pos, Exc_Npy2r_pos)
Exc_Npy1r_Npy2r_pos$IEGs_sum <- colSums(GetAssayData(Exc_Npy1r_Npy2r_pos, slot = "data")[IEG_list, ])

Idents(Exc_Npy1r_Npy2r_pos) <- "stim"
Exc_Npy1r_Npy2r_pos$stim <- factor(Exc_Npy1r_Npy2r_pos$stim, levels = c("noExtinction", "extinction"))
Exc_Npy1r_Npy2r_pos@meta.data %>%
  group_by(stim) %>% 
  group_map(~ t.test(IEGs_sum ~ Npy1r_Npy2r, data = .x), .keep = TRUE)
# [[1]]
# 
# Welch Two Sample t-test
# 
# data:  IEGs_sum by Npy1r_Npy2r
# t = 0.75035, df = 324.48, p-value = 0.4536
# alternative hypothesis: true difference in means between group Npy1r and group Npy2r is not equal to 0
# 95 percent confidence interval:
#   -0.6709618  1.4983601
# sample estimates:
#   mean in group Npy1r mean in group Npy2r 
# 18.97992            18.56622 
# 
# 
# [[2]]
# 
# Welch Two Sample t-test
# 
# data:  IEGs_sum by Npy1r_Npy2r
# t = 2.7898, df = 576.4, p-value = 0.005449
# alternative hypothesis: true difference in means between group Npy1r and group Npy2r is not equal to 0
# 95 percent confidence interval:
#   0.4698063 2.7048969
# sample estimates:
#   mean in group Npy1r mean in group Npy2r 
# 17.93633            16.34898 

pdf("plot/vlnplot_IEGs_sum_Npy1r_neuron_vs_Npy2r_neuron_in_noExt_or_ext.pdf", width = 10, height = 5)
VlnPlot(Exc_Npy1r_Npy2r_pos, features = "IEGs_sum", split.by = "Npy1r_Npy2r", group.by = "stim")
dev.off()