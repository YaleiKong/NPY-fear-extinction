#### Step 7: Inh subcluster ####
Idents(extinction.integrated_data) <- "celltype"
Inh <- subset(extinction.integrated_data, idents = "Inhibitory neuron")
table(Inh$stim)
# extinction noExtinction 
# 1023          278 
Inh <- NormalizeData(object = Inh,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) %>% 
  ScaleData()
Inh <- FindVariableFeatures(Inh)
Inh <- RunPCA(Inh)
ElbowPlot(Inh)
Inh <- RunUMAP(Inh, dims = 1:11)

##### resolution = 0.1 #####
Inh <- FindClusters(Inh, resolution = 0.1)
# modularity.fxn	Modularity function (1 = standard; 2 = alternative).
# Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). 
# Leiden requires the leidenalg python.
Inh$inh_number_renew <- paste0("Inh_sub_", Inh$RNA_snn_res.0.1)
Inh$inh_number_renew <- factor(
  Inh$inh_number_renew,
  levels = c("Inh_sub_0", "Inh_sub_1", "Inh_sub_2", "Inh_sub_3", "Inh_sub_4")
)
Inh$stim <- factor(Inh$stim, levels = c("noExtinction", "extinction"))
save(Inh, file = "rdata/Inh_after_findcluster_res.0.1.rdata")
load("rdata/Inh_after_findcluster_res.0.1.rdata")

table(Inh$inh_number_renew, Inh$stim)
Inh_subcluster_size <- table(Inh$inh_number_renew, Inh$stim) %>% as.data.frame() %>% 
  pivot_wider(
    data = .,
    names_from = Var2,
    values_from = Freq
  ) %>% 
  mutate(extinction_percent = extinction/sum(extinction),
         noExtinction_percent = noExtinction/sum(noExtinction))
# A tibble: 5 × 5
# Var1      noExtinction extinction extinction_percent noExtinction_percent
# <fct>            <int>      <int>              <dbl>                <dbl>
# 1 Inh_sub_0           68        339             0.331                0.245 
# 2 Inh_sub_1          116        266             0.260                0.417 
# 3 Inh_sub_2           49        214             0.209                0.176 
# 4 Inh_sub_3           33        125             0.122                0.119 
# 5 Inh_sub_4           12         79             0.0772               0.0432

###### umap (subcluster) #####
pdf("plot/umap_inh_neuron_subclusters_resolution_0.1.pdf", width = 6, height = 5)
DimPlot(Inh, 
        reduction = "umap", 
        group.by = "inh_number_renew", 
        label = TRUE, 
        pt.size = 1,
        cols = c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462')) +
  ggtitle("Inhibitory Neurons subclusters") + 
  theme(plot.title = element_text(size = 14))  
dev.off()
###### vlnplot (inh marker) (fig s8a) ######
pdf("plot/vlnplot_inh_markers_inh_neuron_subclusters_resolution_0.1.pdf", width = 5, height = 10)
VlnPlot(Inh, 
        features = c("Slc32a1","Gad1","Gad2","Cck","Vip","Npy","Sst","Calb2","Lhx6"),
        group.by = "inh_number_renew", 
        cols = c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462'),
        pt.size = 0.1, ncol=1, combine = T)

VlnPlot(Inh, 
        features = c("Slc32a1","Gad1","Gad2","Cck","Vip","Npy","Sst","Calb2","Lhx6"),
        group.by = "inh_number_renew", 
        cols = c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462'),
        pt.size = 0.1, ncol=1, combine = T) & 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 8)) & 
  xlab(NULL) & 
  ylab(NULL)  

VlnPlot(Inh, 
        features = c("Slc32a1","Gad1","Gad2","Cck","Vip","Npy","Sst","Calb2","Lhx6"),
        group.by = "inh_number_renew", 
        cols = c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462'),
        pt.size = 0, ncol=1, combine = T) & 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 8)) & 
  xlab(NULL) & 
  ylab(NULL) 
dev.off()

###### statistics #####
load("rdata/Inh_after_findcluster_res.0.1.rdata")
meta_df <- Inh@meta.data
Inh_subclusters_Npy_expression <- GetAssayData(Inh, slot = "data")["Npy", ]
meta_df$Npy_expr <- Inh_subclusters_Npy_expression

Inh_subclusters_Npy_expression_t_test <- meta_df %>%
  group_by(inh_number_renew) %>%
  summarise(
    t_test_result = list(t.test(Npy_expr ~ stim)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(t_test_result, ~.$p.value),
    mean_Npy_ext = map_dbl(t_test_result, ~.$estimate[2]),
    mean_Npy_noExt = map_dbl(t_test_result, ~.$estimate[1])
  ) %>%
  select(inh_number_renew, mean_Npy_ext, mean_Npy_noExt, p_value)
# inh_number_renew mean_Npy_ext mean_Npy_noExt    p_value
# <chr>                   <dbl>          <dbl>      <dbl>
#   1 Inh_sub_0               1.08          0.419  0.00000461
# 2 Inh_sub_1               0.338         0.0690 0.0000696 
# 3 Inh_sub_2               1.07          0.795  0.223     
# 4 Inh_sub_3               0.127         0      0.0266    
# 5 Inh_sub_4               0.467         0      0.00105   
write.csv(Inh_subclusters_Npy_expression_t_test, file = "file/Inh_subclusters_Npy_expression_t_test.csv")

##### Exc cluster Npy expression statistics #####
Idents(extinction.integrated_data) <- "celltype_number"
Exc_cluster_subset <- subset(extinction.integrated_data, idents = Exc_clusters)
meta_df <- Exc_cluster_subset@meta.data
Exc_cluster_subset_subclusters_Npy_expression <- GetAssayData(Exc_cluster_subset, slot = "data")["Npy", ]
meta_df$Npy_expr <- Exc_cluster_subset_subclusters_Npy_expression

Exc_cluster_subset_subclusters_Npy_expression_t_test <- meta_df %>%
  group_by(celltype_number) %>%
  summarise(
    t_test_result = list(t.test(Npy_expr ~ stim)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(t_test_result, ~.$p.value),
    mean_Npy_ext = map_dbl(t_test_result, ~.$estimate[1]),
    mean_Npy_noExt = map_dbl(t_test_result, ~.$estimate[2])
  ) %>%
  select(celltype_number, mean_Npy_ext, mean_Npy_noExt, p_value)
# # A tibble: 8 × 4
# celltype_number mean_Npy_ext mean_Npy_noExt    p_value
# <fct>                  <dbl>          <dbl>      <dbl>
#   1 Exc_3                 0.0639         0.0208 0.000295  
# 2 Exc_4                 0.0699         0.0169 0.00000332
# 3 Exc_10                0.0778         0.0460 0.269     
# 4 Exc_16                0.0659         0.0216 0.0260    
# 5 Exc_27                0.0701         0.0497 0.548     
# 6 Exc_28                0.112          0.0131 0.000456  
# 7 Exc_30                0.0861         0.0275 0.143     
# 8 Exc_32                0.106          0.0415 0.209
write.csv(Exc_cluster_subset_subclusters_Npy_expression_t_test, file = "file/Exc_cluster_subset_subclusters_Npy_expression_t_test.csv")
