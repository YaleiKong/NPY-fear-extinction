# note:
# 1.percent.mt 5
# 2.seurat integrate 
# 3.dims=1:30
# 4.resolution 0.5

.libPaths(new = "/lustre/home/acct-medlqian/medlqian-loop3/R/x86_64-redhat-linux-gnu-library/4.1.3")
rm(list=ls())
set.seed(1)

setwd('/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/fear_extinction_10X_20190403/')

#### packages ####
library(Seurat) # V4.3.0
library(MAST) 

library(Matrix)
library(tidyverse)
library(dplyr)

library(ggplot2)
library(ggrepel) # organize the labels nicely using the "ggrepel" package and the geom_text_repel() function #to repel overlapping text 
library(ggpubr) 
library(cowplot) # similar to ggplot2’s theme_classic()
library(RColorBrewer)
library(gridExtra)
library(plotrix) # 3D pie
library(reshape2) 
library(pheatmap)

library(openxlsx)
library(readxl)

library (VennDiagram) 
library(ggvenn)
library(scales)  # for alpha()

library(future)

library(dplyr)
library(purrr)

#### step 1: Create integrated seurat object ####
baseDir <- paste0(getwd(), "/raw_data")
noExtinction_dataDir <- paste(baseDir, "/TX1811XTHJ1X_control", sep="")
extinction_dataDir <- paste(baseDir, "/TX1811XTHJ2X_fear", sep="")

noExtinction.extinction.data <- Read10X(data.dir = paste0(noExtinction_dataDir, "/outs/filtered_feature_bc_matrix"))
extinction.extinction.data <- Read10X(data.dir = paste0(extinction_dataDir, "/outs/filtered_feature_bc_matrix"))

noExtinction.extinction.data <- noExtinction.extinction.data[, Matrix::colSums(noExtinction.extinction.data)>=700 & Matrix::colSums(noExtinction.extinction.data) <=15000]
extinction.extinction.data <- extinction.extinction.data[, Matrix::colSums(extinction.extinction.data)>=700 & Matrix::colSums(extinction.extinction.data) <=15000]

noExtinction.extinction <- CreateSeuratObject(counts = noExtinction.extinction.data,
                                              min.cells = 3,
                                              min.features = 200,
                                              project = "10X_noExtinction")
noExtinction.extinction$stim <- "noExtinction"
noExtinction.extinction
# An object of class Seurat 
# 18900 features across 19556 samples within 1 assay 
# Active assay: RNA (18900 features, 0 variable features)

extinction.extinction <- CreateSeuratObject(counts = extinction.extinction.data,
                                            min.cells = 3,
                                            min.features = 200,
                                            project = "10X_extinction")
extinction.extinction$stim <- "extinction"
extinction.extinction 
# An object of class Seurat 
# 18628 features across 25728 samples within 1 assay 
# Active assay: RNA (18628 features, 0 variable features)

noExtinction.extinction[["percent.mt"]] <- PercentageFeatureSet(object=noExtinction.extinction, pattern='^mt-')
extinction.extinction[["percent.mt"]] <- PercentageFeatureSet(object=extinction.extinction, pattern='^mt-')

plot1 <- FeatureScatter(object = noExtinction.extinction, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = noExtinction.extinction, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

noExtinction.extinction <- subset(x = noExtinction.extinction, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# An object of class Seurat 
# 18900 features across 13281 samples within 1 assay 
# Active assay: RNA (18900 features, 0 variable features)

extinction.extinction <- subset(x = extinction.extinction, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# An object of class Seurat 
# 18628 features across 19811 samples within 1 assay 
# Active assay: RNA (18628 features, 0 variable features)

noExtinction.extinction <- NormalizeData(object = noExtinction.extinction,
                                         normalization.method = "LogNormalize",
                                         scale.factor = 1e4)

extinction.extinction <- NormalizeData(object = extinction.extinction,
                                       normalization.method = "LogNormalize",
                                       scale.factor = 1e4)

noExtinction.extinction <- FindVariableFeatures(object = noExtinction.extinction, selection.method = "vst", nfeatures = 2000)
extinction.extinction <- FindVariableFeatures(object = extinction.extinction, selection.method = "vst", nfeatures = 2000)

extinction.integrated_data.anchors <- FindIntegrationAnchors(object.list = list(noExtinction.extinction, extinction.extinction), dims = 1:30)
extinction.integrated_data <- IntegrateData(anchorset = extinction.integrated_data.anchors, dims = 1:30)
#### step 2: pre-processing ####
DefaultAssay(extinction.integrated_data) <- "integrated"
extinction.integrated_data <- ScaleData(extinction.integrated_data)
extinction.integrated_data <- RunPCA(extinction.integrated_data, npcs = 100)
extinction.integrated_data <- JackStraw(object = extinction.integrated_data, dims = 100, num.replicate = 100)
extinction.integrated_data <- ScoreJackStraw(object = extinction.integrated_data, dims = 1:100)

JackStrawPlot(object = extinction.integrated_data, dims = 1:100)
dev.copy2pdf(width=15, height=7,file= 'plot/JackStrawPlot_50PCs_integrated_data.pdf')
dev.off()

ElbowPlot(object = extinction.integrated_data, ndims = 100)
dev.copy2pdf(width=15, height=7,file= 'plot/ElbowPlot_50PCs_integrated_data.pdf')
dev.off()
#### step 3: Clustering (pca dims = 1:30; resolution = 0.5) ####
extinction.integrated_data <- FindNeighbors(extinction.integrated_data, reduction = "pca", dims = 1:30) 
extinction.integrated_data <- FindClusters(extinction.integrated_data, resolution = 0.5) # adjust the value of resolution bigger to have more communities 
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
extinction.integrated_data <- RunTSNE(extinction.integrated_data, reduction = "pca", dims = 1:30)
extinction.integrated_data <- RunUMAP(extinction.integrated_data, reduction = "pca", dims = 1:30)
save(extinction.integrated_data, file = "extinction.integrated_data_after_findcluster.rdata")
load("extinction.integrated_data_after_findcluster.rdata")

DefaultAssay(extinction.integrated_data) <- "RNA"
integrated.extinction.FindAllMarkers <- FindAllMarkers(object = extinction.integrated_data,
                                                       only.pos = TRUE,
                                                       min.pct = 0.25,
                                                       logfc.threshold = 0.25)
integrated.extinction.FindAllMarkers %>% group_by(cluster) %>% top_n(2, avg_logFC)
#### step 4: Assigning cell type identity to clusters ####
# umap with cluster
pdf("umap_with_cluster.pdf", width=7, height=7)
DimPlot(extinction.integrated_data, reduction = "umap", pt.size = 1, label = T)
dev.off()

# featureplot
list_featurePlot <- list()
allMarkers <- c("Gapdh",
                "Thy1","Slc17a6","Slc17a7",# excitatory neoron marker
                "Slc32a1","Gad1","Gad2","Cck","Vip","Npy","Sst","Calb2",# inhibitory neuron marker
                "Gfap","Aldh1l1","Aldoc","S100b",# Astrocyte marker
                "Olig1","Olig2","Sox10","Opalin","Mbp","Cspg4",# oligodendrocyte marker
                "Pdgfra", # OPC
                "Cx3cr1","Csf1r","Tmem119","P2ry12",# macroglia
                "Vtn","Col1a1","Col1a2","Pdgfrb","Rgs5",#pericyte
                "Cldn5","Cxcl12",#endothelium, smooth muscle
                "Mrc1") # macrophage

pdf("FeaturePlot_all_neuron_markers_umap.pdf", width=7, height=7)
for (i in c(1:35))
{list_featurePlot[[i]] <- FeaturePlot(object = extinction.integrated_data,
                                      features = allMarkers[i],
                                      reduction = "umap")
}
print(list_featurePlot)
dev.off()

# violinPlot (all markers)
pdf("VlnPlot_all_neuron_markers.pdf",width = 35,height = 100)
plot <- VlnPlot(extinction.integrated_data, 
                features = c("Gapdh",#
                             "Thy1","Slc17a6","Slc17a7",# excitatory neoron marker
                             "Slc32a1","Gad1","Gad2","Cck","Vip","Npy","Sst","Calb2",# inhibitory neuron marker
                             "Gfap","Aldh1l1","Aldoc","S100b",# Astrocyte marker
                             "Olig1","Olig2","Sox10","Opalin","Mbp","Cspg4",# oligodendrocyte marker
                             "Pdgfra", # OPC
                             "Cx3cr1","Csf1r","Tmem119","P2ry12",# macroglia
                             "Vtn","Col1a1","Col1a2","Pdgfrb","Rgs5",#pericyte
                             "Cldn5","Cxcl12",#endothelium, smooth muscle
                             "Mrc1"),#macrophage
                pt.size = 0, 
                combine = FALSE, 
                ncol=1)
print(plot)
dev.off()

##### add metadata ##### 
###### celltype ######
extinction.integrated_data <- RenameIdents(extinction.integrated_data, 
                                           '3' = "Excitatory neuron", '4' = "Excitatory neuron", '10' = "Excitatory neuron", '16' = "Excitatory neuron", '27' = "Excitatory neuron", '28'="Excitatory neuron", '30'="Excitatory neuron", '32'="Excitatory neuron", 
                                           '6' = "Inhibitory neuron", 
                                           '1' = "Astrocyte", '5' = "Astrocyte", '26' = "Astrocyte", '31' = "Astrocyte", '34' = "Astrocyte",
                                           '0' = "Oligodendrocyte", '14' = "Oligodendrocyte", '20' = "Oligodendrocyte", '22' = "Oligodendrocyte", '33' = "Oligodendrocyte", 
                                           '8' = "OPC", '19' = "OPC", 
                                           '12' = "Microglia", '13' = "Microglia", '24' = "Microglia", '25'= "Microglia", '29' = "Microglia", 
                                           '7' = "Pericyte",'17' = "Pericyte",'18'="Pericyte", '23' = "Pericyte",
                                           '2' = "Endothelium, smooth muscle", '9'= "Endothelium, smooth muscle", '11' = "Endothelium, smooth muscle", '21' = "Endothelium, smooth muscle", 
                                           '15'= "Macrophage")
extinction.integrated_data$celltype <-Idents(extinction.integrated_data) 
extinction.integrated_data@active.ident <- factor(extinction.integrated_data@active.ident, levels = c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage"))

###### celltype_stim ######
extinction.integrated_data$celltype.stim <- paste(Idents(extinction.integrated_data), extinction.integrated_data$stim, sep = "_") # stored in extinction.integrated_data@meta.data
Idents(extinction.integrated_data) <- "celltype.stim"

###### celltype_number ######
Idents(extinction.integrated_data) <- "seurat_clusters"
extinction.integrated_data <- RenameIdents(extinction.integrated_data, 
                                           '3' = "Exc_3", '4' = "Exc_4", '10' = "Exc_10", '16' = "Exc_16", '27' = "Exc_27", '28'="Exc_28", '30'="Exc_30", '32'="Exc_32", 
                                           '6' = "Inh_6", 
                                           '1' = "Astro_1", '5' = "Astro_5", '26' = "Astro_26", '31' = "Astro_31", '34' = "Astro_34",
                                           '0' = "Olig_0", '14' = "Olig_14", '20' = "Olig_20", '22' = "Olig_22", '33' = "Olig_33", 
                                           '8' = "OPC_8", '19' = "OPC_19", 
                                           '12' = "Micro_12", '13' = "Micro_13", '24' = "Micro_24", '25'= "Micro_25", '29' = "Micro_29", 
                                           '7' = "Pericyte_7",'17' = "Pericyte_17",'18'="Pericyte_18", '23' = "Pericyte_23",
                                           '2' = "Endo_2", '9'= "Endo_9", '11' = "Endo_11", '21' = "Endo_21", 
                                           '15'= "Macrophage_15" )
extinction.integrated_data$celltype_number <- Idents(extinction.integrated_data)
extinction.integrated_data@active.ident <- factor(extinction.integrated_data@active.ident, 
                                                  levels = c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                                                             "Inh_6", 
                                                             "Astro_1", "Astro_5", "Astro_26", "Astro_31", "Astro_34",
                                                             "Olig_0", "Olig_14", "Olig_20", "Olig_22", "Olig_33", 
                                                             "OPC_8", "OPC_19", 
                                                             "Micro_12", "Micro_13", "Micro_24", "Micro_25", "Micro_29", 
                                                             "Pericyte_7", "Pericyte_17", "Pericyte_18", "Pericyte_23",
                                                             "Endo_2", "Endo_9", "Endo_11", "Endo_21", 
                                                             "Macrophage_15"))
###### celltype_number.stim ######
extinction.integrated_data$celltype_number.stim <- paste(Idents(extinction.integrated_data), extinction.integrated_data$stim, sep = "_")

###### seurat_clusters ######
Idents(extinction.integrated_data) <- "seurat_clusters"
extinction.integrated_data@active.ident <- factor(extinction.integrated_data@active.ident, 
                                                  levels = c(3, 4, 10, 16, 27, 28, 30, 32,
                                                             6,
                                                             1, 5, 26, 31, 34, 
                                                             0, 14, 20, 22, 33, 
                                                             8, 19, 
                                                             12, 13, 24, 25, 29, 
                                                             7, 17, 18, 23, 
                                                             2, 9, 11, 21, 
                                                             15))

save(extinction.integrated_data,file="extinction.integrated_data_with_metadata.Rdata")
load("extinction.integrated_data_with_metadata.Rdata")

##### celltype size between ext and noExt (fig s4) ####
load("rdata/extinction.integrated_data_with_metadata.Rdata")
Idents(extinction.integrated_data) <- "stim"
extinction <- subset(extinction.integrated_data, idents = "extinction")
noExtinction <- subset(extinction.integrated_data, idents = "noExtinction")
Idents(extinction) <- "celltype"
Idents(noExtinction) <- "celltype"

celltype <- c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage")
extinction_celltype_cellnumber <- c()
noExtinction_celltype_cellnumber <- c()
for (celltype_order in 1:length(celltype)) {
  extinction_celltype_cellnumber[celltype_order] <- ncol(subset(extinction, idents = celltype[celltype_order]))
  noExtinction_celltype_cellnumber[celltype_order] <- ncol(subset(noExtinction, idents = celltype[celltype_order]))
}
stim_celltype_cellnumber <- data.frame(
  celltype = celltype,
  extinction_celltype_cellnumber = extinction_celltype_cellnumber,
  noExtinction_celltype_cellnumber = noExtinction_celltype_cellnumber
)
stim_celltype_cellnumber_percent <- stim_celltype_cellnumber %>%
  mutate(
    extinction_celltype_cellnumber_percent = extinction_celltype_cellnumber / sum(extinction_celltype_cellnumber),
    noExtinction_celltype_cellnumber_percent = noExtinction_celltype_cellnumber / sum(noExtinction_celltype_cellnumber)
  ) %>%
  select(celltype, extinction_celltype_cellnumber_percent, noExtinction_celltype_cellnumber_percent)

# melt: reshape2
stim_celltype_cell_number_percent_long <- melt(stim_celltype_cellnumber_percent, 
                                               id.vars = "celltype",
                                               variable.name = "condition",
                                               value.name = "cell_number_percent") %>% 
  mutate(celltype = factor(celltype, levels = c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage")),
         condition = factor(condition, levels = c("noExtinction_celltype_cellnumber_percent", "extinction_celltype_cellnumber_percent")))

pdf("plot/stack_plot_stim_celltype_celltype_percent.pdf", width = 5, height = 6)
ggplot(stim_celltype_cell_number_percent_long, aes(x = condition, y = cell_number_percent, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("brown1","sienna","mediumblue", "dodgerblue2","skyblue2", "darkturquoise","springgreen3","seagreen","mediumpurple3")) +
  theme_classic() +
  labs(title = "Stacked Bar Plot of Cell Numbers (Percent)",
       x = "Condition",
       y = "Cell Number Percent (%)",
       fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_x_discrete(labels = c("noExtinction", "extinction")) 
dev.off()

#### step 5: annotation and visualization (umap, violinPlot, dotplot) ####
celltype_color <- c("brown1","sienna","mediumblue", "dodgerblue2","skyblue2", "darkturquoise","springgreen3","seagreen","mediumpurple3")
stim_color <- c("brown1", "royalblue3")
##### umap with celltype and stimulation (fig 1c) #####
pdf(width=16, height=7,file= 'umapPlot_resolution_0.5.pdf')
p1 <- DimPlot(extinction.integrated_data, reduction = "umap", group.by = "celltype", label = TRUE, 
              cols= c("brown1","sienna","mediumblue", "dodgerblue2","skyblue2", "darkturquoise","springgreen3","seagreen","mediumpurple3"),
              pt.size = 1)
p2 <- DimPlot(extinction.integrated_data, reduction = "umap", group.by = "stim", 
              label = TRUE, 
              cols= c("brown1","royalblue3"),
              pt.size = 1,
              shuffle = T)
p <- p1+p2
print(p)
dev.off()

##### violin plot without x_axis_text according to the order of celltype ###### 
Idents(extinction.integrated_data) <- "celltype_number"
pdf(file= "VlnPlot_selected_markers_with_no_axis.pdf",width = 20, height = 15)
plot <- VlnPlot(extinction.integrated_data, 
                features = c("Gapdh",
                             "Thy1", "Slc17a7", #excitatory neuron marker
                             "Slc32a1", "Gad1" , #inhibitory neuron marker
                             "Gfap", #astrocyte marker
                             "Olig1", #oligodendrocyte marker
                             "Pdgfra", #OPC marker
                             "Cx3cr1", #microglia marker
                             "Vtn", #pericyte marker
                             "Cldn5", #endothelium, smooth muscle marker
                             "Mrc1"), #macrophage marker 
                group.by = "celltype_number", 
                cols=c(rep("brown1", 8), rep("sienna", 1), rep("mediumblue", 5), rep("dodgerblue2",5), rep("skyblue2",2), rep("darkturquoise", 5), rep("springgreen3", 4), rep("seagreen", 4), "mediumpurple3"),
                pt.size = 0, ncol=1, combine = T)
plot <- plot & 
  theme(axis.title.x = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x = element_blank()) & 
  xlab(NULL) & 
  ylab(NULL)
print(plot)
dev.off()

##### violin plot with celltype_number name ######
Idents(extinction.integrated_data) <- "celltype_number"
pdf(width=20, height=30,file= "VlnPlot_with_celltype_number_selected_markers.pdf")
plot <- VlnPlot(extinction.integrated_data, 
                features = c("Gapdh",
                             "Thy1", "Slc17a7", #excitatory neuron marker
                             "Slc32a1", "Gad1" , #inhibitory neuron marker
                             "Gfap", #astrocyte marker
                             "Olig1", #oligodendrocyte marker
                             "Pdgfra", #OPC marker
                             "Cx3cr1", #microglia marker
                             "Vtn", #pericyte marker
                             "Cldn5", #endothelium, smooth muscle marker
                             "Mrc1"), #macrophage marker 
                group.by = "celltype_number", 
                cols = c(rep("brown1", 8), rep("sienna", 1), rep("mediumblue", 5), rep("dodgerblue2",5), rep("skyblue2",2), rep("darkturquoise", 5), rep("springgreen3", 4), rep("seagreen", 4), "mediumpurple3"),
                pt.size = 0, ncol=1, combine = T)
plot <- plot & ylab(NULL)
print(plot)
dev.off()

##### dotplot of markers with celltype_number (fig 3a) #####
list_allMarkers <- c("Gapdh",
                     "Thy1","Slc17a7",# excitatory neoron marker
                     "Slc32a1","Gad1","Gad2","Cck", # inhibitory neuron marker
                     "Gfap","Aldoc","S100b",# Astrocyte marker
                     "Olig1","Olig2","Sox10","Opalin","Mbp","Pdgfra", # oligodendrocyte marker
                     "Cx3cr1","Csf1r","Tmem119","P2ry12",# macroglia
                     "Vtn","Col1a1","Col1a2","Rgs5",#pericyte
                     "Cldn5","Cxcl12",#endothelium, smooth muscle
                     "Mrc1")#macrophage

Idents(extinction.integrated_data) <- "celltype_number"
pdf(width=12, height=8,file= "dotplot_with_cluster_selected_markers_resolution_0.5.pdf")
p <- DotPlot(extinction.integrated_data,features =rev(list_allMarkers),
             group.by = "celltype_number",
             cols = c("lightgrey", "dodgerblue3")) + 
  RotatedAxis()+
  coord_flip () 
print(p)
dev.off()
