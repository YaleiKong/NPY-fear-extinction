#### packages ####
library(Seurat) # 4.3.0
library(DESeq2) # 3.14
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

library(xlsx)
library(readxl)

library (VennDiagram) 
library(ggvenn)

#### Create integrated seurat object ####
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

extinction.extinction <- CreateSeuratObject(counts = extinction.extinction.data,
                                            min.cells = 3,
                                            min.features = 200,
                                            project = "10X_extinction")
extinction.extinction$stim <- "extinction"
extinction.extinction 

noExtinction.extinction[["percent.mt"]] <- PercentageFeatureSet(object=noExtinction.extinction, pattern='^mt-')
extinction.extinction[["percent.mt"]] <- PercentageFeatureSet(object=extinction.extinction, pattern='^mt-')

plot1 <- FeatureScatter(object = noExtinction.extinction, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = noExtinction.extinction, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

noExtinction.extinction <- subset(x = noExtinction.extinction, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

extinction.extinction <- subset(x = extinction.extinction, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

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

#### pre-processing ####
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

#### Clustering (pca dims = 1:30; resolution = 0.5) ####
extinction.integrated_data <- FindNeighbors(extinction.integrated_data, reduction = "pca", dims = 1:30) 
extinction.integrated_data <- FindClusters(extinction.integrated_data, resolution = 0.5) # adjust the value of resolution bigger to have more communities 
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

#### Assigning cell type identity to clusters ####

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
extinction.integrated_data@active.ident <- factor(extinction.integrated_data@active.ident, levels = c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
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
extinction.integrated_data@active.ident <- factor(extinction.integrated_data@active.ident, levels = c(3, 4, 10, 16, 27, 28, 30, 32,
                                                                                                      6,
                                                                                                      1, 5, 26, 31, 34, 
                                                                                                      0, 14, 20, 22, 33, 
                                                                                                      8, 19, 
                                                                                                      12, 13, 24, 25, 29, 
                                                                                                      7, 17, 18, 23, 
                                                                                                      2, 9, 11, 21, 
                                                                                                      15))

save(extinction.integrated_data,file="extinction.integrated_data_with_metadata.Rdata")

#### Extended fig 4a: celltype size between ext and noExt ####
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

####fig 1c: umap with celltype and stimulation #####
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

#### Extended fig 3a: dotplot of markers with celltype_number ####
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

#### DE analysis ####
##### DE analysis celltype (MAST) ####
load('extinction.integrated_data_with_metadata.Rdata')
Idents(extinction.integrated_data) <- "celltype"
celltype_subset <- list()
avg_celltype <- list()
DEGs_MAST <- list()
DEGs_MAST_no_filter <- list()

for (i in 1:9){
  celltype <- c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage")
  Idents(extinction.integrated_data) <- "celltype"
  celltype_subset[[i]] <- subset(extinction.integrated_data, idents = celltype[i])
  Idents(celltype_subset[[i]]) <- "stim"
  avg_celltype[[i]] <- as.data.frame(log1p(AverageExpression(celltype_subset[[i]], verbose = FALSE)$RNA))
  avg_celltype[[i]]$gene <- rownames(avg_celltype[[i]])
  
  celltype_noExtinction <- paste0(celltype,"_noExtinction")
  celltype_extinction <- paste0(celltype,"_extinction")
  Idents(extinction.integrated_data) <- "celltype.stim"
  DEGs_MAST[[i]] <- as.data.frame(FindMarkers(extinction.integrated_data, ident.1 =celltype_extinction[i] , ident.2 = celltype_noExtinction[i], test.use ="MAST", verbose = FALSE))
  avg_celltype[[i]]$pct.1 <- DEGs_MAST[[i]][match(rownames(avg_celltype[[i]]), rownames(DEGs_MAST[[i]])), ]$pct.1
  avg_celltype[[i]]$pct.2 <- DEGs_MAST[[i]][match(rownames(avg_celltype[[i]]), rownames(DEGs_MAST[[i]])), ]$pct.2
  
  # filter: logfc.threshold = 0.1, min.pct = 0.01
  avg_celltype[[i]]$p_adj_MAST <- DEGs_MAST[[i]][match(rownames(avg_celltype[[i]]), rownames( DEGs_MAST[[i]])), ]$p_val_adj
  avg_celltype[[i]]$avg_logFC_MAST <- DEGs_MAST[[i]][match(rownames(avg_celltype[[i]]), rownames(DEGs_MAST[[i]])), ]$avg_log2FC
  
  # nofilter: logfc.threshold = 0, min.pct = 0
  DEGs_MAST_no_filter[[i]] <- as.data.frame(FindMarkers(extinction.integrated_data, ident.1 =celltype_extinction[i] , ident.2 = celltype_noExtinction[i] , test.use ="MAST", verbose = FALSE,
                                                        logfc.threshold = 0,min.pct = 0))
  avg_celltype[[i]]$p_adj_MAST_no_filter <- DEGs_MAST_no_filter[[i]][match(rownames(avg_celltype[[i]]), rownames( DEGs_MAST_no_filter[[i]])), ]$p_val_adj
  avg_celltype[[i]]$avg_logFC_MAST_no_filter <- DEGs_MAST_no_filter[[i]][match(rownames(avg_celltype[[i]]), rownames(DEGs_MAST_no_filter[[i]])), ]$avg_log2FC
  
}

save(avg_celltype,file = "avg_celltype.rdata")
save(DEGs_MAST,file = "DEGs_MAST.rdata")
save(DEGs_MAST_no_filter,file = "DEGs_MAST_no_filter.rdata")

celltype <- c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage")
names(avg_celltype) <- paste0(celltype,"_DEGs")
write.xlsx(avg_celltype,file = "./avg_celltype.xlsx")

##### DE analysis celltype_number (MAST) #####
load('extinction.integrated_data_with_metadata.Rdata')

Idents(extinction.integrated_data) <- "celltype_number"

celltype_number_subset <- list()
avg_celltype_number <- list()
DEGs_MAST_celltype_number <- list()
DEGs_MAST_celltype_number_no_filter <- list()

for (i in 1:35){
  celltype_number <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                       "Inh_6", 
                       "Astro_1","Astro_5","Astro_26","Astro_31", 'Astro_34',
                       "Olig_0", "Olig_14","Olig_20","Olig_22","Olig_33",
                       "OPC_8", "OPC_19", 
                       "Micro_12","Micro_13","Micro_24","Micro_25","Micro_29", 
                       "Pericyte_7","Pericyte_17","Pericyte_18","Pericyte_23",
                       "Endo_2","Endo_9","Endo_11","Endo_21", 
                       "Macrophage_15")
  Idents(extinction.integrated_data) <- "celltype_number"
  celltype_number_subset[[i]] <- subset(extinction.integrated_data, idents = celltype_number[i])
  Idents(celltype_number_subset[[i]]) <- "stim"
  avg_celltype_number[[i]] <- as.data.frame(log1p(AverageExpression(celltype_number_subset[[i]], verbose = FALSE)$RNA))
  avg_celltype_number[[i]]$gene <- rownames(avg_celltype_number[[i]])
  
  celltype_number_noExtinction <- paste0(celltype_number,"_noExtinction")
  celltype_number_extinction <- paste0(celltype_number,"_extinction")
  Idents(extinction.integrated_data) <- "celltype_number.stim"
  DEGs_MAST_celltype_number[[i]] <- as.data.frame(FindMarkers(extinction.integrated_data, ident.1 =celltype_number_extinction[i] , ident.2 = celltype_number_noExtinction[i], test.use ="MAST", verbose = FALSE))
  avg_celltype_number[[i]]$pct.1 <- DEGs_MAST_celltype_number[[i]][match(rownames(avg_celltype_number[[i]]), rownames(DEGs_MAST_celltype_number[[i]])), ]$pct.1
  avg_celltype_number[[i]]$pct.2 <- DEGs_MAST_celltype_number[[i]][match(rownames(avg_celltype_number[[i]]), rownames(DEGs_MAST_celltype_number[[i]])), ]$pct.2
  
  # filter: logfc.threshold = 0.1, min.pct = 0.01
  avg_celltype_number[[i]]$p_adj_MAST <- DEGs_MAST_celltype_number[[i]][match(rownames(avg_celltype_number[[i]]), rownames( DEGs_MAST_celltype_number[[i]])), ]$p_val_adj
  avg_celltype_number[[i]]$avg_logFC_MAST <- DEGs_MAST_celltype_number[[i]][match(rownames(avg_celltype_number[[i]]), rownames(DEGs_MAST_celltype_number[[i]])), ]$avg_log2FC
  
  # filter: logfc.threshold = 0, min.pct = 0
  DEGs_MAST_celltype_number_no_filter[[i]] <- as.data.frame(FindMarkers(extinction.integrated_data, ident.1 =celltype_number_extinction[i] , ident.2 = celltype_number_noExtinction[i] , test.use ="MAST", verbose = FALSE,
                                                                        logfc.threshold = 0,min.pct = 0))
  avg_celltype_number[[i]]$p_adj_MAST_no_filter <- DEGs_MAST_celltype_number_no_filter[[i]][match(rownames(avg_celltype_number[[i]]), rownames( DEGs_MAST_celltype_number_no_filter[[i]])), ]$p_val_adj
  avg_celltype_number[[i]]$avg_logFC_MAST_no_filter <- DEGs_MAST_celltype_number_no_filter[[i]][match(rownames(avg_celltype_number[[i]]), rownames(DEGs_MAST_celltype_number_no_filter[[i]])), ]$avg_log2FC
  
}

save(avg_celltype_number,file = "avg_celltype_number.rdata")
save(DEGs_MAST_celltype_number,file = "DEGs_MAST_celltype_number.rdata")
save(DEGs_MAST_celltype_number_no_filter,file = "DEGs_MAST_celltype_number_no_filter.rdata")

load("rdata/DEGs_MAST_celltype_number_no_filter.rdata")
for (i in 1:8) {
  print(DEGs_MAST_celltype_number_no_filter[[i]][c("Npy1r","Npy2r"),])
}


celltype_number <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                     "Inh_6", 
                     "Astro_1","Astro_5","Astro_26","Astro_31", 'Astro_34',
                     "Olig_0", "Olig_14","Olig_20","Olig_22","Olig_33",
                     "OPC_8", "OPC_19", 
                     "Micro_12","Micro_13","Micro_24","Micro_25","Micro_29", 
                     "Pericyte_7","Pericyte_17","Pericyte_18","Pericyte_23",
                     "Endo_2","Endo_9","Endo_11","Endo_21", 
                     "Macrophage_15")
names(avg_celltype_number) <- paste0(celltype_number,"_DEGs")
write.xlsx(avg_celltype_number,file = "./avg_celltype_number.xlsx")

#### fig 1d and Extended fig 6: volcano plot (DEGs) ####
##### celltype #####
celltype_volcanoPlot_list <- list()
for (i in 1:9) {
  celltype <-c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Microglia","Pericyte","Endothelium, smooth muscle","Macrophage")
  celltype_volcanoPlot_list[[i]] <- ggplot(data = avg_celltype[[i]], aes(x = avg_logFC_MAST_no_filter, y = -log10(p_adj_MAST_no_filter), label = gene, colour = ((p_adj_MAST<0.05 & abs(avg_logFC_MAST) > 0.585)==T))) +
    scale_color_manual(values = c("grey50", "brown1")) +
    geom_point() +
    geom_text_repel(data=subset(avg_celltype[[i]], p_adj_MAST<0.05 & abs(avg_logFC_MAST) > 0.585 ), colour='red') + 
    xlab("log2FC") +
    ylab("-log10(FDR)") + 
    ggtitle(celltype[i]) +
    xlim(-3, 3) + 
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1)),
          axis.text = element_text(size = rel(0.75))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
}
save(celltype_volcanoPlot_list, file = "celltype_volcanoPlot_list.rdata")
pdf("volcanol_plot_celltype_DEGs.pdf", width=7, height=7)
print(celltype_volcanoPlot_list)
dev.off()

##### celltype_number #####
celltype_number_volcanoPlot_list <- list()
for (i in 1:35) {
  celltype_number <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                       "Inh_6", 
                       "Astro_1","Astro_5","Astro_26","Astro_31", 'Astro_34',
                       "Olig_0", "Olig_14","Olig_20","Olig_22","Olig_33",
                       "OPC_8", "OPC_19", 
                       "Micro_12","Micro_13","Micro_24","Micro_25","Micro_29", 
                       "Pericyte_7","Pericyte_17","Pericyte_18","Pericyte_23",
                       "Endo_2","Endo_9","Endo_11","Endo_21", 
                       "Macrophage_15")
  celltype_number_volcanoPlot_list[[i]] <- ggplot(data = avg_celltype_number[[i]], aes(x = avg_logFC_MAST_no_filter, y = -log10(p_adj_MAST_no_filter), label = gene, colour = ((p_adj_MAST<0.05 & abs(avg_logFC_MAST) > 0.585)==T))) +
    scale_color_manual(values = c("grey50", "brown1")) +
    geom_point() +
    geom_text_repel(data=subset(avg_celltype_number[[i]], p_adj_MAST<0.05 & abs(avg_logFC_MAST) > 0.585 ), colour='red') + 
    xlab("log2FC") +
    ylab("-log10(FDR)") + 
    ggtitle(celltype_number[i]) +
    xlim(-3, 3) + 
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5)),
          axis.title = element_text(size = rel(1)),
          axis.text = element_text(size = rel(0.75))) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
}
save(celltype_number_volcanoPlot_list, file = "celltype_number_volcanoPlot_list.rdata")
pdf("volcanol_plot_celltype_number_DEGs.pdf", width=7, height=7)
print(celltype_number_volcanoPlot_list)
dev.off()

#### Extended fig 3b: heatmap plot (DEGs) ####
celltype_number <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                     "Inh_6", 
                     "Astro_1","Astro_5","Astro_26","Astro_31", 'Astro_34',
                     "Olig_0", "Olig_14","Olig_20","Olig_22","Olig_33",
                     "OPC_8", "OPC_19", 
                     "Micro_12","Micro_13","Micro_24","Micro_25","Micro_29", 
                     "Pericyte_7","Pericyte_17","Pericyte_18","Pericyte_23",
                     "Endo_2","Endo_9","Endo_11","Endo_21", 
                     "Macrophage_15")

# select DEGs of all celltype_number and deduplicate
load("rdata/DEGs_MAST_celltype_number.rdata")
DEGs_MAST_celltype_number_sig <- lapply(DEGs_MAST_celltype_number, function(df) {
  df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.585, ]
})

names(DEGs_MAST_celltype_number_sig) <- celltype_number
library(openxlsx)
wb <- createWorkbook()
for (i in seq_along(DEGs_MAST_celltype_number_sig)) {
  addWorksheet(wb, sheetName = names(DEGs_MAST_celltype_number_sig)[i])
  writeData(wb, sheet = i, x = DEGs_MAST_celltype_number_sig[[i]], rowNames = TRUE)
}
saveWorkbook(wb, file = "file/DEGs_MAST_celltype_number.xlsx", overwrite = TRUE)

DEGs_list <- list()
for (i in 1:35) {
  DEGs_list[[i]] <- rownames(DEGs_MAST_celltype_number[[i]][which(DEGs_MAST_celltype_number[[i]]$p_val_adj < 0.05 & abs(DEGs_MAST_celltype_number[[i]]$avg_log2FC) > 0.585), ])
}

DEGs_celltype_number <- c(unlist(DEGs_list))
DEGs_celltype_number <- DEGs_celltype_number[!duplicated(DEGs_celltype_number)]

# prepare the data for heatmap
heatmap_preparation_data_list <- list()
for (i in 1:35) {
  heatmap_preparation_data_list[[i]] <- data.frame(DEGs_MAST_celltype_number[[i]][which(DEGs_MAST_celltype_number[[i]]$p_val_adj < 0.05 & abs(DEGs_MAST_celltype_number[[i]]$avg_log2FC) > 0.585), c("avg_log2FC")],
                                                   gene = rownames(DEGs_MAST_celltype_number[[i]])[which(DEGs_MAST_celltype_number[[i]]$p_val_adj < 0.05 & abs(DEGs_MAST_celltype_number[[i]]$avg_log2FC) > 0.585)])
  colnames( heatmap_preparation_data_list[[i]]) <- c(celltype_number[i], "gene")
}


heatmap_preparation_data <- data.frame(DEGs_celltype_number, gene = DEGs_celltype_number)
for (i in 1:35) {
  heatmap_preparation_data <- full_join(heatmap_preparation_data, heatmap_preparation_data_list[[i]], by = "gene")
}
rownames(heatmap_preparation_data) <- heatmap_preparation_data$gene
heatmap_preparation_data <- heatmap_preparation_data[,-(1:2)]
heatmap_data <- tibble::rownames_to_column(heatmap_preparation_data,var = 'gene') %>% 
  pivot_longer( cols = colnames(heatmap_preparation_data),
                names_to = 'celltype_number',
                values_to = 'log2FC')
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data$gene <- factor(heatmap_data$gene, levels = rev(DEGs_celltype_number))
heatmap_data$celltype_number <- factor(heatmap_data$celltype_number, levels = celltype_number)

# plot
pdf('heatmap_celltype_number_DEGs.pdf' , height = 10, width = 10)
p <- ggplot(heatmap_data, aes(x = celltype_number, y = gene)) +
  geom_tile(aes(fill = log2FC), colour = NA) +
  scale_fill_gradient2(low = "blue4", mid = "white", high = "red") +
  scale_x_discrete("celltype_number", labels = c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32', 
                                                 "Inh_6", 
                                                 "Astro_1","Astro_5","Astro_26","Astro_31", 'Astro_34',
                                                 "Olig_0", "Olig_14","Olig_20","Olig_22","Olig_33",
                                                 "OPC_8", "OPC_19", 
                                                 "Micro_12","Micro_13","Micro_24","Micro_25","Micro_29", 
                                                 "Pericyte_7","Pericyte_17","Pericyte_18","Pericyte_23",
                                                 "Endo_2","Endo_9","Endo_11","Endo_21", 
                                                 "Macrophage_15")) +
  ylab(NULL) +
  ggtitle("celltype_number_heatmap") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),) +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1),
                                   vjust = 1,
                                   hjust = 1,
                                   angle = 45)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 
print(p)
dev.off()

#### Npy expression ####
##### fig 1f: violin plot of Npy expression level (inhibitory neuorns expressing Npy) #####
Idents(extinction.integrated_data) <- "celltype"
Inh_neuron <- subset(extinction.integrated_data, idents = "Inhibitory neuron")
Inh_expressing_neuron <- subset(Inh_neuron, subset = Npy > 0)
Idents(Inh_expressing_neuron) <- "stim"
pdf("violinplot_npy_expression_only_inh_expressing_neuron_split_by_stim.pdf", width = 5, height = 5)
VlnPlot(object = Inh_expressing_neuron, 
        group.by = "stim", 
        features = "Npy", 
        pt.size = 0.5, 
        cols = alpha(c("brown1", "royalblue3"), 0.8)) & 
  ylim(1, 6.5) & 
  NoLegend() &
  scale_x_discrete(labels = c("noExt", "ext"))
dev.off()

## statistic analysis (two-side T test)
Inh_Npy_expressing_neuron <- subset(Inh_neuron, subset = Npy > 0)
Idents(Inh_Npy_expressing_neuron) <- "stim"
Inh_Npy_expressing_neuron_extinction <- subset(Inh_Npy_expressing_neuron, idents = "extinction") %>% 
  GetAssayData(slot = "data") %>% 
  as.data.frame() %>%
  filter(row.names(.) %in% c("Npy"))
Inh_Npy_expressing_neuron_noExtinction <- subset(Inh_Npy_expressing_neuron, idents = "noExtinction") %>% 
  GetAssayData(slot = "data") %>% 
  as.data.frame() %>%
  filter(row.names(.) %in% c("Npy"))
t.test(x = as.numeric(Inh_Npy_expressing_neuron_extinction), 
       y = as.numeric(Inh_Npy_expressing_neuron_noExtinction),
       alternative = c("two.sided"), # "less", "greater"
       mu = 0, paired = FALSE, var.equal = T,
       conf.level = 0.95)

##### fig 1e: pie plot of Npy expressing percent #####
Idents(Inh_neuron) <- "stim"
npy_expressing_info <- data.frame(expressing_number = as.numeric(table(Idents(Inh_expressing_neuron))),
                                  all_number = as.numeric(table(Idents(Inh_neuron))),
                                  row.names = unique(Inh_neuron@active.ident)) %>% 
  mutate(expressing_percent = expressing_number/all_number) %>% 
  mutate(no_expressing_percent = 1-expressing_percent)

pdf("pie_plot_npy_expressing_percent.pdf", width = 7, height = 7)
# noExt
pie3D(t(npy_expressing_info)[3:4,1], explode = 0.04, labels = c("Npy+", "Npy-"), 
      col = c("hotpink1", "lavenderblush"), labelcex = 1, theta = pi/3, start = 0, shade = 0.9, height = 0.01, main = "noExtinction")
# ext
pie3D(t(npy_expressing_info)[3:4,2], explode = 0.04, labels = c("Npy+", "Npy-"), 
      col = c("hotpink1", "lavenderblush"), labelcex = 1, theta = pi/3, start = 0, shade = 0.9, height = 0.02, main = "extinction") 

dev.off()

#### Npy expressing percent and expression level in IEGs + Inh neuron ####
##### Extended fig 7a: Npy + percent in IEGs + Inh neuron (above the 90th percentile expression level of all IEGs) (Inh & Npy) #####
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
              "Id2", "Apold1") #128

Idents(extinction.integrated_data) <- "celltype"
Inh_neuron <- subset(extinction.integrated_data, idents = "Inhibitory neuron")
Inh_neuron_data <- GetAssayData(Inh_neuron, assay = "RNA", slot = "data")

# Replace zeros with NA to avoid including them in percentile calculations
IEGs_expression_level_90th_percentile_Inh_neuron <- Inh_neuron_data[IEG_list, ]
IEGs_expression_level_90th_percentile_Inh_neuron[IEGs_expression_level_90th_percentile_Inh_neuron == 0] <- NA

# Calculate the 90th percentile for each gene (row) in the dataset
IEGs_expression_level_90th_percentile_Inh_neuron <- apply(IEGs_expression_level_90th_percentile_Inh_neuron, 
                                                          1, 
                                                          function(x) quantile(x, probs = 0.90, na.rm = TRUE))
IEGs_data_Inh_neuron <- Inh_neuron_data[IEG_list,]
IEGs_data_Inh_neuron[which(IEGs_data_Inh_neuron <= IEGs_expression_level_90th_percentile_Inh_neuron)] <- NA
IEG_number_Inh_neuron <- apply(IEGs_data_Inh_neuron, 2, function(x) sum(x > 0, na.rm = T))
IEG_positive_cells_Inh_neuron <- names(IEG_number_Inh_neuron[IEG_number_Inh_neuron>0])
IEGs_positive_Inh_neuron <- extinction.integrated_data[,IEG_positive_cells_Inh_neuron]
Idents(IEGs_positive_Inh_neuron) <- "stim"
IEGs_positive_Inh_neuron_ext <- subset(IEGs_positive_Inh_neuron, idents = "extinction")
Idents(Inh_neuron) <- "stim"
IEGs_positive_Inh_neuron_noExt <- subset(IEGs_positive_Inh_neuron, idents = "noExtinction")

##### Extended fig 7b: Npy expression level in IEGs + Inh neuron (above the 90th percentile expression level of all IEGs) (Inh & Npy) (ext vs. noExt) ####
Npy_positive_IEGs_positive_Inh_neuron <- subset(IEGs_positive_Inh_neuron, subset = Npy > 0)
Npy_positive_IEGs_positive_Inh_neuron$stim <- factor(Npy_positive_IEGs_positive_Inh_neuron$stim, levels = c("noExtinction", "extinction"))

Npy_positive_IEGs_positive_Inh_neuron_ext <- subset(Npy_positive_IEGs_positive_Inh_neuron, idents = "extinction")
Npy_positive_IEGs_positive_Inh_neuron_noExt <- subset(Npy_positive_IEGs_positive_Inh_neuron, idents = "noExtinction")
Npy_expression_level_Inh_neuron_ext <- Npy_positive_IEGs_positive_Inh_neuron_ext["Npy",] %>% GetAssayData(., assay = "RNA", slot = "data") 
Npy_expression_level_Inh_neuron_noExt <-Npy_positive_IEGs_positive_Inh_neuron_noExt["Npy",] %>% GetAssayData(., assay = "RNA", slot = "data") 

t.test(Npy_expression_level_Inh_neuron_ext, Npy_expression_level_Inh_neuron_noExt)

pdf("plot/vlnplot_Npy_expression_level_90th_IEGs_positive_Inh_neuron_20250108.pdf")
VlnPlot(Npy_positive_IEGs_positive_Inh_neuron, 
        features = "Npy", 
        pt.size = 1, 
        group.by = "stim")
dev.off()


#### Inh subcluster ####
Idents(extinction.integrated_data) <- "celltype"
Inh <- subset(extinction.integrated_data, idents = "Inhibitory neuron")
table(Inh$stim)
Inh <- NormalizeData(object = Inh,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) %>% 
  ScaleData()
Inh <- FindVariableFeatures(Inh)
Inh <- RunPCA(Inh)
ElbowPlot(Inh)
Inh <- RunUMAP(Inh, dims = 1:11)

# Find neighbors and clusters
Inh <- FindNeighbors(Inh, dims = 1:11)

Inh <- FindClusters(Inh, resolution = 0.1)
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

##### Extended fig 8b: umap (subcluster) ####
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

##### findMarker #### 
Inh.FindAllMarkers <- FindAllMarkers(object = Inh,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)

Inh.FindAllMarkers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
save(Inh.FindAllMarkers, file = "rdata/Inh.res0.1.FindAllMarkers.rdata")

Inh.FindAllMarkers[which(Inh.FindAllMarkers$gene == "Npy"),]
# p_val avg_log2FC pct.1 pct.2    p_val_adj cluster gene
# Npy  2.177454e-12  0.6799553 0.339 0.166 4.213591e-08       0  Npy
# Npy1 2.586329e-07  1.2098196 0.331 0.192 5.004806e-03       2  Npy

##### Extended fig 8a:vlnplot (inh marker) ######
pdf("plot/vlnplot_inh_markers_inh_neuron_subclusters_resolution_0.1.pdf", width = 5, height = 10)
VlnPlot(Inh, 
        features = c("Slc32a1","Gad1","Gad2","Cck","Vip","Npy","Sst","Calb2","Lhx6"),
        group.by = "inh_number_renew", 
        cols = c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462'),
        pt.size = 0.1, ncol=1, combine = T)

dev.off()

###### Extended fig 8c: vlnplot (Npy) (split by stim) ######
Inh$Npy_expressionLevel <- GetAssayData(Inh, slot = "data")["Npy", ]

Inh@meta.data %>%
  group_by(inh_number_renew) %>% 
  group_map(~ t.test(Npy_expressionLevel ~ stim, data = .x), .keep = TRUE)

pdf("plot/vlnplot_npy_inh_neuron_subcluster_res0.1_split_by_stim.pdf", width = 10, height = 3)
VlnPlot(Inh, 
        features = "Npy",
        group.by = "inh_number_renew", 
        split.by = "stim", 
        cols = c("royalblue3", "brown1"),
        pt.size = 0.1, ncol=1, combine = T) 
dev.off()

#### Npy receptor expression ####
npy_receptor <- c("Npy1r","Npy2r","Npy5r") # "Npy4r" not found
npy_receptor_color <- c("mediumpurple1", "seagreen3", "dodgerblue")

##### Extended fig 25a: featureplot ("Npy1r","Npy2r","Npy5r") ####
Idents(extinction.integrated_data) <- "celltype"
neuron <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"))

pdf("FeaturePlot_npy_receptor.pdf", width=7, height=7)
p <- list()
for (i in 1:3) {
  p[[i]] <- FeaturePlot(object = neuron,
                        features = npy_receptor[i],
                        reduction = "umap",
                        cols = c("grey95", npy_receptor_color[i]),
                        combine = FALSE,
                        order = TRUE,
                        pt.size = 1)
}
print(p)
dev.off()
##### Extended fig 25b: violin plot (Npyr+ celltype or Npyr+ celltype_number) ####
all_expressing_cell <- list()
violinplot_celltype_plot <- list()
expressing_percent_celltype <- list()

neuron_only_expressing_cell <- list()
violinplot_celltype_number_plot <- list()
expressing_percent_celltype_number <- list()

celltype_abbre <- c("Exc", "Inh", "Astro", "Olig", "OPC", "Micro", "Pericyte", "Endo", "Macrophage")
for (a in 1:3) {
  ### celltype
  # violin plot
  Idents(extinction.integrated_data) <- "celltype"
  all_expressing_cell[[a]] <- subset(extinction.integrated_data, subset = !!as.name(npy_receptor[a]) > 0)
  Idents(all_expressing_cell[[a]]) <- "celltype"
  violinplot_celltype_plot[[a]] <- VlnPlot(object = all_expressing_cell[[a]], group.by = "celltype", 
                                           features = npy_receptor[a], pt.size = 0.5, 
                                           cols = rep(npy_receptor_color[a], 9)) & 
    ylim(1, 3.5) & 
    NoLegend() &
    scale_x_discrete(labels = c("Exc", "Inh", "Astro", "Olig", "OPC", "Micro", "Pericyte", "Endo", "Macrophage")) 
  
  # expressing percent
  expressing_percent_celltype[[a]] <- data.frame(expressing_cellnumber = as.numeric(table(Idents(all_expressing_cell[[a]]))),
                                                 all_cellnumber = as.numeric(table(Idents(extinction.integrated_data))),
                                                 row.names = celltype_abbre) %>% 
    mutate(expressing_percent = expressing_cellnumber/all_cellnumber)
  
  ### celltype_number
  # violin plot
  neuron_only_expressing_cell[[a]] <- subset(neuron, subset = !!as.name(npy_receptor[a]) > 0)
  Idents(neuron_only_expressing_cell[[a]]) <- "celltype_number"
  violinplot_celltype_number_plot[[a]] <- VlnPlot(object = neuron_only_expressing_cell[[a]], 
                                                  group.by = "celltype_number", 
                                                  features = npy_receptor[a], pt.size = 0.5, 
                                                  cols = alpha(rep(npy_receptor_color[a], 9), 0.8)) & 
    ylim(1, 3.5) & 
    NoLegend() 
  
  # expressing percent
  Idents(neuron) <- "celltype_number"
  expressing_percent_celltype_number[[a]] <- data.frame(expressing_cellnumber = as.numeric(table(Idents(neuron_only_expressing_cell[[a]]))),
                                                        all_cellnumber = as.numeric(table(Idents(neuron))),
                                                        row.names = celltype_abbre) %>% 
    mutate(expressing_percent = expressing_cellnumber/all_cellnumber)
  
}

pdf("violinPlot_celltype_npy_receptor_only_expressing_cells.pdf", width = 9, height = 9)
violinplot_celltype_plot <- grid.arrange(grobs=violinplot_celltype_plot,ncol=1) 
print(violinplot_celltype_plot)
dev.off()

names(expressing_percent_celltype) <- npy_receptor
write.xlsx(expressing_percent_celltype, file = "npy_receptor_expressing_percent_celltype.xlsx")

pdf("violinPlot_celltype_number_npy_receptor_only_expressing_cells.pdf", width = 9, height = 9)
violinplot_celltype_number_plot <- grid.arrange(grobs=violinplot_celltype_number_plot,ncol=1) 
print(violinplot_celltype_number_plot)
dev.off()

names(expressing_percent_celltype_number) <- npy_receptor
write.xlsx(expressing_percent_celltype_number, file = "npy_receptor_expressing_percent_celltype_number.xlsx")

##### fig 1b:  featureplot (Npy1r and Npy2r)  ####
Idents(extinction.integrated_data) <- "celltype"
neuron <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"))
only_npy1r_neuron <- subset(neuron, subset = Npy1r > 0 & Npy2r == 0)
only_npy2r_neuron <- subset(neuron, subset = Npy2r > 0 & Npy1r == 0)
npy1r_npy2r_neuron <- subset(neuron, subset = Npy1r > 0 & Npy2r > 0)
neuron@meta.data$Npyr <- ifelse(colnames(neuron) %in% colnames(only_npy1r_neuron), "npy1r",
                                ifelse(colnames(neuron) %in% colnames(only_npy2r_neuron), "npy2r",
                                       ifelse(colnames(neuron) %in% colnames(npy1r_npy2r_neuron), "npy1r_npy2r",
                                              "none")))
neuron@meta.data$Npyr <- factor(neuron@meta.data$Npyr, levels = c("none", "npy2r", "npy1r", "npy1r_npy2r"))

pdf("plot/umap_npy1r_npy2r_neuron.pdf", width=5, height=5)
DimPlot(neuron, reduction = "umap", group.by = "Npyr", label = TRUE, 
        cols= c("grey90", "seagreen3", "mediumpurple1","grey30"),
        pt.size = 0.2,
        order = T,
        label.size = 0
)

FeaturePlot(object = neuron,
            features = "Npy1r",
            reduction = "umap",
            col = c("grey90","mediumpurple1"),
            pt.size = 0.2,
            order = T,
            blend = F) + ggtitle("Npy1r")
FeaturePlot(object = neuron,
            features = "Npy2r",
            reduction = "umap",
            col = c("grey90", "seagreen3"),
            pt.size = 0.2,
            order = T) + ggtitle("Npy2r")
dev.off()


##### Extended fig 25d: Npyr no-overlap bar plot ####
Idents(extinction.integrated_data) <- "celltype"
neuron <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"))
no_npy1r_npy2r_neuron <- subset(neuron, subset = Npy1r == 0 & Npy2r == 0)
only_npy1r_neuron <- subset(neuron, subset = Npy1r > 0 & Npy2r == 0)
only_npy2r_neuron <- subset(neuron, subset = Npy2r > 0 & Npy1r == 0)
npy1r_npy2r_neuron <- subset(neuron, subset = Npy1r > 0 & Npy2r > 0)
npy1r_npy2r_expressing_cellNumber_neuron_mat <- matrix(c(7383, 432, 487, 36), nrow = 2, byrow = TRUE)
colnames(npy1r_npy2r_expressing_cellNumber_neuron_mat) <- c("Npy1r_Negative", "Npy1r_Positive")
rownames(npy1r_npy2r_expressing_cellNumber_neuron_mat)<- c("Npy2r_Negative", "Npy2r_Positive")

chisq.test(npy1r_npy2r_expressing_cellNumber_neuron_mat)

npy1r_npy2r_expressing_cellNumber_neuron_mat_long <- as.data.frame(npy1r_npy2r_expressing_cellNumber_neuron_mat) %>%
  rownames_to_column("Npy2r_status") %>%
  pivot_longer(cols = -Npy2r_status, 
               names_to = "Npy1r_status", 
               values_to = "Cell_Count")

pdf("plot/barplot_Npy1r_npy2r_no_overlap_neuron.pdf")
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




#### fig 4g: DEGs between Npy1r and Npy2r expressing neurons ####
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
diff_Npy1r_Npy2r_nofilter <- list()
avg_Npy1r_Npy2r <- list()

for (i in 1:2){
  Npy1r_only_neuron_list[[i]] <- subset(subset(neuron_list[[i]], subset = Npy1r >0), subset = Npy2r <=0)
  Npy2r_only_neuron_list[[i]] <- subset(subset(neuron_list[[i]], subset = Npy2r >0), subset = Npy1r <=0)
  Idents(Npy1r_only_neuron_list[[i]]) <- "stim"
  Idents(Npy2r_only_neuron_list[[i]]) <- "stim"
  
  diff_Npy1r_Npy2r_filter[[i]] <- FindMarkers(neuron_list[[i]], 
                                              ident.1 = colnames(Npy1r_only_neuron_list[[i]]),
                                              ident.2 = colnames(Npy2r_only_neuron_list[[i]]),
                                              t.use ="MAST", 
                                              verbose = FALSE) 
  
  diff_Npy1r_Npy2r_nofilter[[i]] <- FindMarkers(neuron_list[[i]], 
                                                ident.1 = colnames(Npy1r_only_neuron_list[[i]]),
                                                ident.2 = colnames(Npy2r_only_neuron_list[[i]]),
                                                t.use ="MAST", 
                                                verbose = FALSE,
                                                logfc.threshold = 0,min.pct = 0) 
  
  avg_Npy1r_Npy2r[[i]] <- as.data.frame(log1p(AverageExpression(Npy1r_only_neuron_list[[i]] , verbose = FALSE)$RNA))
  avg_Npy1r_Npy2r[[i]]$Npy2r_all <- as.data.frame(log1p(AverageExpression(Npy2r_only_neuron_list[[i]] , verbose = FALSE)$RNA))[,1]
  avg_Npy1r_Npy2r[[i]]$gene <- rownames(avg_Npy1r_Npy2r[[i]])
  avg_Npy1r_Npy2r[[i]]$pct.1 <- diff_Npy1r_Npy2r_filter[[i]][match(rownames(avg_Npy1r_Npy2r[[i]]), rownames(diff_Npy1r_Npy2r_filter[[i]])), ]$pct.1
  avg_Npy1r_Npy2r[[i]]$pct.2 <- diff_Npy1r_Npy2r_filter[[i]][match(rownames(avg_Npy1r_Npy2r[[i]]), rownames(diff_Npy1r_Npy2r_filter[[i]])), ]$pct.2
  avg_Npy1r_Npy2r[[i]]$p_adj_MAST <- diff_Npy1r_Npy2r_filter[[i]][match(rownames(avg_Npy1r_Npy2r[[i]]), rownames( diff_Npy1r_Npy2r_filter[[i]])), ]$p_val_adj
  avg_Npy1r_Npy2r[[i]]$avg_logFC_MAST <- diff_Npy1r_Npy2r_filter[[i]][match(rownames(avg_Npy1r_Npy2r[[i]]), rownames(diff_Npy1r_Npy2r_filter[[i]])), ]$avg_log2FC
  avg_Npy1r_Npy2r[[i]]$p_adj_MAST_no_filter <- diff_Npy1r_Npy2r_nofilter[[i]][match(rownames(avg_Npy1r_Npy2r[[i]]), rownames( diff_Npy1r_Npy2r_nofilter[[i]])), ]$p_val_adj
  avg_Npy1r_Npy2r[[i]]$avg_logFC_MAST_no_filter <- diff_Npy1r_Npy2r_nofilter[[i]][match(rownames(avg_Npy1r_Npy2r[[i]]), rownames(diff_Npy1r_Npy2r_nofilter[[i]])), ]$avg_log2FC
}

expression_plot <- list()
for (i in 1:2){
  expression_plot[[i]] <- ggplot(data = avg_Npy1r_Npy2r[[i]], aes(x = all, y = Npy2r_all, label = gene, colour = ((p_adj_MAST_no_filter < 0.05 & abs(avg_logFC_MAST_no_filter) > 0.25 )==T))) +
    scale_color_manual(values = c("grey50", "brown1")) +
    geom_point() +
    geom_text_repel(data=subset(avg_Npy1r_Npy2r[[i]], p_adj_MAST_no_filter < 0.05 & abs(avg_logFC_MAST_no_filter) > 0.25 ), colour='red') + 
    xlab("Npy1r +") +
    ylab("Npy2r +") + 
    xlim(0,5) +
    ylim(0,5) +
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

pdf("Nrp1r_Nrp2r_expression_scatter_plot.pdf", width = 5, height = 5)
for (i in 1:2){
  print(expression_plot[[i]])
}
dev.off()

#### Extended fig 27: Npy1r+ vs. Npy2r+ neuron activity (IEGs_sum expression) in ext or noExt conditions ####
neuron_Npy1r_pos <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"), subset = Npy1r > 0)
neuron_Npy1r_pos$Npy1r_Npy2r <- "Npy1r"
neuron_Npy2r_pos <- subset(extinction.integrated_data, idents = c("Excitatory neuron", "Inhibitory neuron"), subset = Npy2r > 0)
neuron_Npy2r_pos$Npy1r_Npy2r <- "Npy2r"
neuron_Npy1r_Npy2r_pos <- merge(neuron_Npy1r_pos, neuron_Npy2r_pos)
neuron_Npy1r_Npy2r_pos$IEGs_sum <- colSums(GetAssayData(neuron_Npy1r_Npy2r_pos, slot = "data")[IEG_list, ])

Idents(neuron_Npy1r_Npy2r_pos) <- "stim"
neuron_Npy1r_Npy2r_pos$stim <- factor(neuron_Npy1r_Npy2r_pos$stim, levels = c("noExtinction", "extinction"))
neuron_Npy1r_Npy2r_pos@meta.data %>%
  group_by(stim) %>% 
  group_map(~ t.test(IEGs_sum ~ Npy1r_Npy2r, data = .x), .keep = TRUE)

pdf("plot/vlnplot_IEGs_sum_Npy1r_neuron_vs_Npy2r_neuron_in_noExt_or_ext.pdf", width = 10, height = 5)
VlnPlot(neuron_Npy1r_Npy2r_pos, features = "IEGs_sum", split.by = "Npy1r_Npy2r", group.by = "stim")
dev.off()

#### Different expression levels of IEGs across different cell types and clusters of No Ext and Ext mice #### 
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
              "Id2", "Apold1") #128
##### Extended fig 5a: heatmap: IEGs log2FC (celltype) ####
load("rdata/DEGs_MAST_no_filter.rdata")
celltype <- c("Exc", "Inh", "Astro", "Olig",  "OPC", "Micro", "Pericyte", "Endo", "Macrophage")

# prepare heatmap data
all_avg_log2FC_IEGs_celltype <- data.frame(row.names = IEG_list)
for (i in 1:9) {
  all_avg_log2FC_IEGs_celltype[celltype[i]] <- DEGs_MAST_no_filter[[i]][IEG_list,"avg_log2FC"]                                  
}

heatmap_data_all_avg_log2FC_IEGs_celltype <- tibble::rownames_to_column(all_avg_log2FC_IEGs_celltype, var = 'gene') %>% 
  pivot_longer( cols = colnames(all_avg_log2FC_IEGs_celltype),
                names_to = 'celltype',
                values_to = 'log2FC')

heatmap_data_all_avg_log2FC_IEGs_celltype$gene <- factor(heatmap_data_all_avg_log2FC_IEGs_celltype$gene, levels = rev(IEG_list))
heatmap_data_all_avg_log2FC_IEGs_celltype$celltype <- factor(heatmap_data_all_avg_log2FC_IEGs_celltype$celltype, levels = celltype)

my_palette <- c(colorRampPalette(c("blue4", "white"))(148), "white",colorRampPalette(c("white", "red"))(106))
zero_color_index <- floor((0 - min(all_avg_log2FC_IEGs_celltype)) / 
                            (max(all_avg_log2FC_IEGs_celltype) - min(all_avg_log2FC_IEGs_celltype)) * 255) + 1
my_palette[zero_color_index] <- "white"

pdf('plot/heatmap_IEGs_log2FC_celltype.pdf' , height = 8, width = 5)
pheatmap(all_avg_log2FC_IEGs_celltype,
         cluster_col = FALSE, 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_palette,
         border_color = NA,
         legend = TRUE,
         legend_breaks = c(-0.5, 0, 0.5),
         legend_labels = c("-0.5", "0", "0.5"))
dev.off()

##### Extended fig 5b: IEGs log2FC (Exc subcluster) ####
load("rdata/DEGs_MAST_celltype_number_no_filter.rdata")

Exc_cluster_number <- c("Exc_3","Exc_4","Exc_10","Exc_16","Exc_27","Exc_28","Exc_30", 'Exc_32')

# prepare heatmap data
all_avg_log2FC_IEGs_Exc <- data.frame(row.names = IEG_list)
for (i in 1:8) {
  all_avg_log2FC_IEGs_Exc[Exc_cluster_number[i]] <- DEGs_MAST_celltype_number_no_filter[[i]][IEG_list,"avg_log2FC"]                                  
}

heatmap_data_all_avg_log2FC_IEGs_Exc <- tibble::rownames_to_column(all_avg_log2FC_IEGs_Exc, var = 'gene') %>% 
  pivot_longer( cols = colnames(all_avg_log2FC_IEGs_Exc),
                names_to = 'Exc',
                values_to = 'log2FC')# heatmap_data_all_avg_log2FC_IEGs_Exc[is.na(heatmap_data_all_avg_log2FC_IEGs_Exc)] <- 0
heatmap_data_all_avg_log2FC_IEGs_Exc$gene <- factor(heatmap_data_all_avg_log2FC_IEGs_Exc$gene, levels = rev(IEG_list))
heatmap_data_all_avg_log2FC_IEGs_Exc$Exc <- factor(heatmap_data_all_avg_log2FC_IEGs_Exc$Exc, levels = Exc_cluster_number)

my_palette <- c(colorRampPalette(c("blue4", "white"))(94), "white",colorRampPalette(c("white", "red"))(160))

zero_color_index <- floor((0 - min(all_avg_log2FC_IEGs_Exc)) / 
                            (max(all_avg_log2FC_IEGs_Exc) - min(all_avg_log2FC_IEGs_Exc)) * 255) + 1
my_palette[zero_color_index] <- "white"

pdf('plot/heatmap_IEGs_log2FC_Exc_subclusters.pdf' , height = 8, width = 5)
pheatmap(all_avg_log2FC_IEGs_Exc,
         cluster_col = FALSE, 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_palette,
         border_color = NA,
         legend = TRUE,
         legend_breaks = c(-0.5, 0, 0.5, 1),
         legend_labels = c("-0.5", "0", "0.5","1"))
dev.off()


