# load into your session
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(SingleR)
library(celldex)

setwd("/Users/angiepoon528/Desktop/SRF/sgRNA-seq workshop/workshop2")

#load data 
pbmc.data <- Read10X_h5(filename="./SC3_v3_NextGem_DI_PBMC_10K_filtered_feature_bc_matrix.h5")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, verbose = FALSE)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca",balanced=TRUE)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Jackstraw plot 
#Jackstraw plot estimates the significance of the structure captured by each principal component
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)

ElbowPlot(pbmc,ndims=50)
###Select dim 1:20 

pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)

### Clustering 
##first computes K nearest neighbors of each cell in PCA space. 
#Jaccard similarity between cells is calculated based on shared neighborhood
pbmc <- FindNeighbors(pbmc, dims = 1:20)

#resolution 
#Small number will generate admixed clusters while large number would give too many clusters which are unlikely to be meaningful.
#trying 3 resolution 
pbmc <- FindClusters(object = pbmc,  resolution = c(0.5, 1, 1.5),  dims.use = 1:10,  save.SNN = TRUE)
saveRDS(pbmc, file = "pbmc_tutorial.rds")
clustree(pbmc)
###3 rows: 3 resolution 
###dot size: # of cells in that cluster 
###Arrows :  show how the cluster assignment changes with increasing resolution
###Some clusters may split into two (or more). That suggests increasing resolution. 
###But if you see a lot of back and forth jumping between clusters, it indicates less stability. 
### better to over cluster, examine genes and then manually choose which clusters to merge.

##resolution =0.5 
Idents(pbmc) <- pbmc$RNA_snn_res.0.5
DimPlot(pbmc, reduction = "umap", label=TRUE)

### Finding cluster enriched marker 
##ROC: classification power of each marker (0 implies random, 1 implies perfect )
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster.all.markers0.5 <- FindAllMarkers(pbmc, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE, min.pct = 0.25) #All clustera 

###heatmap making of highly enriched gene in each cluster 
cluster.all.markers0.5  %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

DoHeatmap(pbmc, features = top5$gene) + NoLegend()

#manual picking of common markers 
#Many of these are from the cluster.all.markers0.5
markers.to.plot <- c("CD3D", "HSPH1", "SELL", "CD14", "LYZ", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "FCGR3A", "MS4A7", "S100A9", "HLA-DQA1","GPR183", "PPBP", "GNG11", "TSPAN13", "IL3RA", "FCER1A", "CST3", "S100A12")
DotPlot(pbmc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +RotatedAxis()
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "GNLY"))
VlnPlot(pbmc, features = c("FCGR3A", "MS4A7"))
VlnPlot(pbmc, features = c("PPBP"))
VlnPlot(pbmc, features = c("FCER1A", "CST3"))   
VlnPlot(pbmc, features = c("CD8A", "CD8B", "CD3D"))
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",  "CD8A"))

###Reference based cluster annotation
##change seurat formt to compatilble format 
sce <- GetAssayData(object = pbmc, assay = "RNA", slot = "data")
##reference: compares different immune cell types using sorted cells.
refMonaco <- MonacoImmuneData()
#information at two levels: main cell types and finer resolution information
prediction_Monaco_main <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc), labels=refMonaco$label.main)
prediction_Monaco_fine <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc), labels=refMonaco$label.fine)
prediction_Monaco_fine <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc), labels=refMonaco$label.fine)
predicted_Monaco <- data.frame(cluster=sort(unique(Idents(pbmc))), Monaco_main= prediction_Monaco_main$labels, Monaco_fine= prediction_Monaco_fine$labels)

##higher resolution 
pbmc1 <- pbmc
Idents(pbmc1) <- pbmc1$RNA_snn_res.1.5
DimPlot(pbmc1, reduction = "umap", label=TRUE)
#C
cluster.all.markers0.51 <- FindAllMarkers(pbmc1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE, min.pct = 0.25) #All clustera 
cluster.all.markers0.51  %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(pbmc1, features = top5$gene) + NoLegend()
prediction_Monaco_main1 <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc1), labels=refMonaco$label.main)
prediction_Monaco_fine1 <- SingleR(test=sce, ref=refMonaco, clusters=Idents(pbmc1), labels=refMonaco$label.fine)
predicted_Monaco1 <- data.frame(cluster=sort(unique(Idents(pbmc1))), Monaco_main= prediction_Monaco_main1$labels, Monaco_fine= prediction_Monaco_fine1$labels)
