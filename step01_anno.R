setwd("/home/disk1/wgs")
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggsci)
library(sctransform)
library(glmGamPoi)
library(pbapply)
library(qs)
library(patchwork)
options(future.globals.maxSize = 2048 * 1024^5)
source("useful_functions.R")
source("/home/disk1/zzm/reference/rscripts tools/Fun_run seurat.R")

##################
### clustering ###
##################
seu <- qread("batch2/data_SCT_integrated.qs")
seu[["sample"]] <- seu$orig.ident

#QC plot
#analyze and plot for Figure S1A, quality control metrics by sample
DefaultAssay(seu) <- "RNA"
seu[["percent.mt"]] <- NULL
seu[["percent_mt"]] <- PercentageFeatureSet(seu, "^mt-")
seu[["percent_ribo"]] <- PercentageFeatureSet(seu, "^Rp[sl]")
seu[["percent_hb"]] <- PercentageFeatureSet(seu, "^Hb[^(p|e|s)]")
seu <- subset(seu,subset = percent_hb<1)
plot_qc_results(seu,islist = F) -> p; plot(p)
savePDF(p, "01_clustering/QC_plot.pdf", width = 12, height = 2)
write.csv(qc_process(seu,islist = F), "01_clustering/QC_info.csv", row.names = F)

#clustering
#analyze for Figure 1A, cluster and annotate major cell types
DefaultAssay(seu) <- "integrated"
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
seu <- FindClusters(seu, resolution = 0.1)
DimPlot(seu,label = T)
DimPlot(seu,group.by = "celltype")
Idents(seu) <- seu$seurat_clusters
seu <- FindSubCluster(seu, c("0"), graph.name = "integrated_snn", resolution = 0.1)
DimPlot(seu,label = T, group.by = "sub.cluster")
Idents(seu) <- seu$sub.cluster
seu <- FindSubCluster(seu, c("0_0"), graph.name = "integrated_snn", resolution = 0.05)
DimPlot(seu,label = T, group.by = "sub.cluster")
Idents(seu) <- seu$sub.cluster
seu <- FindSubCluster(seu, c("3"), graph.name = "integrated_snn", resolution = 0.05)
DimPlot(seu,label = T, group.by = "sub.cluster")
Idents(seu) <- seu$sub.cluster
seu <- FindSubCluster(seu, c("3_0"), graph.name = "integrated_snn", resolution = 0.1)
DimPlot(seu,label = T, group.by = "sub.cluster")
Idents(seu) <- seu$sub.cluster
seu <- FindSubCluster(seu, c("8"), graph.name = "integrated_snn", resolution = 0.05)
DimPlot(seu,label = T, group.by = "sub.cluster")
seu@meta.data$celltype <- "unknown"
seu@meta.data[seu$sub.cluster%in%c("2","6","0_0_0","12","0_5","11","0_4","0_6","3_0_0"),"celltype"] <- "c01:GABAergic Neurons"
seu@meta.data[seu$sub.cluster%in%c("3_1","3_0_1","3_0_2","1","13","0_1","0_0_1","0_3","13","0_2","3_2"),"celltype"] <- "c02:Glutamatergic Neurons"
seu@meta.data[seu$sub.cluster%in%c("8_1"),"celltype"] <- "c03:RGCs"
seu@meta.data[seu$sub.cluster%in%c("4"),"celltype"] <- "c04:Astrocytes"
seu@meta.data[seu$sub.cluster%in%c("9"),"celltype"] <- "c05:OPCs"
seu@meta.data[seu$sub.cluster%in%c("7"),"celltype"] <- "c06:Oligdendrocytes"
seu@meta.data[seu$sub.cluster%in%c("10"),"celltype"] <- "c07:Microglia"
seu@meta.data[seu$sub.cluster%in%c("5"),"celltype"] <- "c08:Tanycytes"
seu@meta.data[seu$sub.cluster%in%c("8_0"),"celltype"] <- "c09:ParsTuber"
seu@meta.data[seu$sub.cluster%in%c("8_2"),"celltype"] <- "c10:Vascular cells"
seu@meta.data[seu$sub.cluster%in%c("8_3"),"celltype"] <- "c11:Pituitary endocrine cells"
seu$cluster <- gsub(".*:","",seu$celltype)
seu$cluster_id <- gsub(":.*","",seu$celltype)
seu <- RunUMAP(seu, dims = 1:30, spread = 1, verbose = FALSE)
DimPlot(seu,group.by = "celltype",label=T)
seu <- PrepSCTFindMarkers(seu)
qsave(seu,"01_clustering/s01_Integration Data.qs") 
seu <- qread("01_clustering/s01_Integration Data.qs")

#07.UMAP plot
#plot for Figure 1A, UMAP of major cell types
umap <- seu@reductions$umap@cell.embeddings %>% as.data.frame()
umap$celltype <- seu$celltype
ggplot(umap)+
  geom_scattermore(aes(UMAP_1,UMAP_2,color=celltype),pointsize = 1.5,pixels = c(1024, 1024))+
  scale_color_manual(values = c(paletteer_d("ggsci::default_jco")[c(3,4,1,6,8,2,9,7,5,10)],"#cc99ff"))+
  theme_bw()+
  theme(
    legend.key.width = unit(0.5,"cm"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5,size=10))+
  labs(x=NULL,y=NULL)-> p; plot(p)
savePDF(p, "01_clustering/p01_all_cells_UMAP.pdf", width = 7.6, height = 5)


#08.marker plot
#plot for Figure S1C, UMAP of cell type marker expression
gene <- c("Tubb3","Slc32a1","Slc17a6","Aldh1l1","Mbp","Cga","Col23a1","Vim","Pdgfra","Cx3cr1","Pdgfrb","Tbx19")
p <- lapply(gene,function(x){
  require(paletteer)
  require(scattermore)
  obj <- seu
  umap <- obj@reductions$umap@cell.embeddings %>% as.data.frame()
  umap$expr <- obj@assays$SCT@data[x,]
  umap$group <- obj$sample_group
  umap <- umap[order(umap$expr),]
  umap$normexpr <- umap$expr/quantile(umap$expr[umap$expr>0],0.9)
  umap$normexpr[umap$normexpr>1] <- 1
  ggplot(umap, aes(UMAP_1,UMAP_2))+
    geom_scattermore(data = umap[umap$expr==0,],
                     pixels = c(1024, 1024),pointsize = 2,color="#d0d0d0")+
    geom_scattermore(data = umap[umap$expr>0,],aes(color=normexpr),
                     pixels = c(1024, 1024),pointsize = 2)+
    geom_scattermore(aes(UMAP_1,UMAP_2,color=normexpr),
                     pixels = c(1024, 1024),pointsize = 2,alpha = 0)+
    scale_color_gradientn(values = c(0,0.05,1),
                          colours = c("#C0C0C0",paletteer_c("grDevices::Blues",30)[c(15,1)]))+
    theme_minimal()+
    theme(
      legend.key.width = unit(0.5,"cm"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5,size=10))+
    labs(x=NULL,y=NULL,title = x)
})
wrap_plots(p,nrow = 3)+plot_layout(guides = "collect") -> plist; plot(plist)
savePDF(plist, "01_clustering/p01_all_cells_UMAP_marker.pdf", width = 15, height = 9)

