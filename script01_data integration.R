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
source("NEW/useful_functions.R")
source("/home/disk1/zzm/reference/rscripts tools/Fun_run seurat.R")


#01.read counts data
files <- data.frame(sample = list.files(".",pattern = "^A"),
                    files = paste(list.files(".",pattern = "^A"),"/outs/filtered_feature_bc_matrix/",sep = ""))
files <- files[c(1:2,5:6),]
mtx <- lapply(1:4, function(x){
  s <- files$sample[x]
  f <- files$files[x]
  mtx <- Seurat::Read10X(f)
})


#02.seurat object and reducing
#数据质控1：所有细胞（sample level，颜色按分组分类）nFeatureRNA nCountRNA percent_mt percent_ribo统计（QC_info.csv）以及图（QC_plot.pdf）
seu <- lapply(1:4,function(x){
  s <- files$sample[x]
  m <- mtx[[x]]
  colnames(m) <- paste(s,colnames(m),sep="_")
  seu <- CreateSeuratObject(m,project = "snRNAseq")
  seu[["sample"]] <- s
  if(s %in% c("AM3N378","AM3N379"))  seu[["sample_group"]] <- "NC" else seu[["sample_group"]] <- "T"
  seu[["percent_mt"]] <- PercentageFeatureSet(seu, "^mt-")
  seu[["percent_ribo"]] <- PercentageFeatureSet(seu, "^Rp[sl]")
  seu[["percent_hb"]] <- PercentageFeatureSet(seu, "^Hb[^(p|e|s)]")
  return(seu)
})
write.csv(qc_process(seu,islist = T), "NEW/QC_info.csv", row.names = F)
plot_qc_results(seu) -> p; plot(p)
savePDF(p, "NEW/QC_plot.pdf", width = 10, height = 2)
seu <- lapply(seu,function(x){
  x <- subset(x,subset = nFeature_RNA < 7500 & percent.mt < 10)
  obj <- SCTransform(x, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30) %>% FindClusters(resolution = 0.4)
})


#03.doublet finding
#使用DoubletFinder去除潜在双细胞的污染
library(DoubletFinder)
seu <- lapply(1:4, function(x){
  annotations <- seu[[x]]$seurat_clusters
  sweep.res.list <- paramSweep_v3(seu[[x]],PCs=1:30,sct=T)
  sweep.stats <- summarizeSweep(sweep.res.list,GT=F)
  bcmvn <- find.pK(sweep.stats)
  pk <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.05*nrow(seu[[x]]@meta.data))
  nExp_adj <- round(nExp_poi*(1-homotypic.prop))
  
  obj <- doubletFinder_v3(seu[[x]],PCs=1:30,pN = 0.25,pK = pk,
                          nExp = nExp_adj,reuse.pANN=F,sct=T)
})
doul <- data.frame(row.names = sapply(seu,function(x){rownames(x@meta.data)})%>%unlist,
                   doubletfinder = sapply(seu,function(x){x@meta.data[,10]})%>%unlist)
qsave(doul,"NEW/s01_doubletfinder results.qs") 
qsave(seu,"NEW/s02_SCTdata sample-level w doubletfinder.qs") 



#04.Integration
#批次效应去除，降维
seu <- qread("NEW/s02_SCTdata sample-level w doubletfinder.qs")
features <- SelectIntegrationFeatures(object.list = seu, nfeatures = 3000)  
seu <- PrepSCTIntegration(object.list = seu, anchor.features = features)
seu <- FindIntegrationAnchors(object.list = seu, normalization.method = "SCT",
                              reference = c(1,3),reduction="rpca",k.anchor = 5,
                              anchor.features = features)
seu <- IntegrateData(anchorset = seu, normalization.method = "SCT")
seu@meta.data <- seu@meta.data[,1:6]
seu <- AddMetaData(seu,qread("NEW/s01_doubletfinder results.qs"))
seu@meta.data <- seu@meta.data %>%
  mutate(sample_group = case_when(orig.ident%in%c("AM3N378","AM3N379") ~ "NC", TRUE ~ "T"))
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:15,spread = 2)
Groupplot(seu, feat = "sample_group", size = 0.3) + scale_color_manual(values = c("#FFB900", "#5773CC")) ->p; plot(p)
savePDF(p, "NEW/QC_Integration.pdf", width = 5, height = 4)
qsave(seu,"NEW/s03_Integration Data k5.qs") 
#04.2 过滤后质控指标
qc <- qc_process(seu, islist = F)
write.csv(qc, "NEW/QC_info2.csv", row.names = F)


#06.cell assignment
seu <- qread("NEW/s03_Integration Data k5.qs")
seu@meta.data <- seu@meta.data[,c(1:8)]
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:15)
seu <- FindClusters(seu, resolution = 1)
DimPlot(seu,label = T)
seu@meta.data$celltype <- "unknown"
seu@meta.data[seu$seurat_clusters%in%c(17,3,7,21,15,14,11,24,0,16),"celltype"] <- "c01:GABAergic Neurons"
seu@meta.data[seu$seurat_clusters%in%c(2,5,6,8,1,13,22),"celltype"] <- "c02:Glutamatergic Neurons"
seu@meta.data[seu$seurat_clusters%in%c(4),"celltype"] <- "c03:Astrocytes"
seu@meta.data[seu$seurat_clusters%in%c(12),"celltype"] <- "c04:Oligdendrocytes"
seu@meta.data[seu$seurat_clusters%in%c(20,29),"celltype"] <- "c05:ParsTuber"
seu@meta.data[seu$seurat_clusters%in%c(9,10),"celltype"] <- "c06:Tanycytes"
seu@meta.data[seu$seurat_clusters%in%c(25),"celltype"] <- "c07:RGCs"
seu@meta.data[seu$seurat_clusters%in%c(18,31),"celltype"] <- "c08:OPCs"
seu@meta.data[seu$seurat_clusters%in%c(19),"celltype"] <- "c09:Microglia"
seu@meta.data[seu$seurat_clusters%in%c(23),"celltype"] <- "c10:Vascular cells"
DimPlot(seu,group.by = "celltype",label=T)
seu <- seu[,!seu$celltype%in%"unknown"]
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:15,spread = 2)
seu@meta.data <- seu@meta.data[,c(1:8)]
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:15)
seu <- FindClusters(seu, resolution = 1)
DimPlot(seu,label = T)
seu@meta.data$celltype <- "unknown"
seu@meta.data[seu$seurat_clusters%in%c(2,16,17,22,12,4,0,14,26,10),"celltype"] <- "c01:GABAergic Neurons"
seu@meta.data[seu$seurat_clusters%in%c(3,5,6,1,8,13,23),"celltype"] <- "c02:Glutamatergic Neurons"
seu@meta.data[seu$seurat_clusters%in%c(15,20),"celltype"] <- "c03:Astrocytes"
seu@meta.data[seu$seurat_clusters%in%c(11,28),"celltype"] <- "c04:Oligdendrocytes"
seu@meta.data[seu$seurat_clusters%in%c(18),"celltype"] <- "c05:OPCs"
seu@meta.data[seu$seurat_clusters%in%c(7,9),"celltype"] <- "c06:Tanycytes"
seu@meta.data[seu$seurat_clusters%in%c(25),"celltype"] <- "c07:RGCs"
seu@meta.data[seu$seurat_clusters%in%c(19),"celltype"] <- "c08:Microglia"
seu@meta.data[seu$seurat_clusters%in%c(21,27),"celltype"] <- "c09:ParsTuber"
seu@meta.data[seu$seurat_clusters%in%c(24),"celltype"] <- "c10:Vascular cells"
seu$cluster <- gsub(".*:","",seu$celltype)
seu$cluster_id <- gsub(":.*","",seu$celltype)
DimPlot(seu,group.by = "celltype",label=T)
Idents(seu) <- as.factor(seu$celltype)
seu <- PrepSCTFindMarkers(seu)
qsave(seu,"NEW/s03_Integration Data mapped k5.qs") 
seu <- qread("NEW/s03_Integration Data mapped k5.qs")
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)

#07.UMAP plot
umap <- seu@reductions$umap@cell.embeddings %>% as.data.frame()
umap$celltype <- seu$celltype
ggplot(umap)+
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color=celltype),raster.dpi = 600, size=0.1)+
  scale_color_manual(values = c(paletteer_d("ggsci::default_jco")[c(3,4,1,6,8,2,9,7,5,10)],"#cc99ff"))+
  theme_bw()+
  theme(
    legend.key.width = unit(0.5,"cm"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5,size=10))+
  labs(x=NULL,y=NULL)
#landscape 12*7

#08.marker plot
gene <- c("Tubb3","Slc32a1","Slc17a6","Aldh1l1","Mbp","Pdgfra","Col23a1","Vim","Cx3cr1","Kcnk3","Epcam","Pdgfrb")
p <- lapply(gene,function(x){
  require(paletteer)
  require(scattermore)
  obj <- seu
  umap <- obj@reductions$umap@cell.embeddings %>% as.data.frame()
  umap$expr <- obj@assays$SCT@data[x,]
  umap$group <- obj$sample_group
  umap <- umap[order(umap$expr),]
  umap$normexpr <- umap$expr/quantile(umap$expr[umap$expr>0],0.75)
  umap$normexpr[umap$normexpr>1] <- 1
  ggplot(umap, aes(UMAP_1,UMAP_2))+
    geom_scattermore(data = umap[umap$expr==0,],
                     pixels = c(1024, 1024),pointsize = 4,color="#C0C0C0")+
    geom_scattermore(data = umap[umap$expr>0,],aes(color=normexpr),
                     pixels = c(1024, 1024),pointsize = 4)+
    geom_scattermore(aes(UMAP_1,UMAP_2,color=normexpr),
                     pixels = c(1024, 1024),pointsize = 4,alpha = 0)+
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
wrap_plots(p,nrow = 3)+plot_layout(guides = "collect")

#08. 数据质控2：各群细胞（sample level）nFeatureRNA nCountRNA percent_mt percent_ribo统计（QC_info.csv）
qc <- lapply(unique(seu$celltype), function(x){
  seux <- subset(seu, subset = celltype %in% x)
  qc_process(seux, islist = F) %>% 
    mutate(celltype = x)
}) %>% purrr::reduce(rbind)
write.csv(qc, "NEW/QC_info3.csv", row.names = F)


