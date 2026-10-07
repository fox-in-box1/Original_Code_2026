setwd("/home/disk1/wgs")
library(Seurat)
library(dplyr)
library(ggplot2)
library(pbapply)
library(qs)
library(ggpointdensity)
library(paletteer)
library(scattermore)
source("useful_functions.R")
source("/home/disk1/zzm/reference/rscripts tools/Fun_run seurat.R")
options(future.globals.maxSize = 200 * 1024^3)

#01.Literature markers
#analyze, load literature-derived neuronal marker sets
seu <- qread("01_clustering/s01_Integration Data.qs")
feat <- readxl::read_xlsx("ref/Literature Markers.xlsx") %>% as.data.frame()
name <- feat$name
feat <- lapply(1:nrow(feat), function(x){
  gene <- feat$genes[x] %>% stringr::str_split(";") %>% unlist %>% .[.%in%rownames(seu@assays$RNA)]
})
names(feat) <- name



#################################
## integrate data & annotation ##
#################################
#01.integrate data
#analyze for Figure 1G, integrate and cluster neurons
neu <- subset(seu,subset = celltype %in% c("c01:GABAergic Neurons","c02:Glutamatergic Neurons"))
mtx <- neu@assays$RNA@counts
meta <- neu@meta.data
neu <- CreateSeuratObject(mtx,min.cells = 10)
neu <- AddMetaData(neu,meta)
neu <- SplitObject(neu,split.by = "orig.ident")
neu <- lapply(neu,function(x){
  obj <- SCTransform(x, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
})
features <- SelectIntegrationFeatures(object.list = neu, nfeatures = 3000)  
neu <- PrepSCTIntegration(object.list = neu, anchor.features = features)
neu <- FindIntegrationAnchors(object.list = neu, normalization.method = "SCT",
                              reduction="rpca",k.anchor = 5,reference = c(1,5),
                              anchor.features = features)
neu <- IntegrateData(anchorset = neu, normalization.method = "SCT")
neu <- subset(neu, subset = percent.mt <=3)
neu <- RunPCA(neu, npcs = 30, verbose = FALSE)
neu <- RunUMAP(neu, dims = 1:15)
neu <- FindNeighbors(neu, reduction = "pca", dims = 1:15)
neu <- FindClusters(neu, resolution = 0.3)
DimPlot(neu,label = T)+theme(legend.position = "")
neu <- FindSubCluster(neu, cluster = "0", graph.name = "integrated_snn", resolution = 0.1)
Idents(neu) <- neu$sub.cluster
neu <- FindSubCluster(neu, cluster = "0_0", graph.name = "integrated_snn", resolution = 0.07)
Idents(neu) <- neu$sub.cluster
neu <- FindSubCluster(neu, cluster = "3", graph.name = "integrated_snn", resolution = 0.1)
Idents(neu) <- neu$sub.cluster
DimPlot(neu,group.by = "sub.cluster",label = T)+theme(legend.position = "")
neu$celltype <- NA
#re-assign celltype manually
idy <- Idents(neu) %in% c("2","8","5","18","17","0_2","0_0_0","13","3_0","10","12","9","16","15")
neu$celltype[idy] <- "c01:GABAergic Neurons"
neu$cluster[idy] <- "GABAergic Neurons"
neu$cluster_id[idy] <- "c01"
idy <- Idents(neu) %in% c("11","6","3_1","0_1","0_0_1","0_3","1","20","4","7","19","14")
neu$celltype[idy] <- "c02:Glutamatergic Neurons"
neu$cluster[idy] <- "Glutamatergic Neurons"
neu$cluster_id[idy] <- "c02"
DimPlot(neu,label = T,group.by = "celltype")+theme(legend.position = "")
neu@meta.data <- neu@meta.data[,1:13]
neu@meta.data$subcluster <- as.character(neu$seurat_clusters) 
neu <- PrepSCTFindMarkers(neu)
qs::qsave(neu,"03_neurons/s01_SEUobj_neurons annotated.qs") ##save step01 results
#02.2 find markers
#analyze for Figure 1G, identify neuronal cluster markers
DefaultAssay(neu) <- "SCT"
mk <- FindAllMarkers(neu, only.pos = T,logfc.threshold = 0.5)
qs::qsave(mk,"03_neurons/re02_markers.qs")

#02.3 cellcluster annotation --> subcluster
#analyze for Figure 1G, annotate neuronal subtypes using marker genes
mk <- qread("03_neurons/re02_markers.qs")
gene <- unique(mk$gene)
cluster <- unique(mk$cluster)
mtx <- matrix(0,nrow = length(gene),ncol = length(cluster)) %>% as.data.frame()
dimnames(mtx) <- c(list(gene),list(cluster))
for (i in cluster) {
  mki <- mk$gene[mk$cluster==i]
  mtx[mki,i] <- 1
}
mki <- vegan::vegdist(t(mtx),method = "jaccard") %>% as.matrix()  # cluster 1 and cluster 17 are similar
mk <- mk %>% arrange(desc(avg_log2FC)) %>% filter(!duplicated(gene)) %>% arrange(desc(pct = pct.1-pct.2))
markers <- c("Htr2c","Pmch","Rorb","Htr2c","Rprm","Slit3","Gpc6","Tac2","Satb2","Hdc","Alk","Unc13c",
             "Sox2ot","Lepr","Sst","Avp","Agrp","Cdh7","Sgcd","Gfra1","Brinp3","Epha6","Pomc","Tac1","Sst","Ghrh")
for(i in 1:26) neu$subcluster[Idents(neu)==sort(levels(Idents(neu)))[i]] <- paste0(markers[i],"+")
neu$subcluster <- paste(neu$cluster_id,neu$subcluster,sep = ":")
DimPlot(neu,label = T,group.by = "subcluster")+theme(legend.position = "")
qs::qsave(neu,"03_neurons/s01_SEUobj_neurons annotated.qs")
neu <- qread("03_neurons/s01_SEUobj_neurons annotated.qs")
DefaultAssay(neu) <- "RNA"
neu <- NormalizeData(neu)


#########################
## Analysis of neurons ##
#########################
#01.1 UMAP plot
#plot for Figure 1G, UMAP of neuronal subtypes
umap <- neu@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  mutate(cluster = neu$subcluster)
lvl <- c("c01:Agrp+","c01:Sst+","c01:Epha6+","c01:Ghrh+","c01:Satb2+","c01:Lepr+",
         "c01:Sgcd+","c01:Gpc6+","c01:Hdc+","c01:Htr2c+","c01:Sox2ot+","c01:Unc13c+",
         "c02:Slit3+","c02:Brinp3+","c02:Tac1+","c02:Cdh7+","c02:Rprm+","c02:Alk+","c02:Rorb+",
         "c02:Pmch+","c02:Avp+","c02:Tac2+","c02:Gfra1+","c02:Pomc+")
umap$cluster <- factor(umap$cluster,levels = lvl)
label_df <- umap %>% filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),.groups = "drop") %>%
  mutate(label = match(as.character(cluster), lvl))
ggplot(umap, aes(UMAP_1, UMAP_2)) +
  geom_point_rast(aes(color = cluster), size = 0.1) +
  scale_color_manual(values = c(paletteer_c("grDevices::Grays", 18)[c(2,4:11,12,14,16)],
                                paletteer_c("grDevices::Reds 3", 18)[c(2:10,12,14,16)])) +
  geom_point(data = label_df,shape = 21, fill = "white", 
             color = "white",size = 5.5, stroke = 0.4,show.legend = FALSE) +
  geom_text(data = label_df,aes(label = label),color = "black", size = 3,show.legend = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  xlab("") + ylab("") -> p; plot(p)
savePDF(p, "03_neurons/p01_neuron_umap.pdf",height = 6, width = 9)


#01.2 Marker genes plot
#analyze and plot for Figure S1J, dot plot of neuronal subtype markers
DefaultAssay(neu) <- "SCT"
Idents(neu) <- factor(neu$subcluster, levels = lvl)
mk <-  FindAllMarkers(neu, logfc.threshold = 0.25, only.pos = T)
qs::qsave(mk,"03_neurons/re02_markers_subclusters.qs")
mk <- mk %>% arrange(desc(avg_log2FC)) %>% filter(!duplicated(gene))
mk_plot <- lapply(lvl,function(x){
  marker <- mk[mk$cluster==x & mk$avg_log2FC > 0.3,] %>% filter(!grepl("^Gm|Rik$",gene))
  m1 <- gsub(".*:(.*)\\+","\\1",x)
  m2 <- marker[,"gene"] %>% .[!.%in%m1] %>% .[1]
  m3 <- marker[order(marker$pct.1-marker$pct.2,decreasing = T),"gene"] %>% .[!.%in%c(m1,m2)] %>% .[1]
  m <- unique(c(m1,m2,m3)) %>% na.omit()
  mpct <- mk[match(m,mk$gene),"avg_log2FC"]
  m <- m[order(mpct)]
})
p <- lapply(1:24,function(i){
  x <- mk_plot[[i]]
  if(length(x)>1){
    expr <- neu@assays$SCT@data[x,]%>%as.matrix()%>%t %>% as.data.frame()
  }else{
    expr <- neu@assays$SCT@data[x,]%>%as.matrix()%>%as.data.frame()
    colnames(expr) <- x
  }
  expr$group <- Idents(neu)
  f <- formula(paste0("cbind(",paste(x,collapse = ","),") ~ group"))
  re <- aggregate(f,expr,mean) %>% as_tibble()
  re[,-1] <- apply(re[,-1],2,scale)
  re <- reshape2::melt(re)
  if(length(x)>1){
    expr <- neu@assays$SCT@data[x,]%>%as.matrix()%>%t %>% as.data.frame()
  }else{
    expr <- neu@assays$SCT@data[x,]%>%as.matrix()%>%as.data.frame()
    colnames(expr) <- x
  }
  expr[expr>0] <- 1
  expr$group <- Idents(neu)
  re2 <- aggregate(f,expr,sum) %>% as_tibble()
  re2[,-1] <- re2[,-1]/as.numeric(table(Idents(neu)))
  re2 <- reshape2::melt(re2)
  mk <- qread("03_neurons/re02_markers_subclusters.qs")
  re3 <- lapply(x,function(xx){
    mki <- mk[mk$gene==xx,]
    re3 <- data.frame(group = levels(Idents(neu)), pv = 0)
    re3$pv[re3$group %in% mki$cluster] <- 1
    return(re3)
  })%>%purrr::reduce(rbind)
  re <- data.frame(re, pct = re2$value, pv = re3$pv)
})%>%purrr::reduce(rbind)
p$variable <- factor(p$variable,levels = unique(p$variable))
p$value[p$value<0] <- 0
p$value[p$value>3] <- 3
ggplot(p,aes(variable,group))+
  geom_point(aes(color=value,size=pct,alpha = pv))+
  scale_alpha_continuous(range = c(0,1))+
  scale_size_continuous(range = c(0,4))+
  scale_color_paletteer_c("grDevices::Reds 3",direction = -1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank())+
  xlab("")+ylab("") ->p; plot(p)
savePDF(p, "03_neurons/p02_marker_plot.pdf",height = 5, width = 12)



#02 miloR
#analyze and plot for Figure 1J, Milo differential abundance of neuronal neighborhoods
neu[["sample"]] <- neu$orig.ident
sel_kp(neu, k = 30, p = 0.1, d = 30)
milore <- run_miloR(neu, k = 30, p = 0.1, d = 30)
qsave(milore,"03_neurons/s02_milore.qs")
p <- plotnodegraph(milore, neu); plot(p)
savePDF(p,"03_neurons/p03.1_neurons_miloR.pdf",width = 5.6, height = 4)



#03. DEGs analysis in c01:GABAergic Neurons based on muscat
#analyze for Figure 1H, identify GABAergic neuronal DEGs with muscat
library(muscat)
neu1 <- subset(neu, subset = celltype == "c01:GABAergic Neurons")
neu1[["batch"]] <- ifelse(grepl("^AM",neu1$sample),"batch1","batch2")
#02.1 DEGs in c01:GABAergic Neurons
seuobj <- neu1
DefaultAssay(seuobj) <- "RNA"
seuobj <- seuobj[rowMeans(seuobj)>0.25,]
message(paste0("detect ", nrow(seuobj), " genes"))
sce <- as.SingleCellExperiment(seuobj)
colData(sce)$sample_id <- colData(sce)[["orig.ident"]]
colData(sce)$group_id  <- colData(sce)[["sample_group"]]
colData(sce)$batch_id  <- colData(sce)[["batch"]]
colData(sce)$cluster_id <- "1"
sce$sample_id <- factor(sce$sample_id)
sce$group_id  <- factor(sce$group_id)
sce$batch_id  <- factor(sce$batch_id)
sce$group_id <- relevel(sce$group_id, ref = "CFD")
sce <- prepSCE(sce,kid = "cluster_id",gid = "group_id",sid = "sample_id")
ps <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
metadata <- data.frame(sample_id = levels(sce$sample_id))
metadata$group_id <- sapply(metadata$sample_id,function(x) as.character(unique(sce$group_id[sce$sample_id == x])))
metadata$batch_id <- sapply(metadata$sample_id,function(x) as.character(unique(sce$batch_id[sce$sample_id == x])))
metadata$group_id <- factor(metadata$group_id)
metadata$batch_id <- factor(metadata$batch_id)
metadata$group_id <- relevel(metadata$group_id, ref = "CFD")
rownames(metadata) <- metadata$sample_id
design <- model.matrix(~ batch_id + group_id, data = metadata)
res_DS <- muscat::pbDS(ps, method = "edgeR", design = design, coef = 3, filter = "none")
pd <- muscat::aggregateData(sce, assay = "counts", fun = "num.detected", by = c("cluster_id", "sample_id"))
res_DD <- muscat::pbDS(pd, method = "DD", design = design, coef = 3,filter = "none")
res <- muscat::stagewise_DS_DD(res_DS,res_DD,sce = sce)
res <- data.frame(gene = res[[1]][[1]]$gene, p_adj = res[[1]][[1]]$p_adj,
                  DSp = res[[1]][[1]]$res_DS$p_val, DSpadj = res[[1]][[1]]$res_DS$p_adj.loc, DSfc = res[[1]][[1]]$res_DS$logFC,
                  DDp = res[[1]][[1]]$res_DD$p_val, DDpadj = res[[1]][[1]]$res_DD$p_adj.loc, DDfc = res[[1]][[1]]$res_DD$logFC,
                  row.names = res[[1]][[1]]$gene)
write.csv(res %>% arrange(p_adj), "03_neurons/muscat_DEG_detect_GABAergic Neurons.csv", row.names = F)
#02.2 DEGs in different subclusters 
#02.2.1 up DEGs 
#analyze and plot for Figure S1L, UMAP of upregulated gene scores in GABAergic neurons
degup <- res %>% filter((DSp < 0.05 & DSfc > 0.25) | (DDp < 0.05 & DDfc > 0.25)) %>% pull(gene)
neu1 <- AddModuleScore(neu1, list(degup), name = "response_up")
Featplot(neu1, "response_up1", group = "sample_group", size = 0.25) -> p; plot(p)
savePDF(p,"03_neurons/p07.1_response_to_treatment_up_UMAP_GABAergic Neurons.pdf",width = 5,height = 2.6)
#plot for Figure S1N, violin plots of upregulated gene scores in GABAergic subtypes
ggplot(neu1@meta.data, aes(subcluster, response_up1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.3, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  scale_fill_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"03_neurons/p07.1_response_to_treatment_up_violin_GABAergic Neurons.pdf",width = 7.5,height = 2.6)
#analyze and plot for Figure 1H, identify and map upward-shifted GABAergic neurons
model <- lmer(response_up1 ~ sample_group * subcluster + (1|batch), data = neu1@meta.data)
stat <- treat_lmer(model)
neu1[["shift_state"]] <- calc_shift(neu1, degup, group = "sample_group")
x1 <- (table(neu1$subcluster[neu1$sample_group=="HFD"], neu1$shift_state[neu1$sample_group=="HFD"]) %>% prop.table(margin = 1))*100
x2 <- (table(neu1$subcluster[neu1$sample_group=="CFD"], neu1$shift_state[neu1$sample_group=="CFD"]) %>% prop.table(margin = 1))*100
stat$effect_df$`shifted_cells(HFD)(%)` <- x1[,1]
stat$effect_df$`shifted_cells(CFD)(%)` <- x2[,1]
openxlsx::write.xlsx(stat, "03_neurons/p07.1_response_to_treatment_up_stat_GABAergic Neurons.xlsx",rowNames = F)
DimPlot(neu1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void() -> p; plot(p)
savePDF(p,"03_neurons/p07.1_response_to_treatment_up_shifted_state_GABAergic Neurons.pdf",width = 4.8, height = 2.6)
#02.2.2 dn DEGs 
#analyze and plot for Figure S1L, UMAP of downregulated gene scores in GABAergic neurons
degdn <- res %>% filter((DSp < 0.05 & DSfc < (-0.3)) | (DDp < 0.05 & DDfc < (-0.3))) %>% pull(gene)
neu1 <- AddModuleScore(neu1, list(degdn), name = "response_dn")
Featplot(neu1, "response_dn1", group = "sample_group", size =0.25) ->p; plot(p)
savePDF(p,"03_neurons/p07.2_response_to_treatment_dn_UMAP_GABAergic Neurons.pdf",width = 5,height = 2.6)
#plot for Figure S1M, violin plots of downregulated gene scores in GABAergic subtypes
ggplot(neu1@meta.data, aes(subcluster, response_dn1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.1, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  scale_fill_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"03_neurons/p07.2_response_to_treatment_dn_violin_GABAergic Neurons.pdf",width = 7.5,height = 2.6)
#analyze and plot for Figure 1H, identify and map downward-shifted GABAergic neurons
model <- lmer(response_dn1 ~ sample_group * subcluster + (1|batch), data = neu1@meta.data)
stat <- treat_lmer(model)
neu1[["shift_state"]] <- calc_shift(neu1, degdn, orient = -1, group = "sample_group",thres = 0.99)
x1 <- (table(neu1$subcluster[neu1$sample_group=="HFD"], neu1$shift_state[neu1$sample_group=="HFD"]) %>% prop.table(margin = 1))*100
x2 <- (table(neu1$subcluster[neu1$sample_group=="CFD"], neu1$shift_state[neu1$sample_group=="CFD"]) %>% prop.table(margin = 1))*100
stat$effect_df$`shifted_cells(HFD)(%)` <- x1[,1]
stat$effect_df$`shifted_cells(CFD)(%)` <- x2[,1]
openxlsx::write.xlsx(stat, "03_neurons/p07.2_response_to_treatment_dn_stat_GABAergic Neurons.xlsx",rowNames = F)
DimPlot(neu1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"03_neurons/p07.2_response_to_treatment_dn_shifted_state_GABAergic Neurons.pdf",width = 4.8,height = 2.6)



#04. DEGs analysis in c02:Glutamatergic Neurons based on muscat 
#analyze for Figure 1I, identify glutamatergic neuronal DEGs with muscat
neu2 <- subset(neu, subset = celltype == "c02:Glutamatergic Neurons")
neu2[["batch"]] <- ifelse(grepl("^AM",neu2$sample),"batch1","batch2")
#02.1 DEGs in c02:Glutamatergic Neurons
seuobj <- neu2
DefaultAssay(seuobj) <- "RNA"
seuobj <- seuobj[rowMeans(seuobj)>0.25,]
message(paste0("detect ", nrow(seuobj), " genes"))
sce <- as.SingleCellExperiment(seuobj)
colData(sce)$sample_id <- colData(sce)[["orig.ident"]]
colData(sce)$group_id  <- colData(sce)[["sample_group"]]
colData(sce)$batch_id  <- colData(sce)[["batch"]]
colData(sce)$cluster_id <- "1"
sce$sample_id <- factor(sce$sample_id)
sce$group_id  <- factor(sce$group_id)
sce$batch_id  <- factor(sce$batch_id)
sce$group_id <- relevel(sce$group_id, ref = "CFD")
sce <- prepSCE(sce,kid = "cluster_id",gid = "group_id",sid = "sample_id")
ps <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
metadata <- data.frame(sample_id = levels(sce$sample_id))
metadata$group_id <- sapply(metadata$sample_id,function(x) as.character(unique(sce$group_id[sce$sample_id == x])))
metadata$batch_id <- sapply(metadata$sample_id,function(x) as.character(unique(sce$batch_id[sce$sample_id == x])))
metadata$group_id <- factor(metadata$group_id)
metadata$batch_id <- factor(metadata$batch_id)
metadata$group_id <- relevel(metadata$group_id, ref = "CFD")
rownames(metadata) <- metadata$sample_id
design <- model.matrix(~ batch_id + group_id, data = metadata)
res_DS <- muscat::pbDS(ps, method = "edgeR", design = design, coef = 3, filter = "none")
pd <- muscat::aggregateData(sce, assay = "counts", fun = "num.detected", by = c("cluster_id", "sample_id"))
res_DD <- muscat::pbDS(pd, method = "DD", design = design, coef = 3,filter = "none")
res <- muscat::stagewise_DS_DD(res_DS,res_DD,sce = sce)
res <- data.frame(gene = res[[1]][[1]]$gene, p_adj = res[[1]][[1]]$p_adj,
                  DSp = res[[1]][[1]]$res_DS$p_val, DSpadj = res[[1]][[1]]$res_DS$p_adj.loc, DSfc = res[[1]][[1]]$res_DS$logFC,
                  DDp = res[[1]][[1]]$res_DD$p_val, DDpadj = res[[1]][[1]]$res_DD$p_adj.loc, DDfc = res[[1]][[1]]$res_DD$logFC,
                  row.names = res[[1]][[1]]$gene)
write.csv(res %>% arrange(p_adj), "03_neurons/muscat_DEG_detect_Glutamatergic Neurons.csv", row.names = F)
#02.2 DEGs in different subclusters 
#02.2.1 up DEGs 
#analyze and plot for Figure S1O, UMAP of upregulated gene scores in glutamatergic neurons
degup <- res %>% filter((DSp < 0.05 & DSfc > 0.25) | (DDp < 0.05 & DDfc > 0.25)) %>% pull(gene)
neu2 <- AddModuleScore(neu2, list(degup), name = "response_up")
Featplot(neu2, "response_up1", group = "sample_group", size = 0.25) -> p; plot(p)
savePDF(p,"03_neurons/p08.1_response_to_treatment_up_UMAP_Glutamatergic Neurons.pdf",width = 5,height = 2.6)
#plot for Figure S1Q, violin plots of upregulated gene scores in glutamatergic subtypes
ggplot(neu2@meta.data, aes(subcluster, response_up1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.1, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  scale_fill_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"03_neurons/p08.1_response_to_treatment_up_violin_Glutamatergic Neurons.pdf",width = 6.5,height = 2.6)
#analyze and plot for Figure 1I, identify and map upward-shifted glutamatergic neurons
model <- lmer(response_up1 ~ sample_group * subcluster + (1|batch), data = neu2@meta.data)
stat <- treat_lmer(model)
neu2[["shift_state"]] <- calc_shift(neu2, degup, group = "sample_group")
x1 <- table(neu2$subcluster[neu2$sample_group=="HFD"], neu2$shift_state[neu2$sample_group=="HFD"]) %>% prop.table(margin = 1)*100
x2 <- table(neu2$subcluster[neu2$sample_group=="CFD"], neu2$shift_state[neu2$sample_group=="CFD"]) %>% prop.table(margin = 1)*100
stat$effect_df$`shifted_cells(HFD)(%)` <- x1[,1]
stat$effect_df$`shifted_cells(CFD)(%)` <- x2[,1]
openxlsx::write.xlsx(stat, "03_neurons/p08.1_response_to_treatment_up_stat_Glutamatergic Neurons.xlsx",rowNames = F)
DimPlot(neu2, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void() -> p; plot(p)
savePDF(p,"03_neurons/p08.1_response_to_treatment_up_shifted_state_Glutamatergic Neurons.pdf",width = 4.8, height = 2.6)
#02.2.2 dn DEGs 
#analyze and plot for Figure S1O, UMAP of downregulated gene scores in glutamatergic neurons
degdn <- res %>% filter((DSp < 0.05 & DSfc < (-0.25)) | (DDp < 0.05 & DDfc < (-0.25))) %>% pull(gene)
neu2 <- AddModuleScore(neu2, list(degdn), name = "response_dn")
Featplot(neu2, "response_dn1", group = "sample_group", size =0.25) ->p; plot(p)
savePDF(p,"03_neurons/p08.2_response_to_treatment_dn_UMAP_Glutamatergic Neurons.pdf",width = 5,height = 2.6)
#plot for Figure S1P, violin plots of downregulated gene scores in glutamatergic subtypes
ggplot(neu2@meta.data, aes(subcluster, response_dn1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.1, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  scale_fill_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"03_neurons/p08.2_response_to_treatment_dn_violin_Glutamatergic Neurons.pdf",width = 6.5,height = 2.6)
#analyze and plot for Figure 1I, identify and map downward-shifted glutamatergic neurons
model <- lmer(response_dn1 ~ sample_group * subcluster + (1|batch), data = neu2@meta.data)
stat <- treat_lmer(model)
neu2[["shift_state"]] <- calc_shift(neu2, degdn, orient = -1, group = "sample_group")
x1 <- (table(neu2$subcluster[neu2$sample_group=="HFD"], neu2$shift_state[neu2$sample_group=="HFD"]) %>% prop.table(margin = 1))*100
x2 <- (table(neu2$subcluster[neu2$sample_group=="CFD"], neu2$shift_state[neu2$sample_group=="CFD"]) %>% prop.table(margin = 1))*100
stat$effect_df$`shifted_cells(HFD)(%)` <- x1[,1]
stat$effect_df$`shifted_cells(CFD)(%)` <- x2[,1]
openxlsx::write.xlsx(stat, "03_neurons/p08.2_response_to_treatment_dn_stat_Glutamatergic Neurons.xlsx",rowNames = F)
DimPlot(neu2, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"03_neurons/p08.2_response_to_treatment_dn_shifted_state_Glutamatergic Neurons.pdf",width = 4.8,height = 2.6)



#03.1 density plot # dont use that!
#analyze and plot, neuronal density on UMAP
umap <- data.frame(Embeddings(neu,reduction = "umap"),
                   sample_group = neu$sample_group,
                   cluster = neu$subcluster)
umap[umap$sample_group=="CFD","densi"] <- num_density(umap[umap$sample_group=="CFD",1:2],dist=0.2)
umap[umap$sample_group=="HFD","densi"] <- num_density(umap[umap$sample_group=="HFD",1:2],dist=0.2)
q <- min(quantile(umap[umap$sample_group=="CFD","densi"],0.99),
         quantile(umap[umap$sample_group=="HFD","densi"],0.99))
p <- lapply(c("CFD","HFD"),function(x){
  plotdata <- umap[umap$sample_group==x,]
  plotdata$densi[plotdata$densi>q] <- q
  p <- ggplot(plotdata,mapping = aes(UMAP_1,UMAP_2))+
    geom_point_rast(aes(color=densi),shape=".",raster.dpi = 600)+
    scale_color_paletteer_c("grDevices::Grays",direction = -1)+
    theme_void()+
    xlab("")+ylab("")
})
p <- wrap_plots(p,nrow = 2,guides = "collect"); plot(p)
savePDF(p,"03_neurons/p03_neurons_density.pdf",width = 4,height = 5)
#03.2 density chiqtest # dont use that
#analyze and plot for Figure S1R, chi-square residuals of neuronal subtype abundance
chitable <- table(umap$sample_group,umap$cluster)
chire <- chisq.test(chitable)
p <- data.frame(re = chire$stdres[1,], cluster = names(chire$stdres[1,])) %>% arrange(re)
p$cluster <- factor(p$cluster,levels = p$cluster)
col <- c(paletteer_c("grDevices::Grays", 18)[c(2,4:11,12,14,16)],paletteer_c("grDevices::Reds 3", 18)[c(2:10,12,14,16)])
names(col) <- lvl
ggplot(p,aes(cluster,-re,fill=cluster))+
  geom_col(position = "dodge")+
  scale_y_continuous(limits = c(-25,25))+
  scale_fill_manual(values = col)+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank())+
  labs(title = "X2test standardized residuals", x=NULL, y=NULL) -> p; plot(p)
savePDF(p,"03_neurons/p03_neurons_x2test.pdf",width = 5,height = 3)



#04 Fgfr2 expression #Dont use that!
#&neu@reductions$umap@cell.embeddings[,2]<(-10)

#plot for Figure 2V, Fgfr2 expression in Agrp-positive neurons
neu1 <- neu[,neu$subcluster=="c01:Agrp+" & neu@reductions$umap@cell.embeddings[,2]>(5) &
              neu@reductions$umap@cell.embeddings[,1]<(4)]
Geneplot(neu1[,neu1$sample_group=="CFD"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
savePDF(p,"03_neurons/p04_Agrp_Fgfr2_NC.pdf",width = 3,height = 2.4)
Geneplot(neu1[,neu1$sample_group=="HFD"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
savePDF(p,"03_neurons/p04_Agrp_Fgfr2_T.pdf",width = 3,height = 2.4)
#analyze and plot for Figure 2W, violin plot of Fgfr2 expression in Agrp-positive neurons
umap <- data.frame(Embeddings(neu1,reduction = "umap"),
                   cluster = neu1$subcluster,
                   sample_group = neu1$sample_group,
                   value = neu1@assays$SCT@data["Fgfr2",])
wilcox.test(umap$value[umap$sample_group=="CFD"], umap$value[umap$sample_group=="HFD"])
ggplot(umap,aes("c01:Agrp+",value,fill = sample_group))+
  geom_jitter(aes(color = sample_group),alpha = 0.6,size=0.5,
              position = position_jitterdodge(jitter.width = 0.3,dodge.width = 0.75))+
  geom_violin(width = 0.7,scale = "area",alpha = 0.7,adjust=1)+
  scale_y_continuous(limits = c(0,3.2))+
  scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
  scale_fill_manual(values = c("#fd8d3c","#4d8f87"))+
  theme_bw()+
  labs(title = "Agrp_Fgfr2: p < 2.2e-16 ", x = NULL)-> p; plot(p)
savePDF(p,"03_neurons/p04_Agrp_Fgfr2_violin.pdf",width = 3.5,height = 2.2)

#analyze and plot, Fgfr2 expression in Pomc-positive neurons
neu1 <- neu[,neu$subcluster=="c02:Pomc+"&
              neu@reductions$umap@cell.embeddings[,1]<(3.5) & neu@reductions$umap@cell.embeddings[,1]>(-1.5) &
              neu@reductions$umap@cell.embeddings[,2]<(6) & neu@reductions$umap@cell.embeddings[,2]>(1)]
Geneplot(neu1[,neu1$sample_group=="CFD"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
savePDF(p,"03_neurons/p04_Pomc_Fgfr2_NC.pdf",width = 3,height = 2.4)
Geneplot(neu1[,neu1$sample_group=="HFD"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
savePDF(p,"03_neurons/p04_Pomc_Fgfr2_T.pdf",width = 3,height = 2.4)
umap <- data.frame(Embeddings(neu1,reduction = "umap"),
                   cluster = neu1$subcluster,
                   sample_group = neu1$sample_group,
                   value = neu1@assays$SCT@data["Fgfr2",])
wilcox.test(umap$value[umap$sample_group=="CFD"], umap$value[umap$sample_group=="HFD"])
ggplot(umap,aes("c02:Pomc+",value,fill = sample_group))+
  geom_jitter(aes(color = sample_group),alpha = 0.6,size=0.5,
              position = position_jitterdodge(jitter.width = 0.3,dodge.width = 0.75))+
  geom_violin(width = 0.7,scale = "area",alpha = 0.7,adjust=1)+
  scale_y_continuous(limits = c(0,3.2))+
  scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
  scale_fill_manual(values = c("#fd8d3c","#4d8f87"))+
  theme_bw()+
  labs(title = "Pomc_Fgfr2: p = 3.985e-13 ", x = NULL)-> p; plot(p)
savePDF(p,"03_neurons/p04_Pomc_Fgfr2_violin.pdf",width = 3.5,height = 2.2)















