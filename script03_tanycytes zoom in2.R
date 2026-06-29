setwd("/home/disk1/wgs")
library(Seurat)
library(monocle3)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggsci)
library(sctransform)
library(glmGamPoi)
library(pbapply)
library(qs)
library(patchwork)
library(paletteer)
library(ComplexHeatmap)
library(circlize)
library(lme4)
options(future.globals.maxSize = 2048 * 1024^5)
source("useful_functions.R")
source("/home/disk1/zzm/reference/rscripts tools/Fun_run seurat.R")


#########################
## Tanycytes extract ####
#########################
#01 integrate data with reference
seu <- qread("NEW/s03_Integration Data mapped k5.qs")
#re-cluster
tany <- subset(seu,subset = celltype == "c06:Tanycytes")
tany <- reCreateSeuratObject(tany,dims = 10)
DimPlot(tany)
# We find a cluster which has high nFeature_RNA and splits from major cell cluster, so we remove it.
tany <- tany[,tany@reductions$umap@cell.embeddings[,2]<=5]
tany <- reCreateSeuratObject(tany,dims = 10)
DimPlot(tany)



########################
## Map into reference ##
########################
# Dont use this
#01.1 use FindTransferAnchors + MapQuery 
# ref <- qread("NEW/ref/Tanycytes.qs")
# anchors <- FindTransferAnchors(reference = ref, query = tany, dims = 1:20,
#                                reference.reduction = "pca",k.anchor = 20)
# tany <- MapQuery(anchorset = anchors, reference = ref, query = tany,
#                  refdata = list(celltype = "cluster"), reference.reduction = "pca", reduction.model = "umap")
# tany <- RunUMAP(tany,reduction = "ref.pca",dims = 1:20,min.dist = 0.2,n.neighbors = 20)
# DimPlot(tany,reduction = "ref.umap",group.by = "predicted.celltype")
# DimPlot(tany,reduction = "umap",group.by = "predicted.celltype")
# Geneplot(tany,"Slc1a2",q=0.9,size = 0.5)
# tany[["subcluster"]] <- tany$predicted.celltype
# Idents(tany) <- as.factor(tany$subcluster)
# qsave(tany,"NEW/Tanycytes/re01_Tanseu_MQ.qs")
# #1.1.1 trajectory
# cds <- create_cds(tany)
# cds <- Transfer_cds(tany,cds)
# cds <- learn_graph_cds(cds,branch = 15)
# cds <- order_cells(cds)
# plot_cells(cds,color_cells_by = "pseudotime")
# qsave(cds,"NEW/Tanycytes/re01_TanCDS_MQ.qs")


#01.2 use IntegrateData
ref <- qread("NEW/ref/Tanycytes.qs")
ref@meta.data[,6:11] <- NULL
ref[["sample_group"]] <- ref[["dataset"]] <- "reference"
tany[["dataset"]] <- "nfyy"
tany <- merge(tany,ref)
tany <- SplitObject(tany,split.by = "dataset")
tany <- lapply(tany,function(x){
  obj <- SCTransform(x, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = tany, nfeatures = 2000)  
tany <- PrepSCTIntegration(object.list = tany, anchor.features = features)
tany <- FindIntegrationAnchors(object.list = tany, normalization.method = "SCT",
                               reduction="rpca",k.anchor = 20, anchor.features = features)
tany <- IntegrateData(anchorset = tany, normalization.method = "SCT")
tany <- RunPCA(tany, npcs = 50, verbose = FALSE)
tany <- RunUMAP(tany, dims = 1:15)
DimPlot(tany,group.by = "cluster",split.by = "sample_group",ncol = 2,label = T)
Geneplot(tany,"Slc1a2",q=0.9,size = 0.5)
tany <- FindNeighbors(tany, reduction = "pca", dims = 1:15)
tany <- FindClusters(tany, resolution = 0.4)
tany$subcluster <- tany$seurat_clusters %>% as.character()
cluster_name <- c("A2","B2","B1","B2","A1","Epen","Epen","A1")
for(i in 0:7) tany$subcluster[tany$subcluster==i] <- cluster_name[i+1]
Idents(tany) <- tany$subcluster
DimPlot(tany,label = T)
tany <- PrepSCTFindMarkers(tany)
qsave(tany,"NEW/Tanycytes/re01_Tanseu_Ig.qs")
#1.2.1 trajectory
cds <- create_cds(tany)
cds <- Transfer_cds(tany,cds)
cds <- learn_graph_cds(cds, branch = 20)
cds <- order_cells(cds)
plot_cells(cds,color_cells_by = "pseudotime")
qsave(cds,"NEW/Tanycytes/re01_TanCDS_Ig.qs")



###########################
## Analysis of tanycytes ##
###########################
tany <- qread("NEW/Tanycytes/re01_Tanseu_Ig.qs")
DefaultAssay(tany) <- "RNA"
tany0 <- subset(tany, subset = dataset=="reference")
tany1 <- subset(tany, subset = dataset=="nfyy")
cds <- qread("NEW/Tanycytes/re01_TanCDS_Ig.qs")

#01.markers
#01.1 subclusters cellMarker
Idents(tany1) <- factor(tany1$subcluster)
marker <- FindAllMarkers(tany1, only.pos = T)
openxlsx::write.xlsx(marker,"NEW/Tanycytes/stat01_Tanycytes_cell_subclusters_markers.xlsx")
#01.2 sample_group degs
Idents(tany1) <- factor(tany1$sample_group)
deg <- FindMarkers(tany1,ident.1 = c("T"),logfc.threshold = 0) %>% mutate(gene = rownames(.))
qsave(deg,"NEW/Tanycytes/re03_DEGs_SampleLevel.qs")
openxlsx::write.xlsx(deg,"NEW/Tanycytes/stat01_Tanycytes_cell_sampleGroup_DEGs.xlsx")
#01.3 sample_group & subclusters degs
Idents(tany1) <- factor(tany1$sample_group)
deg <- pblapply(c("B2","B1","A2","A1"), function(x){
  n <- subset(tany1,subset = subcluster%in%x)
  m <- FindMarkers(n,ident.1 = c("T"),logfc.threshold = 0) %>% mutate(gene = rownames(.)) %>%
    filter(p_val_adj < 0.05)
})
names(deg) <- c("B2","B1","A2","A1")
openxlsx::write.xlsx(deg,"NEW/Tanycytes/stat01_Tanycytes_cell_sampleGroup_subcluster_DEGs.xlsx")

#02.DEGs analysis based on muscat
library(muscat)
#02.1 DEGs in all tanycytes
deg <- findDEG_muscat(tany1, sample = "orig.ident", group = "sample_group")
plot_res <- deg %>%
  mutate(x = -log10(DSpadj), y = -log10(DDpadj), x1 = abs(DSfc), y1 = abs(DDfc),
         col = case_when(p_adj < 0.1 & abs(DDfc) >0.3 & abs(DSfc) >0.5 & DSpadj <0.1 & DDpadj <0.1 ~ "1-both",
                         p_adj < 0.1 & abs(DSfc) >0.5 & DSpadj <0.1 ~ "2-sig. in DS",
                         p_adj < 0.1 & abs(DDfc) >0.3 & DDpadj <0.1 ~ "2-sig. in DD",
                         T ~ "not sig"))
ggplot(plot_res,aes(x, y, col = col))+
  geom_point(size = 0.75)+
  geom_hline(yintercept = 1, lty = "dashed")+
  geom_vline(xintercept = 1, lty = "dashed")+
  ggrepel::geom_text_repel(data = plot_res[plot_res$gene=="Fgf10",], aes(label = gene),
                           box.padding = 2)+
  scale_color_manual(values = c("#CC0C00","#FD8CC1","#FD7446","#D6D6CE"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x = "adjusted.pvalue in differential state analysis",
       y = "adjusted.pvalue in differential detection analysis") -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p05_muscat_DEG_detect.pdf",width = 4,height = 2.6)
write.csv(deg, "NEW/Tanycytes/p05_muscat_DEG_detect.csv", row.names = F)
#02.2 DEGs in different subclusters 
#02.2.1 up DEGs 
degup <- deg %>% filter(p_adj<0.1 & (DSfc > 0.5 | DDfc > 0.3)) %>% pull(gene)
tany1 <- AddModuleScore(tany1, list(degup), name = "response_up")
Featplot(tany1, "response_up1", group = "sample_group", size =0.8) ->p; plot(p)
savePDF(p,"NEW/Tanycytes/p02.1_response_to_treatment_up_UMAP.pdf",width = 5,height = 2.6)
ggplot(tany1@meta.data[tany1$subcluster!="Epen",], aes(subcluster, response_up1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.3, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  scale_fill_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"NEW/Tanycytes/p02.1_response_to_treatment_up_violin.pdf",width = 5,height = 2.6)
model <- lmer(response_up1 ~ sample_group * subcluster + (1|orig.ident), data = tany1@meta.data[tany1$subcluster!="Epen",])
stat <- treat_lmer(model)
openxlsx::write.xlsx(stat, "NEW/Tanycytes/p02.1_response_to_treatment_up_stat.xlsx",rowNames = F)
tany1[["shift_state"]] <- calc_shift(tany1, degup, group = "sample_group")
table(tany1$subcluster[tany1$sample_group=="T"], tany1$shift_state[tany1$sample_group=="T"]) %>% prop.table(margin = 1)
DimPlot(tany1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"NEW/Tanycytes/p02.1_response_to_treatment_up_shiftedstate_UMAP.pdf",width = 4.8,height = 2.6)
#02.2.2 dn DEGs 
degdn <- deg %>% filter(p_adj<0.1 & (DSfc < (-0.5) | DDfc < (-0.3))) %>% pull(gene)
tany1 <- AddModuleScore(tany1, list(degdn), name = "response_dn")
Featplot(tany1, "response_dn1", group = "sample_group", size =0.8) ->p; plot(p)
savePDF(p,"NEW/Tanycytes/p02.1_response_to_treatment_dn_UMAP.pdf",width = 5,height = 2.6)
ggplot(tany1@meta.data[tany1$subcluster!="Epen",], aes(subcluster, response_dn1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.3, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  scale_fill_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"NEW/Tanycytes/p02.1_response_to_treatment_dn_violin.pdf",width = 5,height = 2.6)
model <- lmer(response_dn1 ~ sample_group * subcluster + (1|orig.ident), data = tany1@meta.data[tany1$subcluster!="Epen",])
stat <- treat_lmer(model)
openxlsx::write.xlsx(stat, "NEW/Tanycytes/p02.1_response_to_treatment_dn_stat.xlsx",rowNames = F)
tany1[["shift_state"]] <- calc_shift(tany1, degdn, orient = -1, group = "sample_group")
table(tany1$subcluster[tany1$sample_group=="T"], tany1$shift_state[tany1$sample_group=="T"]) %>% prop.table(margin = 1)
DimPlot(tany1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"NEW/Tanycytes/p02.1_response_to_treatment_dn_shiftedstate_UMAP.pdf",width = 4.8,height = 2.6)



#03. umap and marker gene
umap <- data.frame(Embeddings(tany,reduction = "umap"),
                   dataset = tany$dataset,
                   cluster = tany$subcluster,
                   sample_group = tany$sample_group)
lvl <- c("Epen","A1","A2","B1","B2")
col <- c("#C7522B","#D2854A","#EDD9A9","#91A17B","#486349")
umap$cluster <- factor(umap$cluster,levels = lvl)
ggplot(umap,aes(UMAP_1,UMAP_2))+
  geom_point(data = umap[umap$dataset%in%"reference",],color="#E7E7E7",size=0.1,alpha=0.6)+
  geom_point(data = umap[umap$dataset%in%"nfyy",],aes(color=cluster),size=0.1)+
  # geom_point(aes(color=cluster),size=0.1)+
  scale_color_manual(values = col)+
  theme_void()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("")+ylab("") -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p01_tanycytes.pdf",width = 4,height = 3.2)
#markers: Slc1a2 Vcan Frzb Col25a1
Geneplot(tany0,"Fcgr2b",size = 0.1)[[1]] + scale_color_paletteer_c("grDevices::Oslo") -> p1; plot(p1)
Geneplot(tany1,"Pcp4",size = 0.1)[[1]] + scale_color_paletteer_c("grDevices::Oslo") -> p2;  plot(p2)
p <- p1+p2; plot(p)
savePDF(p,"NEW/Tanycytes/p01_tanycytes_Col25a1.pdf",width = 4,height = 2)



#04. pseudospatial
umap <- data.frame(Embeddings(tany,reduction = "umap"),
                   dataset = tany$dataset,
                   cluster = tany$subcluster,
                   sample_group = tany$sample_group,
                   ps = pseudotime(cds))
p <- plot_cells(cds,color_cells_by = "cluster",cell_size = 0,
                label_branch_points = F,label_leaves = F,label_roots = F); plot(p)
p1 <- ggplot(umap)+
  geom_point(aes(UMAP_1,UMAP_2,color=ps),size=0.05,alpha=0.6)+
  scale_size_continuous(range = c(0.05,2))+
  scale_color_paletteer_c("viridis::magma")+
  theme_void()+
  theme(legend.key.width = unit(0.5,"cm"))+
  labs(x="",y="")
p1$layers <- c(p1$layers,p$layers[3]); plot(p1)
savePDF(p1,"NEW/Tanycytes/p02_tanycytes_pseudospatial.pdf",width = 4,height = 3.2)
p2 <- ggplot(umap,aes(UMAP_1,UMAP_2,color=ps))+
  geom_point(data = umap[umap$sample_group=="reference",],color="#E7E7E7",size=0.01)+
  geom_point(data = umap[umap$sample_group=="NC",],size=0.05)+
  scale_size_continuous(range = c(0.05,2))+
  scale_color_paletteer_c("viridis::magma")+
  theme_void()+
  theme(legend.key.width = unit(0.5,"cm"))+
  labs(x="",y="")
p2$layers <- c(p2$layers,p$layers[3]); plot(p2)
savePDF(p2,"NEW/Tanycytes/p02_tanycytes_pseudospatial_NC.pdf",width = 4,height = 3.2)
#04.1 pseudospatial heatmap
cds <- cds[,cds@colData$dataset %in% "nfyy"]
gene_fit <- graph_test(cds, neighbor_graph="principal_graph", cores=30) %>% na.omit()
qsave(gene_fit,"NEW/Tanycytes/re02_tany_genefit.qs")
genes <- rownames(gene_fit)[gene_fit$morans_I>0.1]
ps <- sort(pseudotime(cds))
mtx <- cds@assays@data$counts[genes,names(ps)] %>% as.matrix()
mtx <- apply(mtx,1,function(x){smooth.spline(x,df=3)$y}) %>% scale %>% t
col <- colorRamp2(c(seq(2,0.5,-0.1),seq(0.5,-0.5,-0.075),seq(-0.5,-2,-0.5)), 
                  c(paletteer_c("grDevices::Blues 2", 30)[1:16],
                    paletteer_c("grDevices::Blues 2", 30)[17:30],
                    paletteer_c("grDevices::Grays", 30)[c(25,29,30,30)]))
genes_order <- sapply(genes,function(x){
  ps <- floor2(norm_scale(ps))*10 
  data <- data.frame(ps = as.character(ps))
  dum <- model.matrix(~ ps,data)[,-1] %>% as.data.frame()
  dum$ps0 <- ifelse(rowSums(dum)==0,1,0)
  dum <- dum[,paste0("ps",0:9)]
  colnames(dum)[which.max(cor(as.numeric(mtx[x,]),dum))]
})
genes_order <- factor(genes_order,levels = paste0("ps",0:9))
meta <- cds@colData$subcluster[order(pseudotime(cds))]
colanno <- HeatmapAnnotation(A1 = ifelse(meta=="A1","yes","no"),
                             A2 = ifelse(meta=="A2","yes","no"), 
                             B1 = ifelse(meta=="B1","yes","no"), 
                             B2 = ifelse(meta=="B2","yes","no"), 
                             col = list(
                               A1 = c("yes" = "#D2854A", "no" = "white"),
                               A2 = c("yes" = "#EDD9A9", "no" = "white"),
                               B1 = c("yes" = "#91A17B", "no" = "white"),
                               B2 = c("yes" = "#486349", "no" = "white")),border = T)
Heatmap(mtx, col=col, border = F, use_raster=T, 
        row_split = genes_order, row_gap = unit(0,"mm"), 
        cluster_rows = T, clustering_method_rows = "ward.D2",
        cluster_columns = F, top_annotation = colanno,
        show_row_names = F,show_column_names = F,show_row_dend = F,show_column_dend = F, width = unit(6, "cm"))




#05 density
#05.1 density umap
umap[umap$sample_group=="NC","densi"] <- num_density(umap[umap$sample_group=="NC",1:2],dist=0.7)
umap[umap$sample_group=="T","densi"] <- num_density(umap[umap$sample_group=="T",1:2],dist=0.7)
q <- min(quantile(umap[umap$sample_group=="NC","densi"],0.99),
         quantile(umap[umap$sample_group=="T","densi"],0.99))
p <- lapply(c("NC","T"),function(x){
  plotdata <- umap[umap$sample_group==x,]
  plotdata$densi[plotdata$densi>q] <- q
  p <- ggplot(plotdata,mapping = aes(UMAP_1,UMAP_2))+
    geom_point(aes(color=densi),size=0.1)+
    scale_color_paletteer_c("grDevices::Grays",direction = -1)+
    theme_void()+
    xlab("")+ylab("")
})
p <- wrap_plots(p,nrow = 2,guides = "collect"); plot(p)
savePDF(p,"NEW/Tanycytes/p02_tanycytes_density.pdf",width = 3,height = 5)
#05.2 density line
ggplot(umap,aes(ps,densi,color=sample_group))+
  geom_point(data = umap[umap$sample_group=="NC",],color="orange",alpha=0.5,size=0.2)+
  geom_point(data = umap[umap$sample_group=="T",],color="darkseagreen3",alpha=0.5,size=0.2)+
  geom_smooth(method = "gam",se=F,lwd=1.5)+
  scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
  theme_classic() -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p02_tanycytes_density_line.pdf",width = 4,height = 2)
#05.3 density chiqtest
chitable <- table(umap$sample_group,umap$cluster)[-2,]
chire <- chisq.test(chitable)
p <- data.frame(re = chire$stdres[1,], cluster = names(chire$stdres[1,]))
p$cluster <- factor(p$cluster,levels = lvl)
ggplot(p,aes(cluster,-re,fill=cluster))+
  geom_col(position = "dodge",width=1.01)+
  scale_y_continuous(limits = c(-20,20))+
  scale_fill_manual(values = col)+
  theme_minimal()+
  theme(panel.grid.major.x = element_blank())+
  labs(title = "X2test standardized residuals", x=NULL, y=NULL) -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p02_tanycytes_x2test.pdf",width = 3,height = 3)
#05.4
tany1[["sample"]] <- tany1$orig.ident
sel_kp(tany1, k = 10, p = 0.1, d = 15)
milore <- run_miloR(tany1, k = 10, p = 0.1, d = 15)
p <- plotnodegraph(milore); plot(p)
savePDF(p,"NEW/Tanycytes/p02.2_tanycytes_miloR2.pdf",width = 5, height = 4)



#06 Fgf10 expression
Geneplot(tany1[,tany1$sample_group=="NC"],"Fgf10",size = 0.01)[[1]] -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p04_tanycytes_Fgf10_NC.pdf",width = 3,height = 2.4)
Geneplot(tany1[,tany1$sample_group=="T"][,sample(1:4225,2000)],"Fgf10",size = 0.01)[[1]] -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p04_tanycytes_Fgf10_T.pdf",width = 3,height = 2.4)
umap <- data.frame(Embeddings(tany,reduction = "umap"),
                   cluster = tany$subcluster,
                   sample_group = tany$sample_group,
                   value = tany@assays$SCT@data["Fgf10",]) %>% 
  filter(sample_group!="reference" & cluster != "Epen")
ggplot(umap,aes("tanycytes",value,fill = sample_group))+
  geom_jitter(aes(color = sample_group),alpha = 0.6,size=0.5,
              position = position_jitterdodge(jitter.width = 0.3,dodge.width = 0.75))+
  geom_violin(width = 0.7,scale = "area",alpha = 0.7,adjust=2)+
  scale_y_continuous(limits = c(0,3.2))+
  scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
  scale_fill_manual(values = c("#fd8d3c","#4d8f87"))+
  theme_bw() -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p04_tanycytes_Fgf10_violin.pdf",width = 3.5,height = 2.2)
#06.2 Fgf10 smoothed expression and positive rate
tany1[["Fgf10_dens"]] <- gene_density(tany1, "Fgf10", group = "orig.ident", dims = 15, neighbors = 10)
Featplot(tany1,"Fgf10_dens",group = "orig.ident", size = 1, raster = T) + scale_color_paletteer_c("grDevices::Blue-Red") -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p05_Fgf10_smoothed expression.pdf",width = 9,height = 2.6)
tany1[["Fgf10_pos"]] <- positive_rate(tany1, "Fgf10", group = "orig.ident", dims = 15, neighbors = 10)
Featplot(tany1,"Fgf10_pos",group = "orig.ident", size = 1, raster = T) + scale_color_paletteer_c("grDevices::Blue-Red") -> p; plot(p)
savePDF(p,"NEW/Tanycytes/p05_Fgf10_positive rate.pdf",width = 9,height = 2.6)
# 57.1% cells' positive rate > 0.5 in NC while 45.0% in T




#07 literature geneset
tany1 <- subset(tany, subset = dataset=="nfyy")
gs <- readxl::read_xlsx("NEW/41467_2024_50913_MOESM6_ESM.xlsx") %>% filter(Condition == "Fed")
gs <- split(gs$Gene,gs$Classification)
gs <- sapply(gs, function(x) x[x%in%rownames(tany1)])
gs <- gs[c(9:11,1:3,5:8)]
tany1 <- AddModuleScore(tany1, gs)
colnames(tany1@meta.data)[17:ncol(tany1@meta.data)] <- names(gs)
Featplot(tany1,names(gs), size = 0.3, q = 0.95, nrow = 3)->p;plot(p)
savePDF(p, "NEW/Tanycytes/p01_location_score.pdf",width = 6.5,height = 4)



#08. dropout
#question: How would the authors know that the speudospatial gradient (Fig. 1C) is not a dropout artifact without performing anatomical stratification (in situ hybridization on candidate genes in control vs. HFD animals) using the most significant anchoring genes and quantify cell numbers.
#作为替代，我们进行了额外的in silico dropout模拟分析。我们首先将Tanycytes细胞根据测得基因数的中位数分为两组。
#对于基因数较低的组，我们假定其代表了真实测序过程中更高的dropout rate群体，即high droupout rate group
#反之，基因数较高的组为low droupout rate group. 
#我们通过对low droupout rate group中细胞添加随机掩码来模拟不同dropout rate（0.1, 0.3, 0.5, 0.8），并评估其对细胞分群和pseudospatial gradient分数的影响。
#我们发现当masked rate达到0.3时，即每个细胞随机丢失了30%的检测到的基因，low droupout rate group与High droupout rate group间中位基因数相近。
#在此模拟水平下，经掩码处理的细胞在分群结果和pseudospatial gradient分数上与原始数据具有高度一致性。
#事实上，在0.1至0.5的模拟范围内，细胞分群及pseudospatial gradient均表现出高度的稳健性，提示我们的单细胞数据中的基因覆盖度足以构建稳定的pseudospatial gradient。
#尽管计算模拟无法替代原位杂交的解剖学验证方法，但上述结果表明，我们观察到的梯度不太可能由随机技术性dropout所驱动。
#01. masked expression
mask_expr <- function(seuobj, rate = 0.1){
  idy <- seuobj$nFeature_RNA > median(seuobj$nFeature_RNA)
  d_rate <- rate
  mtxb <- mtx <- seuobj@assays$RNA@counts
  mtxb[mtxb>0] <- 1
  set.seed(123)
  mtxb@x[sample(1:length(mtx@x),d_rate*length(mtx@x))] <- 0
  mtx <- as.matrix(mtx * mtxb)
  mtx[,!idy] <- as.matrix(seuobj@assays$RNA@counts[,!idy])
  return(mtx)
}
map_with_ref <- function(mtx){
  seud <- CreateSeuratObject(mtx)
  ref <- qread("NEW/ref/Tanycytes.qs")
  ref@meta.data[,6:11] <- NULL
  ref[["sample_group"]] <- ref[["dataset"]] <- "reference"
  seud[["dataset"]] <- "nfyy"
  seud <- merge(seud,ref)
  seud <- SplitObject(seud,split.by = "dataset")
  seud <- lapply(seud,function(x){
    obj <- SCTransform(x, vst.flavor = "v2", verbose = FALSE) %>%
      RunPCA(npcs = 50, verbose = FALSE)
  })
  features <- SelectIntegrationFeatures(object.list = seud, nfeatures = 2000)  
  seud <- PrepSCTIntegration(object.list = seud, anchor.features = features)
  seud <- FindIntegrationAnchors(object.list = seud, normalization.method = "SCT",
                                 reduction="rpca",k.anchor = 20, anchor.features = features)
  seud <- IntegrateData(anchorset = seud, normalization.method = "SCT")
  seud <- RunPCA(seud, npcs = 50, verbose = FALSE)
  seud <- RunUMAP(seud, dims = 1:15, seed = 123)
  seud <- FindNeighbors(seud, reduction = "pca", dims = 1:15)
  seud <- FindClusters(seud, resolution = 0.4)
  seud[["subcluster_raw"]] <- tany$subcluster[colnames(seud)]
  return(seud)
}
make_pseudospatial_gradient <- function(seuobj){
  seuobj[["subcluster"]] <- seuobj[["subcluster_raw"]]
  cdsx <- create_cds(seuobj)
  cdsx <- Transfer_cds(seuobj,cdsx)
  cdsx <- learn_graph_cds(cdsx, branch = 20)
}
cells <- colnames(tany1)[tany1$nFeature_RNA > median(tany1$nFeature_RNA)]
tany[["masked"]] <- ifelse(colnames(tany)%in%cells,"yes","no")
tany[["ps"]] <- pseudotime(cds)
seud_0.1 <- mask_expr(tany1, rate = 0.1) %>% map_with_ref
seud_0.3 <- mask_expr(tany1, rate = 0.3) %>% map_with_ref
seud_0.5 <- mask_expr(tany1, rate = 0.5) %>% map_with_ref
seud_0.8 <- mask_expr(tany1, rate = 0.8) %>% map_with_ref
cds_0.1 <- make_pseudospatial_gradient(seud_0.1) %>% order_cells
cds_0.3 <- make_pseudospatial_gradient(seud_0.3) %>% order_cells
cds_0.5 <- make_pseudospatial_gradient(seud_0.5) %>% order_cells
cds_0.8 <- make_pseudospatial_gradient(seud_0.8) %>% order_cells

#violin
library(ComplexHeatmap)
library(circlize)
tanyx <- seud_0.8
df_nf <- rbind(data.frame(nF = tany$nFeature_RNA, masked = tany$masked, group = "Before masking")[tany$dataset=="nfyy",],
               data.frame(nF = tanyx$nFeature_RNA, masked = tany$masked, group = "After masking")[tany$dataset=="nfyy",])
df_nf$group <- factor(df_nf$group, levels = c("Before masking", "After masking"))
ggplot(df_nf, aes(group, nF))+
  geom_point(aes(color = masked), size = 0.2, alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(aes(fill = masked), width = 0.7, alpha = 0.3, 
              color = "#616161",linewidth = 0.4, adjust = 2)+
  scale_color_manual(values = c(yes = "#CD534C", no = "#79AF97"))+
  scale_fill_manual(values = c(yes = "#F1F1F1", no = "#F1F1F1"))+
  theme_classic()+
  labs(x = NULL, y = "nFeature_RNA", title = "Masking rate = 0.8")-> p; plot(p)
savePDF(p, "NEW/Tanycytes/p06_dropout_0.8_nFeature_RNA.pdf",width = 4, height = 2.5)
#heatmap
df_class <- table(tanyx$seurat_clusters[tany$masked=="yes"], tanyx$subcluster_raw[tany$masked=="yes"]) %>% 
  as.data.frame.matrix() %>% .[,c("A1","A2","B1","B2")]
Heatmap(df_class, col = colorRamp2(seq(0,200,10),paletteer_c("grDevices::Blue-Yellow", 21)[21:1]),
        cluster_columns = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          text_col <- ifelse(df_class[i, j] > 150, "white", "black")
          grid.text(sprintf("%.0f", df_class[i, j]), x, y, gp = gpar(fontsize = 10, col = text_col))})
#p06_dropout_0.1_classification.pdf 4*3
#pseudospatial gradient
cdsx <- cds_0.8
df_ps <- data.frame(ps1 = pseudotime(cds), ps2 = pseudotime(cdsx),
                    masked = tany$masked) %>% filter(masked == "yes")
lm_fit <- lm(ps2 ~ ps1, data = df_ps)
r2 <- summary(lm_fit)$r.squared
RMSE <- sqrt(mean(residuals(lm_fit)^2))
ggplot(df_ps, aes(ps1, ps2))+
  geom_point(color = "#919191")+
  geom_smooth(method = "lm", se = F, color = "black", lty = "dashed")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  annotate("text", x = 5, y = Inf, 
           label = paste0("R² = ", round(r2,3), "\nRMSE = ", round(RMSE,3)), 
           hjust = 0, vjust = 1.2, size = 4)+
  theme_classic()+
  labs(x = "pseudospatial gradient score \n (Before masking)", 
       y = "pseudospatial gradient score \n (After masking)", 
       title = "Masking rate = 0.8") -> p; plot(p)
savePDF(p, "NEW/Tanycytes/p06_dropout_0.8_pseudospatial gradient.pdf",width = 3, height = 2.5)

