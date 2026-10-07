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
options(warn=0)
source("useful_functions.R")
source("/home/disk1/zzm/reference/rscripts tools/Fun_run seurat.R")


#########################
## Tanycytes extract ####
#########################
#01 integrate data with reference
#analyze for Figure 1B, extract and recluster tanycytes
seu <- qread("01_clustering/s01_Integration Data.qs")
#re-cluster
tany <- subset(seu,subset = celltype == "c08:Tanycytes")
tany <- reCreateSeuratObject(tany,dims = 10)
tany[["dataset"]] <- "nfyy_batch1"
tany$dataset[tany$sample%in%c("HX327","HX328","HX329","HX330")] <- "nfyy_batch2"
DimPlot(tany,group.by = "dataset")
#remove cluster3/4 because they express neuron markers: Tubb3, Snap25, Nxph1, Tenm1, Adarb2...
tany <- FindNeighbors(tany, reduction = "pca", dims = 1:10)
tany <- FindClusters(tany, resolution = 0.1)
DimPlot(tany,label = T)
tany <- subset(tany, subset = seurat_clusters %in% c(0,1,2))
DimPlot(tany,group.by = "dataset")


#02.combine with reference
#analyze for Figure 1B, integrate with reference and annotate tanycyte subtypes
ref <- qread("ref/Tanycytes.qs")
ref@meta.data[,6:11] <- NULL
ref[["sample_group"]] <- ref[["dataset"]] <- "reference"
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
tany <- RunUMAP(tany, dims = 1:20)
DimPlot(tany,group.by = "cluster",split.by = "sample_group",ncol = 2,label = T)
tany <- FindNeighbors(tany, reduction = "pca", dims = 1:20)
tany <- FindClusters(tany, resolution = 0.4)
DimPlot(tany,split.by = "dataset",label = T)
tany$subcluster <- tany$seurat_clusters %>% as.character()
cluster_name <- c("B2","A2","B1","B2","A1","Epen","Epen","A1")
for(i in 0:7) tany$subcluster[tany$subcluster==i] <- cluster_name[i+1]
Idents(tany) <- tany$subcluster
tany <- FindSubCluster(tany, c("A2"), graph.name = "integrated_snn", resolution = 0.2)
DimPlot(tany,label = T, group.by = "sub.cluster")
tany$subcluster[tany$sub.cluster%in%c("A2_1","A2_0")] <- "A2"
tany$subcluster[tany$sub.cluster%in%c("A2_2")] <- "A1"
Idents(tany) <- tany$subcluster
DimPlot(tany,label = T)
tany <- PrepSCTFindMarkers(tany)
qsave(tany,"02_subclustering/s01_Tanycytes_combined_w_ref.qs")

#1.2.1 trajectory
#analyze for Figure 1D, infer the tanycyte pseudospatial gradient
cds <- create_cds(tany)
cds <- Transfer_cds(tany,cds)
cds <- learn_graph_cds(cds, branch = 20)
cds <- order_cells(cds)
plot_cells(cds,color_cells_by = "pseudotime")
qsave(cds,"02_subclustering/s01_Tanycytes_combined_w_ref_cds.qs")


###########################
## Analysis of tanycytes ##
###########################
tany <- qread("02_subclustering/s01_Tanycytes_combined_w_ref.qs")
DefaultAssay(tany) <- "RNA"
tany0 <- subset(tany, subset = dataset=="reference")
tany1 <- subset(tany, subset = dataset!="reference")
cds <- qread("02_subclustering/s01_Tanycytes_combined_w_ref_cds.qs")

#01.markers
#01.1 subclusters cellMarker
#analyze, identify tanycyte subtype markers
Idents(tany1) <- factor(tany1$subcluster)
marker <- FindAllMarkers(tany1, only.pos = T)
openxlsx::write.xlsx(marker,"02_subclustering/stat01_Tanycytes_cell_subclusters_markers.xlsx")
#01.2 sample_group degs
#analyze, identify condition-associated DEGs in tanycytes
Idents(tany1) <- factor(tany1$sample_group)
deg <- FindMarkers(tany1,ident.1 = c("HFD"),logfc.threshold = 0) %>% mutate(gene = rownames(.))
openxlsx::write.xlsx(deg,"02_subclustering/stat02_Tanycytes_cell_sampleGroup_DEGs.xlsx")
#01.3 sample_group & subclusters degs
#analyze, identify condition-associated DEGs within tanycyte subtypes
Idents(tany1) <- factor(tany1$sample_group)
deg <- pblapply(c("B2","B1","A2","A1"), function(x){
  n <- subset(tany1,subset = subcluster%in%x)
  m <- FindMarkers(n,ident.1 = c("HFD"),logfc.threshold = 0) %>% mutate(gene = rownames(.)) %>%
    filter(p_val_adj < 0.05)
})
names(deg) <- c("B2","B1","A2","A1")
openxlsx::write.xlsx(deg,"02_subclustering/stat03_Tanycytes_cell_sampleGroup_subcluster_DEGs.xlsx")

#02.DEGs analysis based on muscat
#analyze and plot for Figure S2F, pseudobulk differential state and detection tests
require(muscat)
require(SingleCellExperiment)
seuobj <- tany1
seuobj <- subset(tany1,subset = subcluster != "Epen")
DefaultAssay(seuobj) <- "RNA"
seuobj <- seuobj[rowMeans(seuobj)>0.25,]
seuobj[["batch"]] <- seuobj$dataset
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
plot_res <- res %>%
  mutate(x = -log10(DSpadj), y = -log10(DDpadj), x1 = DSfc, y1 = DDfc, z = p_adj,
         col = case_when(p_adj < 0.1 & DSfc > 0  & DDfc > 0 ~ "upregulated in HFD",
                         p_adj < 0.1 & DSfc < 0  & DDfc < 0 ~ "downregulted in HFD",
                         T ~ "not sig"),
         z = case_when(p_adj<0.1 ~ "padj<0.1", T~"not significant")) %>% arrange(desc(p_adj))
ggplot(plot_res,aes(x1, y1, col = col))+
  geom_point(aes(size = z))+
  geom_hline(yintercept = c(-0.3,0.3), lty = "dashed")+
  geom_vline(xintercept = c(-0.4,0.4), lty = "dashed")+
  ggrepel::geom_text_repel(data = plot_res[plot_res$gene=="Fgf10",], aes(label = gene),
                           box.padding = 2)+
  scale_color_manual(values = c("blue","#D6D6CE","#CC0C00"))+
  scale_size_manual(values = c(0.9,1.5))+
  scale_x_continuous(expand = c(0,0),limits = c(-1.5,1.2))+
  scale_y_continuous(expand = c(0,0),limits = c(-1.1,1.2))+
  theme_classic()+
  labs(x = "Log2FC in Differentail state test",
       y = "Log2FC in differential detection analysis") -> p; plot(p)
savePDF(p,"02_subclustering/p05_muscat_DEG_detect.pdf",width = 5,height = 3)
write.csv(res %>% arrange(p_adj), "02_subclustering/p05_muscat_DEG_detect.csv", row.names = F)
#analyze and plot for Figure S2G, overlap of DEGs between batches
deg1 <- findDEG_muscat(tany1[,tany1$dataset=="nfyy_batch1"], sample = "orig.ident", group = "sample_group")
g1 <- deg1 %>% filter(DSp<0.05 & DDp<0.05 & abs(DSfc)>0.4 & abs(DDfc)>0.3) %>% pull(gene)
deg2 <- findDEG_muscat(tany1[,tany1$dataset=="nfyy_batch2"], sample = "orig.ident", group = "sample_group")
g2 <- deg2 %>% filter(DSp<0.05 & DDp<0.05 & abs(DSfc)>0.4 & abs(DDfc)>0.3) %>% pull(gene)
ggvenn(data = list(`DEGs in batch1` = g1, `DEGs in batch2` = g2)) ->p; plot(p)
savePDF(p, "02_subclustering/p05.1_venn.pdf", width = 5,height = 3)
df <- data.frame(gene = c(g1,g2), dataset = c(rep("DEGs in batch1", length(g1)),rep("DEGs in batch2", length(g2))))
write.csv(df, "02_subclustering/p05.1_venn_source.csv", row.names = F)

#02.2 DEGs in different subclusters 
#02.2.1 up DEGs 
#analyze and plot for Figure S1F, UMAP of upregulated gene module scores
degup <- res %>% filter(p_adj<0.15 & (DSfc > 0 | DDfc > 0)) %>% pull(gene)
tany1 <- AddModuleScore(tany1, list(degup), name = "response_up")
Featplot(tany1, "response_up1", group = "sample_group", size =0.8) ->p; plot(p)
savePDF(p,"02_subclustering/p02.1_response_to_treatment_up_UMAP.pdf",width = 5,height = 2.6)
#plot for Figure S1H, violin plots of upregulated gene module scores
ggplot(tany1@meta.data[tany1$subcluster!="Epen",], aes(subcluster, response_up1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.3, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  scale_fill_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"02_subclustering/p02.1_response_to_treatment_up_violin.pdf",width = 5,height = 2.6)
#analyze and plot for Figure 1E, identify and map upward-shifted tanycytes
model <- lme4::lmer(response_up1 ~ sample_group * subcluster + (1|dataset), data = tany1@meta.data[tany1$subcluster!="Epen",])
stat <- treat_lmer(model)
tany1[["shift_state"]] <- calc_shift(tany1, degup, group = "sample_group")
x1 <- (table(tany1$subcluster[tany1$sample_group=="HFD"], tany1$shift_state[tany1$sample_group=="HFD"]) %>% prop.table(margin = 1))*100
x2 <- (table(tany1$subcluster[tany1$sample_group=="CFD"], tany1$shift_state[tany1$sample_group=="CFD"]) %>% prop.table(margin = 1))*100
stat$effect_df$`shifted_cells(HFD)(%)` <- x1[1:4,1]
stat$effect_df$`shifted_cells(CFD)(%)` <- x2[1:4,1]
openxlsx::write.xlsx(stat, "02_subclustering/p02.1_response_to_treatment_up_stat.xlsx",rowNames = F)
DimPlot(tany1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"02_subclustering/p02.1_response_to_treatment_up_shiftedstate_UMAP.pdf",width = 4.8,height = 2.6)
#02.2.2 dn DEGs 
#analyze and plot for Figure S1F, UMAP of downregulated gene module scores
degdn <- res %>% filter(p_adj<0.15 & (DSfc < 0 | DDfc < 0)) %>% pull(gene)
tany1 <- AddModuleScore(tany1, list(degdn), name = "response_dn")
Featplot(tany1, "response_dn1", group = "sample_group", size =0.8) ->p; plot(p)
savePDF(p,"02_subclustering/p02.1_response_to_treatment_dn_UMAP.pdf",width = 5,height = 2.6)
#plot for Figure S1G, violin plots of downregulated gene module scores
ggplot(tany1@meta.data[tany1$subcluster!="Epen",], aes(subcluster, response_dn1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.3, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  scale_fill_manual(values = c(CFD = "#FFB900", HFD = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"02_subclustering/p02.1_response_to_treatment_dn_violin.pdf",width = 5,height = 2.6)
#analyze and plot for Figure 1E, identify and map downward-shifted tanycytes
model <- lme4::lmer(response_dn1 ~ sample_group * subcluster + (1|dataset), data = tany1@meta.data[tany1$subcluster!="Epen",])
stat <- treat_lmer(model)
tany1[["shift_state"]] <- calc_shift(tany1, degdn, orient = -1, group = "sample_group")
x1 <- (table(tany1$subcluster[tany1$sample_group=="HFD"], tany1$shift_state[tany1$sample_group=="HFD"]) %>% prop.table(margin = 1))*100
x2 <- (table(tany1$subcluster[tany1$sample_group=="CFD"], tany1$shift_state[tany1$sample_group=="CFD"]) %>% prop.table(margin = 1))*100
stat$effect_df$`shifted_cells(HFD)(%)` <- x1[1:4,1]
stat$effect_df$`shifted_cells(CFD)(%)` <- x2[1:4,1]
openxlsx::write.xlsx(stat, "02_subclustering/p02.1_response_to_treatment_dn_stat.xlsx",rowNames = F)
DimPlot(tany1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"02_subclustering/p02.1_response_to_treatment_dn_shiftedstate_UMAP.pdf",width = 4.8,height = 2.6)


#03. umap and marker gene
#plot for Figure 1B, UMAP of integrated tanycyte subtypes
umap <- data.frame(Embeddings(tany,reduction = "umap"),
                   dataset = tany$dataset,
                   cluster = tany$subcluster,
                   sample_group = tany$sample_group)
lvl <- c("Epen","A1","A2","B1","B2")
col <- c("#CC0C00","#D2854A","#EDD9A9","#91A17B","#486349")
umap$cluster <- factor(umap$cluster,levels = lvl)
ggplot(umap,aes(UMAP_1,UMAP_2))+
  # geom_point(data = umap[umap$dataset%in%"reference",],color="#E7E7E7",size=0.1,alpha=0.6)+
  # geom_point(data = umap[!umap$dataset%in%"reference",],aes(color=cluster),size=0.1)+
  geom_point(aes(color=cluster),size=0.1)+
  scale_color_manual(values = col)+
  theme_void()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("")+ylab("") -> p; plot(p)
savePDF(p,"02_subclustering/p01_tanycytes_allcells.pdf",width = 4,height = 3.2)
#markers: Slc1a2 Vcan Frzb Col25a1
#plot for Figure 1C, marker expression in reference and study tanycytes
Geneplot(tany0,"Col25a1",size = 0.1,stroke=0.2)[[1]] + scale_color_paletteer_c("grDevices::Oslo") -> p1; plot(p1)
Geneplot(tany1,"Col25a1",size = 0.1,stroke=0.2)[[1]] + scale_color_paletteer_c("grDevices::Oslo") -> p2;  plot(p2)
p <- p1+p2; plot(p)
savePDF(p,"02_subclustering/p01_tanycytes_Col25a1.pdf",width = 4,height = 2)


#04. pseudospatial
#plot for Figure 1D, tanycyte pseudospatial gradient on UMAP
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
savePDF(p1,"02_subclustering/p02_tanycytes_pseudospatial.pdf",width = 4,height = 3.2)
p2 <- ggplot(umap,aes(UMAP_1,UMAP_2,color=ps))+
  geom_point(data = umap[umap$sample_group=="reference",],color="#E7E7E7",size=0.01)+
  geom_point(data = umap[umap$sample_group=="HFD",],size=0.05)+
  scale_size_continuous(range = c(0.05,2))+
  scale_color_paletteer_c("viridis::magma")+
  theme_void()+
  theme(legend.key.width = unit(0.5,"cm"))+
  labs(x="",y="")
p2$layers <- c(p2$layers,p$layers[3]); plot(p2)
savePDF(p2,"02_subclustering/p02_tanycytes_pseudospatial_HFD.pdf",width = 4,height = 3.2)



#05 density
#05.1 density umap
#analyze and plot, tanycyte density on UMAP
umap[umap$sample_group=="CFD","densi"] <- num_density(umap[umap$sample_group=="CFD",1:2],dist=0.7)
umap[umap$sample_group=="HFD","densi"] <- num_density(umap[umap$sample_group=="HFD",1:2],dist=0.7)
q <- min(quantile(umap[umap$sample_group=="NC","densi"],0.99),
         quantile(umap[umap$sample_group=="T","densi"],0.99))
p <- lapply(c("CFD","HFD"),function(x){
  plotdata <- umap[umap$sample_group==x,]
  plotdata$densi[plotdata$densi>q] <- q
  p <- ggplot(plotdata,mapping = aes(UMAP_1,UMAP_2))+
    geom_point(aes(color=densi),size=0.1)+
    scale_color_paletteer_c("grDevices::Grays",direction = -1)+
    theme_void()+
    xlab("")+ylab("")
})
p <- wrap_plots(p,nrow = 2,guides = "collect"); plot(p)
savePDF(p,"02_subclustering/p02_tanycytes_density.pdf",width = 3,height = 5)
#05.2 density line
#plot, cell density along the pseudospatial gradient
ggplot(umap,aes(ps,densi,color=sample_group))+
  geom_point(data = umap[umap$sample_group=="CFD",],color="orange",alpha=0.5,size=0.2)+
  geom_point(data = umap[umap$sample_group=="HFD",],color="darkseagreen3",alpha=0.5,size=0.2)+
  geom_smooth(method = "gam",se=F,lwd=1.5)+
  scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
  theme_classic() -> p; plot(p)
savePDF(p,"02_subclustering/p02_tanycytes_density_line.pdf",width = 4,height = 2)
#05.3 density chiqtest
#analyze and plot for Figure S1I, chi-square residuals of tanycyte subtype abundance
chitable <- table(umap$sample_group,umap$cluster)[-3,]
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
savePDF(p,"02_subclustering/p02_tanycytes_x2test.pdf",width = 3,height = 3)
#05.4
#analyze and plot for Figure 1F, Milo differential abundance of tanycyte neighborhoods
tany1[["sample"]] <- tany1$orig.ident
sel_kp(tany1, k = 10, p = 0.1, d = 15)
milore <- run_miloR(tany1, k = 10, p = 0.1, d = 15)
p <- plotnodegraph(milore,tany1); plot(p)
savePDF(p,"02_subclustering/p02.2_tanycytes_miloR2.pdf",width = 5, height = 4)



#06 Fgf10 expression
#plot for Figure 2R, Fgf10 expression in tanycytes by condition
Geneplot(tany1[,tany1$sample_group=="CFD"],"Fgf10",size = 0.01)[[1]] -> p; plot(p)
savePDF(p,"02_subclustering/p04_tanycytes_Fgf10_CFD.pdf",width = 3,height = 2.4)
Geneplot(tany1[,tany1$sample_group=="HFD"][,sample(1:4225,2000)],"Fgf10",size = 0.01)[[1]] -> p; plot(p)
savePDF(p,"02_subclustering/p04_tanycytes_Fgf10_HFD.pdf",width = 3,height = 2.4)
#analyze and plot for Figure 2S, violin plot of tanycyte Fgf10 expression
umap <- data.frame(Embeddings(tany,reduction = "umap"),
                   cluster = tany$subcluster,
                   sample_group = tany$sample_group,
                   value = tany@assays$SCT@data["Fgf10",]) %>% 
  filter(sample_group!="reference" & cluster != "Epen")
wilcox.test(umap$value[umap$sample_group=="CFD"], umap$value[umap$sample_group=="HFD"])
ggplot(umap,aes("tanycytes",value,fill = sample_group))+
  geom_jitter(aes(color = sample_group),alpha = 0.6,size=0.5,
              position = position_jitterdodge(jitter.width = 0.3,dodge.width = 0.75))+
  geom_violin(width = 0.7,scale = "area",alpha = 0.7,adjust=2)+
  scale_y_continuous(limits = c(0,3.2))+
  scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
  scale_fill_manual(values = c("#fd8d3c","#4d8f87"))+
  theme_bw()+
  labs(title = "Tanycytes Fgf10: p = 7.284e-13", x = NULL)-> p; plot(p)-> p; plot(p)
savePDF(p,"02_subclustering/p04_tanycytes_Fgf10_violin.pdf",width = 3.5,height = 2.2)
#06.2 Fgf10 smoothed expression and positive rate
tany1.1 <- tany[,tany$dataset=="nfyy_batch1"]
#analyze and plot for Figure S2E, smoothed Fgf10 expression in batch 1
tany1.1[["Fgf10_dens"]] <- gene_density(tany1.1, "Fgf10", group = "orig.ident", dims = 15, neighbors = 10)
Featplot(tany1.1,"Fgf10_dens",group = "orig.ident", size = 1, raster = T) + scale_color_paletteer_c("grDevices::Blue-Red") -> p; plot(p)
savePDF(p,"02_subclustering/p05_Fgf10_smoothed expression_batch1.pdf",width = 9,height = 2.6)
#analyze and plot for Figure S2D, local Fgf10-positive cell rates in batch 1
tany1.1[["Fgf10_pos"]] <- positive_rate(tany1.1, "Fgf10", group = "orig.ident", dims = 15, neighbors = 10)
Featplot(tany1.1,"Fgf10_pos",group = "orig.ident", size = 1, raster = T) + scale_color_paletteer_c("grDevices::Blue-Red") -> p; plot(p)
savePDF(p,"02_subclustering/p05_Fgf10_positive rate_batch1.pdf",width = 9,height = 2.6)

tany1.2 <- tany[,tany$dataset=="nfyy_batch2"]
#analyze and plot for Figure S2E, smoothed Fgf10 expression in batch 2
tany1.2[["Fgf10_dens"]] <- gene_density(tany1.2, "Fgf10", group = "orig.ident", dims = 15, neighbors = 10)
Featplot(tany1.2,"Fgf10_dens",group = "orig.ident", size = 1, raster = T) + scale_color_paletteer_c("grDevices::Blue-Red") -> p; plot(p)
savePDF(p,"02_subclustering/p05_Fgf10_smoothed expression_batch2.pdf",width = 9,height = 2.6)
#analyze and plot for Figure S2D, local Fgf10-positive cell rates in batch 2
tany1.2[["Fgf10_dens"]] <- positive_rate(tany1.2, "Fgf10", group = "orig.ident", dims = 15, neighbors = 10)
Featplot(tany1.2,"Fgf10_dens",group = "orig.ident", size = 1, raster = T) + scale_color_paletteer_c("grDevices::Blue-Red") -> p; plot(p)
savePDF(p,"02_subclustering/p05_Fgf10_positive rate_batch2.pdf",width = 9,height = 2.6)


#07 literature geneset
#analyze and plot for Figure S1E, UMAP of literature-derived spatial gene set scores
gs <- readxl::read_xlsx("NEW/41467_2024_50913_MOESM6_ESM.xlsx") %>% filter(Condition == "Fed")
gs <- split(gs$Gene,gs$Classification)
gs <- sapply(gs, function(x) x[x%in%rownames(tany1)])
gs <- gs[c(9:11,1:3,5:8)]
tany1 <- AddModuleScore(tany1, gs)
colnames(tany1@meta.data)[26:ncol(tany1@meta.data)] <- names(gs)
Featplot(tany1,names(gs), size = 0.3, q = 0.95, nrow = 3)->p;plot(p)
savePDF(p, "02_subclustering/p01_location_score.pdf",width = 6.5,height = 4)



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
#analyze, simulate dropout and reconstruct clustering and pseudospatial gradients
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
  ref <- qread("ref/Tanycytes.qs")
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
#plot, detected gene counts before and after masking
library(ComplexHeatmap)
library(circlize)
tanyx <- seud_0.1
df_nf <- rbind(data.frame(nF = tany$nFeature_RNA, masked = tany$masked, group = "Before masking")[tany$dataset!="reference",],
               data.frame(nF = tanyx$nFeature_RNA, masked = tany$masked, group = "After masking")[tany$dataset!="reference",])
df_nf$group <- factor(df_nf$group, levels = c("Before masking", "After masking"))
ggplot(df_nf, aes(group, nF))+
  geom_point(aes(color = masked), size = 0.2, alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(aes(fill = masked), width = 0.7, alpha = 0.3, 
              color = "#616161",linewidth = 0.4, adjust = 2)+
  scale_color_manual(values = c(yes = "#CD534C", no = "#79AF97"))+
  scale_fill_manual(values = c(yes = "#F1F1F1", no = "#F1F1F1"))+
  theme_classic()+
  labs(x = NULL, y = "nFeature_RNA", title = "Masking rate = 0.5")-> p; plot(p)
savePDF(p, "02_subclustering/p06_dropout_0.5_nFeature_RNA.pdf",width = 4, height = 2.5)
#heatmap
#analyze and plot, cluster concordance after masking
tanyx <- FindClusters(tanyx, resolution = 0.2)
df_class <- table(tanyx$seurat_clusters[tany$masked=="yes"], tanyx$subcluster_raw[tany$masked=="yes"]) %>% 
  as.data.frame.matrix()
Heatmap(df_class, col = colorRamp2(seq(0,200,10),paletteer_c("grDevices::Blue-Yellow", 21)[21:1]),
        cluster_columns = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          text_col <- ifelse(df_class[i, j] > 150, "white", "black")
          grid.text(sprintf("%.0f", df_class[i, j]), x, y, gp = gpar(fontsize = 10, col = text_col))})
#p06_dropout_0.1_classification.pdf 4*3
#pseudospatial gradient
#analyze and plot, pseudospatial gradient concordance before and after masking
cdsx <- cds_0.1
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
       title = "Masking rate = 0.5") -> p; plot(p)
savePDF(p, "02_subclustering/p06_dropout_0.8_pseudospatial gradient.pdf",width = 3, height = 2.5)










