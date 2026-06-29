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

#01.Literature markers
seu <- qread("NEW/s03_Integration Data mapped k5.qs")
feat <- readxl::read_xlsx("step03_cluster zoom in/Neurons/Literature Markers.xlsx") %>% as.data.frame()
feat <- lapply(1:nrow(feat), function(x){
  gene <- feat$genes[x] %>% stringr::str_split(";") %>% unlist %>% .[.%in%rownames(seu@assays$RNA)]
})
names(feat) <- name


#################################
## integrate data & annotation ##
#################################
#01.integrate data
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
                              reduction="rpca",k.anchor = 5,reference = c(1,3),
                              anchor.features = features)
neu <- IntegrateData(anchorset = neu, normalization.method = "SCT")
neu <- RunPCA(neu, npcs = 30, verbose = FALSE)
neu <- RunUMAP(neu, dims = 1:20,spread = 2)
DimPlot(neu,group.by = "sample_group")
neu <- FindNeighbors(neu, reduction = "pca", dims = 1:20)
neu <- FindClusters(neu, resolution = 0.4)
DimPlot(neu,label = T)+theme(legend.position = "")
neu$celltype <- NA
#re-assign celltype manually
idy <- neu$seurat_clusters %in% c("11","16","18","7","2","15","6","0","12","19","10","8")
neu$celltype[idy] <- "c01:GABAergic Neurons"
neu$cluster[idy] <- "GABAergic Neurons"
neu$cluster_id[idy] <- "c01"
idy <- neu$seurat_clusters %in% c("1","17","3","5","13","4","9","14")
neu$celltype[idy] <- "c02:Glutamatergic Neurons"
neu$cluster[idy] <- "Glutamatergic Neurons"
neu$cluster_id[idy] <- "c02"
DimPlot(neu,label = T,group.by = "celltype")+theme(legend.position = "")
neu@meta.data <- neu@meta.data[,1:13]
neu@meta.data$subcluster <- as.character(neu$seurat_clusters) 
neu <- PrepSCTFindMarkers(neu)
qs::qsave(neu,"NEW/Neurons/re01_SEUobj_neurons annotated.qs") ##save step01 results
#02.2 find markers
DefaultAssay(neu) <- "SCT"
mk <- FindAllMarkers(neu, only.pos = T)
qs::qsave(mk,"NEW/Neurons/re02_markers.qs")

#02.3 cellcluster annotation --> subcluster
mk <- qread("NEW/Neurons/re02_markers.qs")
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
markers <- c("Agrp","Gpr149","Epha6","Tcf4","Pomc","Tmem163",
             "Sst","Tac1","Sv2c","Tac2","Ndst4",
             "Ghrh","Satb2","Gfra1","Avp","Kcnq5",
             "Utrn","Gpr149","Lef1","Gpc5")
for(i in 0:19) neu$subcluster[neu$subcluster==i] <- paste0(markers[i+1],"+")
neu$subcluster <- paste(neu$cluster_id,neu$subcluster,sep = ":")
DimPlot(neu,label = T,group.by = "subcluster")+theme(legend.position = "")
qs::qsave(neu,"NEW/Neurons/re01_SEUobj_neurons annotated.qs")
neu <- qread("NEW/Neurons/re01_SEUobj_neurons annotated.qs")
DefaultAssay(neu) <- "RNA"
neu <- NormalizeData(neu)


#########################
## Analysis of neurons ##
#########################
#01.1 UMAP plot
umap <- neu@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  mutate(cluster = neu$subcluster)
lvl <- c("c01:Agrp+","c01:Sst+","c01:Kcnq5+","c01:Epha6+","c01:Gpc5+",
         "c01:Tac1+","c01:Ndst4+","c01:Satb2+","c01:Sv2c+","c01:Lef1+","c01:Utrn+",
         "c01:Ghrh+","c02:Gpr149+","c02:Tac2+","c02:Avp+","c02:Tcf4+","c02:Tmem163+",
         "c02:Gfra1+","c02:Pomc+")
umap$cluster <- factor(umap$cluster,levels = lvl)
ggplot(umap,aes(UMAP_1,UMAP_2))+
  geom_point_rast(aes(color=cluster),size=0.1)+
  scale_color_manual(values = c(paletteer_c("grDevices::Grays", 18)[c(2,4:14)],
                                paletteer_c("grDevices::Reds 3", 12)[3:9]))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("")+ylab("")

#01.2 Marker genes plot
DefaultAssay(neu) <- "SCT"
Idents(neu) <- factor(neu$subcluster, levels = lvl)
mk <-  FindAllMarkers(neu, logfc.threshold = 0.25)
qs::qsave(mk,"NEW/Neurons/re02_markers_subclusters.qs")
mk <- mk %>% arrange(desc(avg_log2FC)) %>% filter(!duplicated(gene))
mk_plot <- lapply(lvl,function(x){
  marker <- mk[mk$cluster==x & mk$avg_log2FC > 0.5,] %>% filter(!grepl("^Gm|Rik$",gene))
  m1 <- gsub(".*:(.*)\\+","\\1",x)
  m2 <- marker[,"gene"] %>% .[!.%in%m1] %>% .[1:2]
  m3 <- marker[order(marker$pct.1-marker$pct.2,decreasing = T),"gene"] %>% .[!.%in%c(m1,m2)] %>% .[1]
  m <- unique(c(m1,m2,m3))
  mpct <- mk[match(m,mk$gene),"avg_log2FC"]
  m <- m[order(mpct)]
})
p <- lapply(1:19,function(i){
  x <- mk_plot[[i]]
  expr <- neu@assays$SCT@data[x,]%>%as.matrix()%>%t %>% as.data.frame()
  expr$group <- Idents(neu)
  f <- formula(paste0("cbind(",paste(x,collapse = ","),") ~ group"))
  re <- aggregate(f,expr,mean)
  re[,-1] <- apply(re[,-1],2,scale)
  re <- reshape2::melt(re)
  expr <- neu@assays$RNA@data[x,]%>%as.matrix()%>%t %>% as.data.frame()
  expr[expr>0] <- 1
  expr$group <- Idents(neu)
  re2 <- aggregate(f,expr,sum)
  re2[,-1] <- re2[,-1]/as.numeric(table(Idents(neu)))
  re2 <- reshape2::melt(re2)
  mk <- qread("NEW/Neurons/re02_markers_subclusters.qs")
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
  xlab("")+ylab("")


#02 miloR
neu[["sample"]] <- neu$orig.ident
sel_kp(neu, k = 30, p = 0.1, d = 30)
milore <- run_miloR(neu, k = 30, p = 0.1, d = 30)
p <- plotnodegraph(milore, neu); plot(p)
savePDF(p,"NEW/Neurons/p03.1_neurons_miloR.pdf",width = 5.6, height = 4)



#03. DEGs analysis in c01:GABAergic Neurons based on muscat
library(muscat)
neu1 <- subset(neu, subset = celltype == "c01:GABAergic Neurons")
#02.1 DEGs in c01:GABAergic Neurons
deg <- findDEG_muscat(neu1, sample = "orig.ident", group = "sample_group")
write.csv(deg, "NEW/Neurons/muscat_DEG_detect_GABAergic Neurons.csv", row.names = F)
#02.2 DEGs in different subclusters 
#02.2.1 up DEGs 
degup <- deg %>% filter(p_adj < 0.1 & (DSfc > 0.5 | DDfc > 0.3)) %>% pull(gene)
neu1 <- AddModuleScore(neu1, list(degup), name = "response_up")
Featplot(neu1, "response_up1", group = "sample_group", size = 0.5) -> p; plot(p)
savePDF(p,"NEW/Neurons/p07.1_response_to_treatment_up_UMAP_GABAergic Neurons.pdf",width = 5,height = 2.6)
ggplot(neu1@meta.data, aes(subcluster, response_up1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.3, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  scale_fill_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"NEW/Neurons/p07.1_response_to_treatment_up_violin_GABAergic Neurons.pdf",width = 7,height = 2.6)
model <- lmer(response_up1 ~ sample_group * subcluster + (1|orig.ident), data = neu1@meta.data)
stat <- treat_lmer(model)
openxlsx::write.xlsx(stat, "NEW/Neurons/p07.1_response_to_treatment_up_stat_GABAergic Neurons.xlsx",rowNames = F)
neu1[["shift_state"]] <- calc_shift(neu1, degup, group = "sample_group")
table(neu1$subcluster[neu1$sample_group=="T"], neu1$shift_state[neu1$sample_group=="T"]) %>% prop.table(margin = 1)
DimPlot(neu1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void() -> p; plot(p)
savePDF(p,"NEW/Neurons/p07.1_response_to_treatment_up_shifted_state_GABAergic Neurons.pdf",width = 4.8, height = 2.6)
#02.2.2 dn DEGs 
degdn <- deg %>% filter(p_adj<0.1 & (DSfc < (-0.5) | DDfc < (-0.3))) %>% pull(gene)
neu1 <- AddModuleScore(neu1, list(degdn), name = "response_dn")
Featplot(neu1, "response_dn1", group = "sample_group", size =0.8) ->p; plot(p)
savePDF(p,"NEW/Neurons/p07.2_response_to_treatment_dn_UMAP_GABAergic Neurons.pdf",width = 5,height = 2.6)
ggplot(neu1@meta.data, aes(subcluster, response_dn1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.1, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  scale_fill_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"NEW/Neurons/p07.2_response_to_treatment_dn_violin_GABAergic Neurons.pdf",width = 7,height = 2.6)
model <- lmer(response_dn1 ~ sample_group * subcluster + (1|orig.ident), data = neu1@meta.data)
stat <- treat_lmer(model)
openxlsx::write.xlsx(stat, "NEW/Neurons/p07.2_response_to_treatment_dn_stat_GABAergic Neurons.xlsx",rowNames = F)
neu1[["shift_state"]] <- calc_shift(neu1, degdn, orient = -1, group = "sample_group")
table(neu1$subcluster[neu1$sample_group=="T"], neu1$shift_state[neu1$sample_group=="T"]) %>% prop.table(margin = 1)
DimPlot(neu1, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"NEW/Neurons/p07.2_response_to_treatment_dn_shifted_state_GABAergic Neurons.pdf",width = 4.8,height = 2.6)



#04. DEGs analysis in c02:Glutamatergic Neurons based on muscat 
neu2 <- subset(neu, subset = celltype == "c02:Glutamatergic Neurons")
#02.1 DEGs in c02:Glutamatergic Neurons
deg <- findDEG_muscat(neu2, sample = "orig.ident", group = "sample_group")
write.csv(deg, "NEW/Neurons/muscat_DEG_detect_Glutamatergic Neurons.csv", row.names = F)
#02.2 DEGs in different subclusters 
#02.2.1 up DEGs 
degup <- deg %>% filter(p_adj < 0.1 & (DSfc > 0.4 | DDfc > 0.2)) %>% pull(gene)
neu2 <- AddModuleScore(neu2, list(degup), name = "response_up")
Featplot(neu2, "response_up1", group = "sample_group", size = 0.5) -> p; plot(p)
savePDF(p,"NEW/Neurons/p08.1_response_to_treatment_up_UMAP_Glutamatergic Neurons.pdf",width = 5,height = 2.6)
ggplot(neu2@meta.data, aes(subcluster, response_up1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.1, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  scale_fill_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"NEW/Neurons/p08.1_response_to_treatment_up_violin_Glutamatergic Neurons.pdf",width = 6,height = 2.6)
model <- lmer(response_up1 ~ sample_group * subcluster + (1|orig.ident), data = neu2@meta.data)
stat <- treat_lmer(model)
openxlsx::write.xlsx(stat, "NEW/Neurons/p08.1_response_to_treatment_up_stat_Glutamatergic Neurons.xlsx",rowNames = F)
neu2[["shift_state"]] <- calc_shift(neu2, degup, group = "sample_group")
table(neu2$subcluster[neu2$sample_group=="T"], neu2$shift_state[neu2$sample_group=="T"]) %>% prop.table(margin = 1)
DimPlot(neu2, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void() -> p; plot(p)
savePDF(p,"NEW/Neurons/p08.1_response_to_treatment_up_shifted_state_Glutamatergic Neurons.pdf",width = 4.8, height = 2.6)
#02.2.2 dn DEGs 
degdn <- deg %>% filter(p_adj<0.1 & (DSfc < (-0.4) | DDfc < (-0.2))) %>% pull(gene)
neu2 <- AddModuleScore(neu2, list(degdn), name = "response_dn")
Featplot(neu2, "response_dn1", group = "sample_group", size =0.8) ->p; plot(p)
savePDF(p,"NEW/Neurons/p08.2_response_to_treatment_dn_UMAP_Glutamatergic Neurons.pdf",width = 5,height = 2.6)
ggplot(neu2@meta.data, aes(subcluster, response_dn1, fill = sample_group))+
  geom_point(aes(group = sample_group), size = 0.1, alpha = 1, color = "#757575",
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(width = 0.7, alpha = 0.5, color = "#F1F1F1",linewidth = 0.3)+
  scale_color_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  scale_fill_manual(values = c(NC = "#FFB900", `T` = "#5773CC"))+
  theme_classic() ->p; plot(p)
savePDF(p,"NEW/Neurons/p08.2_response_to_treatment_dn_violin_Glutamatergic Neurons.pdf",width = 6,height = 2.6)
model <- lmer(response_dn1 ~ sample_group * subcluster + (1|orig.ident), data = neu2@meta.data)
stat <- treat_lmer(model)
openxlsx::write.xlsx(stat, "NEW/Neurons/p08.2_response_to_treatment_dn_stat_Glutamatergic Neurons.xlsx",rowNames = F)
neu2[["shift_state"]] <- calc_shift(neu2, degdn, orient = -1, group = "sample_group")
table(neu2$subcluster[neu2$sample_group=="T"], neu2$shift_state[neu2$sample_group=="T"]) %>% prop.table(margin = 1)
DimPlot(neu2, group.by = "shift_state", split.by = "sample_group", pt.size = 0.1) +
  scale_color_manual(values = c(stable = "#d7d7d7", shifted = "#1B1919")) +
  theme_void()->p; plot(p)
savePDF(p,"NEW/Neurons/p08.2_response_to_treatment_dn_shifted_state_Glutamatergic Neurons.pdf",width = 4.8,height = 2.6)





#02.1 density plot # dont use that!
# neu <- qread("NEW/Neurons/re01_SEUobj_neurons annotated.qs")
# umap <- data.frame(Embeddings(neu,reduction = "umap"),
#                    sample_group = neu$sample_group,
#                    cluster = neu$subcluster)
# umap[umap$sample_group=="NC","densi"] <- num_density(umap[umap$sample_group=="NC",1:2],dist=0.5)
# umap[umap$sample_group=="T","densi"] <- num_density(umap[umap$sample_group=="T",1:2],dist=0.5)
# q <- min(quantile(umap[umap$sample_group=="NC","densi"],0.99),
#          quantile(umap[umap$sample_group=="T","densi"],0.99))
# p <- lapply(c("NC","T"),function(x){
#   plotdata <- umap[umap$sample_group==x,]
#   plotdata$densi[plotdata$densi>q] <- q
#   p <- ggplot(plotdata,mapping = aes(UMAP_1,UMAP_2))+
#     geom_point_rast(aes(color=densi),shape=".",raster.dpi = 600)+
#     scale_color_paletteer_c("grDevices::Grays",direction = -1)+
#     theme_void()+
#     xlab("")+ylab("")
# })
# p <- wrap_plots(p,nrow = 2,guides = "collect"); plot(p)
# savePDF(p,"NEW/Neurons/p03_neurons_density.pdf",width = 4,height = 5)
#02.2 density chiqtest # dont use that
# chitable <- table(umap$sample_group,umap$cluster)
# chire <- chisq.test(chitable)
# p <- data.frame(re = chire$stdres[1,], cluster = names(chire$stdres[1,])) %>% arrange(re)
# p$cluster <- factor(p$cluster,levels = p$cluster)
# col <- c(paletteer_c("grDevices::Grays", 18)[c(2,4:14)], paletteer_c("grDevices::Reds 3", 12)[3:9])
# names(col) <- lvl
# ggplot(p,aes(cluster,-re,fill=cluster))+
#   geom_col(position = "dodge")+
#   scale_y_continuous(limits = c(-25,25))+
#   scale_fill_manual(values = col)+
#   theme_minimal()+
#   theme(panel.grid.major.x = element_blank())+
#   labs(title = "X2test standardized residuals", x=NULL, y=NULL) -> p; plot(p)
# savePDF(p,"NEW/Neurons/p03_neurons_x2test.pdf",width = 4,height = 3)


#03.1 DEG calculate --> in different sample group  #dont use that!
# neus <- SplitObject(neu)
# DEG <- pblapply(lvl, function(x){
#   n <- neus[[x]]
#   Idents(n) <- factor(n$sample_group)
#   m <- FindMarkers(n,ident.1 = "NC",ident.2 = "T",recorrect_umi=F,logfc.threshold = 0,assay = "SCT") %>%
#     mutate(comparsion = "NC vs T",
#            p_val_adj = p.adjust(p_val,method = "BH"))
#   m$gene <- rownames(m);m
# })
# names(DEG) <- gsub(":","_",lvl)
# openxlsx::write.xlsx(DEG,"NEW/Neurons/stat01_DEGs_SampleLevel_allSubclusters.xlsx")
# names(DEG) <- lvl
# qsave(DEG,"NEW/Neurons/re03_DEGs_SampleLevel.qs")
# #03.2 DEG map
# umap <- data.frame(Embeddings(neu,reduction = "umap"),
#                    cluster = neu$subcluster)
# degmap <- lapply(c("NC vs T"), function(x){
#   re <- lapply(lvl,function(l){
#     n_deg <- sum(DEG[[l]]$p_val_adj<0.05 & DEG[[l]]$avg_log2FC>0.1)
#     if(n_deg < 50){n_deg <- 50}
#     # if(n_deg > 150){n_deg <- 150}
#     umap <- umap[umap$cluster==l,] %>%
#       mutate(ndeg = n_deg)
#   })%>%purrr::reduce(rbind)
#   # x1 <- re[!duplicated(re$cluster),]
#   ggplot(re,mapping = aes(UMAP_1,UMAP_2))+
#     geom_scattermore(aes(color=ndeg),pixels = c(1024,1024),pointsize = 3)+
#     scale_color_gradient(low = "grey",high = "#F72735")+
#     theme_classic()+
#     theme(axis.ticks = element_blank(),
#           axis.line = element_blank(),
#           axis.text = element_blank(),
#           panel.grid = element_blank())+
#     xlab("")+ylab("") -> p
# })
# plist <- patchwork::wrap_plots(degmap,ncol = 1,guides = "collect"); plist
# savePDF(plist,"NEW/Neurons/p04_DEG_umap.pdf",width = 4,height = 3)
# #03.3 enrichment analysis and plot
# ORA <- list.files("NEW/Neurons/",pattern = "^MS",full.names = T)
# ORA <- lapply(ORA,function(f){
#   re <- readxl::read_xlsx(f,sheet = 2)
#   idy <- re$sel %in% "YES"
#   re <- re[idy,c(2,3,4,5,8,9)] %>% 
#     mutate(Description = stringr::str_to_sentence(gsub("- M.*","",Description)),
#            subtpye = gsub(".*MS_(.*)_.*.xlsx","\\1",f),
#            # vs = gsub(".*[Agrp|Pomc]_(.*)_.*","\\1",f),
#            ori = gsub(".*[AGRP|POMC]_(.*)\\.xlsx","\\1",f),
#            count = gsub("/.*","",InTerm_InList)%>%as.numeric())
#   if(nrow(re)>7){re <- re[1:7,]} ; re
# })%>%purrr::reduce(rbind)
# op <- ORA[ORA$subtpye=="AGRP",] %>% 
#   group_by(paste(ori,Description)) %>% 
#   summarise(Description = unique(Description),
#             LogP = mean(LogP),
#             ori = unique(ori),
#             # vs = paste(vs,collapse = ";"),
#             count = mean(count)) %>% 
#     arrange(LogP) %>%
#     mutate(Description = factor(Description,levels = unique(Description)%>%rev))
# ggplot(op)+
#   geom_segment(aes(x=0,xend=count,y=Description,yend=Description,color=paste(ori)),
#                size=8,alpha=0.8)+
#   geom_point(aes(x=count,y=Description,size=(-LogP),color=paste(ori)),shape=1)+
#   scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
#   scale_size_continuous(range = c(2,7))+
#   scale_x_continuous(expand = c(0,0),limits = c(0,35))+
#   facet_grid(ori~.,scales = "free_y",space="free_y")+
#   coord_cartesian(clip = "off")+
#   theme_classic()+
#   theme(panel.spacing.y = unit(0.2,"cm"),
#         axis.text.y = element_text(margin = margin(r=-5,unit = "cm"),
#                                    color="black",hjust = 1),
#         axis.ticks.y = element_blank())+
#   ylab("")+xlab("") -> p; plot(p)
# savePDF(p,"NEW/Neurons/p05_ORA_AGRP.pdf",width = 5.5,height = 5)
# rm(op)
# #03.4 DEG volcano plot
# vc <- DEG$`c02:Pomc+` %>%
#   filter(comparsion %in% c("NC vs T")) %>%
#   mutate(avg_log2FC = -avg_log2FC)
# vc$avg_log2FC[vc$avg_log2FC>1.05] <- 1.05
# vc1 <- vc[abs(vc$avg_log2FC)>0.1&vc$p_val_adj<0.05,] %>% 
#   mutate(col = paste(comparsion,sign(avg_log2FC)))
# ggplot(vc,aes(avg_log2FC,-log10(p_val_adj),shape=comparsion))+
#   geom_point(color="#A1A1A1")+
#   geom_point(data = vc1,aes(color=col,size=abs(avg_log2FC)))+
#   # scale_color_manual(values = c("#A852CC","#6B990F","#0F6B99","#FF2B00"))+
#   scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
#   scale_size_continuous(range = c(1,3))+
#   scale_x_continuous(limits = c(-0.6,0.6),expand = c(0,0))+
#   geom_hline(yintercept = 1.31,lty="dashed")+
#   geom_vline(xintercept = c(-0.1,0.1),lty="dashed")+
#   theme_classic()+
#   xlab("Log2(Fold Change)")+ylab("-Log10(adj.P)")->p; plot(p)
# savePDF(p,"NEW/Neurons/p05_DEGvol_POMC.pdf",width = 4.5,height = 2.5)

#04 Fgfr2 expression #Dont use that!
# neu1 <- neu[,neu$subcluster=="c01:Agrp+"&neu@reductions$umap@cell.embeddings[,2]<(-10)]
# Geneplot(neu1[,neu1$sample_group=="NC"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
# savePDF(p,"NEW/Neurons/p06_Agrp_Fgfr2_NC.pdf",width = 3,height = 2.4)
# Geneplot(neu1[,neu1$sample_group=="T"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
# savePDF(p,"NEW/Neurons/p06_Agrp_Fgfr2_T.pdf",width = 3,height = 2.4)
# neu1 <- neu[,neu$subcluster=="c02:Pomc+"&
#               neu@reductions$umap@cell.embeddings[,1]<(0) & neu@reductions$umap@cell.embeddings[,1]>(-10) &
#               neu@reductions$umap@cell.embeddings[,2]<(-2) & neu@reductions$umap@cell.embeddings[,2]>(-10)]
# Geneplot(neu1[,neu1$sample_group=="NC"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
# savePDF(p,"NEW/Neurons/p06_Pomc_Fgfr2_NC.pdf",width = 3,height = 2.4)
# Geneplot(neu1[,neu1$sample_group=="T"],"Fgfr2",size = 0.75,raster = F)[[1]] -> p; plot(p)
# savePDF(p,"NEW/Neurons/p06_Pomc_Fgfr2_T.pdf",width = 3,height = 2.4)
# 
# umap <- data.frame(Embeddings(neu1,reduction = "umap"),
#                    cluster = neu1$subcluster,
#                    sample_group = neu1$sample_group,
#                    value = neu1@assays$SCT@data["Fgfr2",])
# ggplot(umap,aes("c01:Agrp+",value,fill = sample_group))+
#   geom_jitter(aes(color = sample_group),alpha = 0.6,size=0.5,
#               position = position_jitterdodge(jitter.width = 0.3,dodge.width = 0.75))+
#   geom_violin(width = 0.7,scale = "area",alpha = 0.7,adjust=1)+
#   scale_y_continuous(limits = c(0,3.2))+
#   scale_color_manual(values = c("#fd8d3c","#4d8f87"))+
#   scale_fill_manual(values = c("#fd8d3c","#4d8f87"))+
#   theme_bw() -> p; plot(p)
# savePDF(p,"NEW/Neurons/p06_Agrp_Fgfr2_violin.pdf",width = 3.5,height = 2.2)





df <- data.frame(expr = seu@assays$RNA@data["Fgfr2",],
                 sample_group = seu$sample_group,
                 celltype = seu$celltype)
ggplot(df, aes(celltype, expr))+
  geom_point(aes(color = sample_group), size = 0.2, alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_violin(aes(fill = sample_group), width = 0.7, alpha = 0.3, 
              color = "#616161",linewidth = 0.4, adjust = 2)+
  scale_color_manual(values = c(yes = "#CD534C", no = "#79AF97"))+
  scale_fill_manual(values = c(yes = "#F1F1F1", no = "#F1F1F1"))+
  theme_classic()+
  labs(x = NULL, y = "nFeature_RNA", title = "Masking rate = 0.8")-> p; plot(p)


df <- data.frame(expr = neu@assays$RNA@data["Fgfr2",],
                 sample_group = neu$sample_group,
                 celltype = neu$subcluster)
ggplot(df, aes(celltype, expr))+
  geom_point(aes(color = sample_group), size = 0.2, alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.24, dodge.width = 0.7)) +
  geom_boxplot(aes(fill = sample_group))+
  geom_violin(aes(fill = sample_group), width = 0.7, alpha = 0.3, 
              color = "#616161",linewidth = 0.4, adjust = 2)+
  scale_color_manual(values = c(yes = "#CD534C", no = "#79AF97"))+
  scale_fill_manual(values = c(yes = "#F1F1F1", no = "#F1F1F1"))+
  theme_classic()+
  labs(x = NULL, y = "nFeature_RNA", title = "Masking rate = 0.8")-> p; plot(p)











