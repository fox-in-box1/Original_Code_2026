setwd("/home/disk1/wgs")
library(CellChat)
library(patchwork)
library(Seurat)
library(qs)
library(ComplexHeatmap)
source("NEW/FUN02_run_cc.R")
source("NEW/useful_functions.R")


##################
## data prepare ##
##################
seu <- qread("NEW/s03_Integration Data mapped k5.qs")
neu <- qread("NEW/Neurons/re01_SEUobj_neurons annotated.qs")
lvl <- unique(neu$subcluster) %>% .[grep("^c0",.)]
tany <- qread("NEW/Tanycytes/re01_Tanseu_Ig.qs")
cell <- rbind(neu@meta.data[,c("orig.ident","subcluster")],
              tany@meta.data[,c("orig.ident","subcluster")])
seu <- AddMetaData(seu,cell)
seu$subcluster[is.na(seu$subcluster)] <- seu$celltype[which(is.na(seu$subcluster))]


#######################
## cellchat analysis ##
#######################
#data for cellchat
CellChatDB <- CellChatDB.mouse
ccobj <- lapply(c("NC","T"),function(x){
  idy <- seu$subcluster %in% c("A1","A2","B1","B2","c01:Agrp+","c02:Pomc+") & seu$sample_group %in% x
  sdata <- seu[,idy]
  meta <- data.frame(labels = sdata$subcluster, row.names = colnames(sdata))
  cc <- run_cellchat(sdata@assays$RNA@data,meta) 
})
names(ccobj) <- c("NC","T")
col <- c("#D2854A","#EDD9A9","#91A17B","#486349","#373737FF","#B01526FF")
qsave(ccobj,"NEW/cellchat/re01_ccobj.qs")
ccobj <- qread("NEW/cellchat/re01_ccobj.qs")
cellchat <- mergeCellChat(ccobj, add.names = names(ccobj))

#plot01 interactions overview
weight.max <- getMaxWeight(ccobj, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(ccobj)) {
  netVisual_circle(ccobj[[i]]@net$weight, weight.scale = T, label.edge= F, 
                   color.use = col,
                   edge.weight.max = weight.max[2]/2, edge.width.max = 12, 
                   title.name = paste0("Weight of interactions - ", names(ccobj)[i]))
}
#landscape 6*4 plot/p01_overview_weigt.pdf
weight.max <- getMaxWeight(ccobj, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(ccobj)) {
  netVisual_circle(ccobj[[i]]@net$count, weight.scale = T, label.edge= F, 
                   color.use = col,
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Count of interactions - ", names(ccobj)[i]))
}
#landscape 6*4 plot/p01_overview_count.pdf

#plot02 outgoing & incoming overview and differences
# combining all the identified signaling pathways from different datasets 
pathway.union <- c(ccobj[[1]]@netP$pathways, ccobj[[2]]@netP$pathways) %>% unique
ht1 = netAnalysis_signalingRole_heatmap(ccobj[[1]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(ccobj)[1], width = 3, height = 6,color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(ccobj[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(ccobj)[2], width = 3, height = 6,color.heatmap = "OrRd")
# ht3 = ht2@matrix-ht1@matrix
# ht3 = netAnalysis_signalingRole_heatmap(ccobj[[3]], pattern = "outgoing", signaling = pathway.union,
#                                         title = names(ccobj)[3], width = 3, height = 6,color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) -> p
savePDF(p,"NEW/cellchat/p02_outgoing.pdf",width = 6,height = 4)

ht1 = netAnalysis_signalingRole_heatmap(ccobj[[1]], pattern = "incoming", signaling = pathway.union,
                                        title = names(ccobj)[1], width = 3, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(ccobj[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(ccobj)[2], width = 3, height = 6)
# ht3 = netAnalysis_signalingRole_heatmap(ccobj[[3]], pattern = "incoming", signaling = pathway.union,
#                                         title = names(ccobj)[3], width = 3, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) -> p
savePDF(p,"NEW/cellchat/p02_incoming.pdf",width = 6,height = 4)

#differences:01_outgoing/incoming differences
num.link <- sapply(ccobj, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(ccobj)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(ccobj[[i]], title = names(ccobj)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
data <- base::rbind(gg[[1]]$data %>% mutate(sg = "NC"),
                    gg[[2]]$data %>% mutate(sg = "T"))
ggplot(data[data$sg%in%c("NC","T")&data$labels%in%c("c01:Agrp+","c02:Pomc+"),],aes(x,y))+
  geom_point(aes(color=labels,size=Count,shape=sg))+
  scale_color_manual(values = col[5:6])+
  scale_size_continuous(range = c(3,6))+
  scale_shape_manual(values = c(20,1))+
  theme_classic()+
  ylab("Incoming interaction strength")+xlab("Outgoing interaction strength")->p; plot(p)
savePDF(p,"NEW/cellchat/p03_celltypes_overall_diff1.pdf",width = 4.2,height = 3)
ggplot(data[data$sg%in%c("NC","T")&data$labels%in%c("A1","A2","B1","B2"),],aes(x,y))+
  geom_point(aes(color=labels,size=Count,shape=sg))+
  scale_color_manual(values = col)+
  scale_size_continuous(range = c(3,6))+
  scale_shape_manual(values = c(20,1))+
  theme_classic()+
  ylab("Incoming interaction strength")+xlab("Outgoing interaction strength")->p; plot(p)
savePDF(p,"NEW/cellchat/p03_celltypes_overall_diff2.pdf",width = 4.2,height = 3)
rm(num.link,weight.MinMax,gg,data)
#landscape 4.2*3

cellchat <- mergeCellChat(ccobj, add.names = names(ccobj))
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2))->p; plot(p)
savePDF(p,"NEW/cellchat/p02_information_flow.pdf",width = 3.5,height = 5)
# rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,3))
#portait 5*3.5

#differences:02_signaling differences in neurons (Pathway level)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "c01:Agrp+");plot(gg1)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "c02:Pomc+");plot(gg2)
data <- rbind(gg1$data %>% mutate(celltype="Agrp"),
              gg2$data %>% mutate(celltype="Pomc"))
data$max <- apply(data[,1:2],1,function(x) x[which.max(abs(x))])
data <- data %>% arrange(celltype,max)
data2 <- data[c(1:6,19:23),]
data2[data2$outgoing<(-0.02),"outgoing"] <- (-0.021)
data2[data2$incoming<(-0.02),"incoming"] <- (-0.021)
ggplot(data2,aes(outgoing,incoming))+
  geom_point(aes(shape=celltype),size=2,color="#C1C1C1")+
  geom_point(data=data2,aes(color=celltype,shape=celltype),size=4)+
  scale_color_manual(values = c("#373737FF","#B01526FF"))+
  scale_shape_manual(values = c(20,17))+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_vline(xintercept = 0,lty="dashed")+
  theme_classic()+
  ylab("Differential incoming strength")+xlab("Differential outgoing strength")+
  ggrepel::geom_text_repel(data = data2, 
                           mapping = aes(label = labels, colour = celltype), 
                           size = 3, show.legend = F, segment.size = 0.2, 
                           segment.alpha = 0.5)->p; plot(p)
savePDF(p,"NEW/cellchat/p04_signing_diff_Neurons.pdf",width = 4.2,height = 2.8)

#differences:02_signaling differences in tanycytes (Pathway level)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B2",comparison = c(1,2));plot(gg1)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B1",comparison = c(1,2));plot(gg2)
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "A1",comparison = c(1,2));plot(gg3)
data <- rbind(gg1$data %>% mutate(celltype="B2"),
              gg2$data %>% mutate(celltype="B1"),
              gg3$data %>% mutate(celltype="A1"))
data$max <- apply(data[,1:2],1,function(x) x[which.max(abs(x))])
data <- data %>% arrange(celltype,max)
data2 <- data[c(1:3,12:14,23:25),]
# data2[data2$outgoing<(-0.02),"outgoing"] <- (-0.021)
# data2[data2$incoming<(-0.02),"incoming"] <- (-0.021)
ggplot(data2,aes(outgoing,incoming))+
  geom_point(aes(shape=celltype),size=2,color="#C1C1C1")+
  geom_point(data=data2,aes(color=celltype,shape=celltype),size=4)+
  scale_color_manual(values = c("#D2854A","#91A17B","#486349"))+
  scale_shape_manual(values = c(18,15,5))+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_vline(xintercept = 0,lty="dashed")+
  theme_classic()+
  ylab("Differential incoming strength")+xlab("Differential outgoing strength")+
  ggrepel::geom_text_repel(data = data2, 
                           mapping = aes(label = labels, colour = celltype), 
                           size = 3, show.legend = F, segment.size = 0.2, 
                           segment.alpha = 0.5)->p; plot(p)
savePDF(p,"NEW/cellchat/p04_signing_diff_tanycytes.pdf",width = 4.2,height = 2.8)


#Differential L-R pairs
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:6),  comparison = c(1,2), angle.x = 45)->p1
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(5:6),  comparison = c(1,2), angle.x = 45)->p2
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(5:6),  comparison = c(1,2), angle.x = 45)->p3
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:6),  comparison = c(1,2), angle.x = 45)->p4
p1 + p2 + p3 + p4 + plot_layout(guides = "collect", nrow = 1) -> plist1; plot(plist1)
savePDF(plist1, "NEW/cellchat/p04.1_Tanycyte_to_Neuron_signaling_dot.pdf", width = 9, height = 4)
netVisual_bubble(cellchat, sources.use = c(5:6), targets.use = 1,  comparison = c(1,2), angle.x = 45)->p5
netVisual_bubble(cellchat, sources.use = c(5:6), targets.use = 2,  comparison = c(1,2), angle.x = 45)->p6
netVisual_bubble(cellchat, sources.use = c(5:6), targets.use = 3,  comparison = c(1,2), angle.x = 45)->p7
netVisual_bubble(cellchat, sources.use = c(5:6), targets.use = 4,  comparison = c(1,2), angle.x = 45)->p8
p5 + p6 + p7 + p8 + plot_layout(guides = "collect", nrow = 1) -> plist2; plot(plist2)
savePDF(plist2, "NEW/cellchat/p04.2_Neuron_to_Tanycyte_signaling_dot.pdf", width = 9, height = 4)

difflr1 <- rbind(p1$data,p2$data,p3$data,p4$data) %>% na.omit()
difflr1 <- reshape2::dcast(difflr1, interaction_name_2 + pathway_name + annotation + source + target + ligand + receptor ~ dataset, value.var = "prob.original")
difflr1[is.na(difflr1)] <- 0
difflr1$diff_prob <- difflr1$`T` - difflr1$NC
difflr2 <- rbind(p5$data,p6$data,p7$data,p8$data) %>% na.omit()
difflr2 <- reshape2::dcast(difflr2, interaction_name_2 + pathway_name + annotation + source + target + ligand + receptor ~ dataset, value.var = "prob.original")
difflr2[is.na(difflr2)] <- 0
difflr2$diff_prob <- difflr2$`T` - difflr2$NC
difflr <- rbind(difflr1,difflr2) %>% arrange(desc(abs(diff_prob)))
write.csv(difflr,"NEW/cellchat/p04.1_04.2_diff_LRs.csv", row.names = F)



#all celltypes chord plot
#选择Tanycytes与Neurons间通信概率在T组和NC组之间差异最大的前10个受配体对
#Select the top 10 ligand-receptor pairs showing the largest changes in communication probability between Tanycytes and Neurons across different feeding condition groups.
#为了评估这10对配体-受体的特异性，我们进一步检测了它们在所有细胞类型中的表达模式及通讯概率。
#01. expression patterns
extexpr_summary <- function(seu.obj, gene){
  df <- data.frame(expr = seu.obj@assays[["RNA"]]@data[gene,],
                   celltype = seu.obj$celltype) %>%
    mutate(expr_b = case_when(expr > 0 ~ 1, T ~ 0)) %>%
    group_by(celltype) %>%
    summarise(mexpr = mean(expr), prop = sum(expr_b)*100 / n()) %>%
    mutate(rexpr = scale(mexpr), gene = gene)
}
sel_pair <- difflr %>% filter(!duplicated(interaction_name_2)) %>% .[1:10,]
df <- lapply(unique(sel_pair$pathway_name),function(x){
  subpair <-  sel_pair %>% filter(pathway_name == x)
  genes_l <- unique(subpair$ligand)
  genes_r <- unique(subpair$receptor)
  df1 <- lapply(genes_l, extexpr_summary, seu.obj = seu) %>% purrr::reduce(rbind) %>% 
    mutate(pathway_name = x, type = "ligand")
  df2 <- lapply(genes_r, extexpr_summary, seu.obj = seu) %>% purrr::reduce(rbind) %>% 
    mutate(pathway_name = x, type = "receptor")
  df <- rbind(df1,df2) %>% arrange(type)
}) %>% purrr::reduce(rbind) %>% filter(!duplicated(paste(gene, celltype))) 
df$gene <- factor(df$gene, levels = unique(df$gene))
n_genes <- length(levels(df$gene))
grid_breaks <- c(2.5, 6.5, 9.5, 10.5, 11.5, 13.5)
ggplot(df, aes(x = as.numeric(gene), y = celltype))+
  geom_vline(xintercept = grid_breaks, lty = "dashed")+ 
  geom_point(aes(size = prop, color = rexpr))+
  scale_color_paletteer_c("ggthemes::Gold-Purple Diverging", direction = 1)+
  scale_x_continuous(breaks = 1:n_genes, labels = levels(df$gene))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = NULL, y = NULL) -> p; plot(p)
savePDF(p, "NEW/cellchat/p05_Diff_LRs_heatmap.pdf", width = 9, height = 3)
genes <- levels(df$gene)
plist <- Geneplot(seu, genes, q = 0.99, size = 0.01, nrow = 3) + plot_layout(guides = "collect"); plot(plist)
savePDF(plist, "NEW/cellchat/p05_Diff_LRs_allcells_UMAP.pdf", width = 7, height = 4)

#02. communication prob.
meta <- data.frame(labels = seu$celltype, row.names = colnames(seu))
cc <- run_cellchat(seu@assays$RNA@data, meta) 
sel_pair <- difflr$interaction_name_2 %>% as.character() %>% unique() %>% .[1:10]
schat <- netVisual_bubble(cc, sources.use = 1:10, targets.use = 1:10, angle.x = 45, return.data=T)$communication %>%
  filter(interaction_name_2 %in% sel_pair)
#prepare data
cell_col <- paletteer_d("ggsci::default_jco")[c(3,4,1,6,8,2,9,7,5,10)]
names(cell_col) <- levels(Idents(seu))
df1 <- schat %>% na.omit() %>%
  mutate(ligand = paste(ligand, source, sep = "_"), receptor = paste(receptor, target, sep = "_"),
         value = abs(prob)) %>%
  select(ligand, receptor, value, source, target) %>%
  mutate(linecol = adjustcolor(cell_col[match(source,names(cell_col))], alpha = 0.5))
df1$value[df1$value<0.25] <- 0.05
# -> script06_plot_chord.R



#Fgf10 - Fgfr2 signal in all celltypes
seu$subcluster[seu$celltype%in%c("c06:Tanycytes")] <- "c06:Tanycytes"
ccobj <- lapply(c("NC","T"),function(x){
  idy <- seu$sample_group %in% x
  sdata <- seu[,idy]
  meta <- data.frame(labels = sdata$subcluster, row.names = colnames(sdata))
  cc <- run_cellchat(sdata@assays$RNA@data,meta) 
})
names(ccobj) <- c("NC","T")
cellchat <- mergeCellChat(ccobj, add.names = names(ccobj))
LR <- netVisual_bubble(cellchat, sources.use = 23, targets.use = c(1:23),  comparison = c(1,2),return.data=T)$communication %>%
  filter(interaction_name_2 %in% c("Fgf10  - Fgfr2"))
LR <- reshape2::dcast(LR, interaction_name_2 + pathway_name + annotation + source + target + ligand + receptor ~ dataset, value.var = "prob.original") %>%
  mutate(diff = NC - `T`) %>% arrange(desc(diff))
write.csv(LR, "NEW/cellchat/Fgf10_Fgfr2_all_celltype_diff.csv", row.names = F)

