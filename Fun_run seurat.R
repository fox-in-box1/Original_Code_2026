
#function2: automatic QC process --> Part2 Visualization function 5: QC visualization
qc_process  <- function(seuobj, islist = T){
  if(islist){
    qc <- lapply(seuobj, function(x){
      data.frame(sample = unique(x$sample), ncell = ncol(x), 
                 nCounts = median(x$nCount_RNA), nFeature = median(x$nFeature_RNA),
                 percent_mt_0.95 = quantile(x$percent_mt,0.95), 
                 percent_ribo_0.05 = quantile(x$percent_ribo,0.05),
                 percent_hb_0.99 = quantile(x$percent_hb,0.99))
    }) %>% purrr::reduce(rbind)
  }else{
    qc <- seuobj@meta.data %>%
      group_by(sample) %>%
      summarise(ncell = n(), nCounts = median(nCount_RNA), nFeature = median(nFeature_RNA),
                percent_mt_0.95 = quantile(percent_mt,0.95), 
                percent_ribo_0.05 = quantile(percent_ribo,0.05),
                percent_hb_0.99 = quantile(percent_hb,0.99))
  }
  return(qc)
}


#function5: QC visualization
#colored by sample_group. 
plot_qc_results <- function(seu, islist = T, feats = NULL, cols = "default"){
  require(ggrastr)
  if(islist){
    df <- purrr::reduce(lapply(seu,function(x)x@meta.data), rbind)
  }else{
    df <- seu@meta.data
  }
  if(!"sample" %in% colnames(df)) stop("Error: There is no column named 'sample'")
  if(!"sample_group" %in% colnames(df)) stop("Error: There is no column named 'sample_group'")
  df[["sample"]] <- paste(df[["sample_group"]], df[["sample"]], sep = "_")
  if(cols == "default"){
    cols <- c("#FFB900", "#5773CC", "#CD534C", "#868686", "#79AF97")
  }
  if(is.null(feats)){
    feats <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb")
    feats <- feats[feats %in% colnames(df)]
    message(paste0(paste(feats,collapse = "; ")," exist in seurat object"))
  }
  p.list <- lapply(feats,function(f){
    df <- df[,c(f,"sample_group","sample")]
    ggplot(df, aes(sample, df[[f]], fill = sample_group))+
      geom_jitter_rast(aes(color = sample_group), alpha = 0.03, shape=16, size=0.4)+
      geom_violin(scale="width", alpha=0.7, adjust=2)+
      scale_color_manual(values = cols)+
      scale_fill_manual(values = cols)+
      theme_classic()+
      labs(title = f,y=NULL)
  })
  patchwork::wrap_plots(p.list, nrow = ceiling(length(feats)/5)) + 
    patchwork::plot_layout(guide="collect")
}
# ----------- Part1 end -----------------



# ----------- Part2 Visualization -----------------
# part1: Visualization

Geneplot <- function(seu_obj, gene, assay="RNA", reduction="umap", 
                     binary = F, q = 0.95, group = NULL, order = T,
                     raster = T, size = 1, stroke = NULL, nrow = 1){
  require(scattermore)
  require(ggsci)
  require(ggplot2)
  require(paletteer)
  require(patchwork)
  idy <- gene %in% rownames(seu_obj@assays[[assay]]@data)
  if(sum(idy)==0) {message("!!!!!Error: No genes were detected in the object"); return(NULL)}
  if(any(idy == F)) message(paste0("!!!Warning: ", paste(gene[!idy],collapse = " ~ "), "were not detected in the object"))
  gene <- gene[idy]
  p <- lapply(gene,function(x){
    message(paste0("PLOT: ",x))
    umap <- Embeddings(seu_obj,reduction = reduction) %>% as.data.frame()
    colnames(umap) <- c("UMAP_1","UMAP_2")
    umap$expr <- seu_obj@assays[[assay]]@data[x,]
    if(binary) umap$expr <- ifelse(umap$expr>0, 1, 0)
    if(!is.null(group)){
      if(length(group)==1) umap$g <- unlist(seu_obj[[group]]) else umap$g <- group
    }
    if(is.null(stroke)){
      if(nrow(umap)>20000) stroke <- 0 else stroke <- 0.4
    }
    if(order) umap <- umap[order(umap$expr,decreasing = F),]
    q <- quantile(umap$expr[umap$expr>0],q)
    umap$expr[umap$expr>q] <- q
    if(raster){
      require(ggrastr)
      ggplot(umap, aes(UMAP_1, UMAP_2))+
        geom_point_rast(data = umap[umap$expr<=0,], color = "#e1e1e1", raster.dpi = 300, 
                        stroke = stroke, size = size, alpha =0.7)+
        geom_point_rast(data = umap[umap$expr>0,], aes(color = expr), raster.dpi = 300, 
                        stroke = stroke, size = size)+
        scale_color_gradientn(values = c(0,0.1,0.3,0.5,0.7,0.9,1),
                              colors = paletteer_c("grDevices::OrRd",30)[c(28,24,18,12,6,1)]) -> p
    }else{
      ggplot(umap, aes(UMAP_1, UMAP_2))+
        geom_point(data = umap[umap$expr<=0,], color = "#c1c1c1",  
                   shape = 20, size = size)+
        geom_point(data = umap[umap$expr>0,], aes(color = expr), 
                   shape = 20, size = size)+
        scale_color_gradientn(values = c(0,0.1,0.3,0.5,0.7,0.9,1),
                              colors = paletteer_c("grDevices::OrRd",30)[c(28,24,18,12,6,1)]) -> p
    }
    message(paste0("DONE"))
    if(!is.null(group)){
      p <- p +
        facet_wrap(facets = . ~ g, nrow = 1)+
        theme(strip.background = element_rect(fill = "#f0f0f0"))
    }
    p <- p +
      theme_void()+
      # theme_bw()+
      theme(legend.key.width = unit(0.5,"cm"))+
      labs(x=NULL,y=NULL,title = x)
  })
  plist <- wrap_plots(p, nrow = nrow)
  return(plist)
}


Featplot <- function(seu_obj, feat, reduction="umap", 
                     binary = F, q = 0.99, group = NULL,
                     raster = T, size = 0.1, nrow = 1){
  require(ggsci)
  require(ggplot2)
  require(paletteer)
  require(patchwork)
  p <- lapply(feat,function(x){
    message(paste0("PLOT: ",x))
    umap <- Embeddings(seu_obj,reduction = reduction) %>% as.data.frame()
    colnames(umap) <- c("UMAP_1","UMAP_2")
    umap$expr <- unlist(seu_obj[[x]])
    if(binary) umap$expr <- ifelse(umap$expr>0, 1, 0)
    if(!is.null(group)){
      if(length(group)==1) umap$g <- unlist(seu_obj[[group]]) else umap$g <- group
    }
    if(nrow(umap)<20000) stroke = 0.1 else stroke = 0
    # umap <- umap[order(umap$expr,decreasing = F),]
    if(!is.null(q)){
      q1 <- quantile(umap$expr[umap$expr>0],q)
      q2 <- quantile(umap$expr[umap$expr<0],1-q)
      umap$expr[umap$expr>q1] <- q1
      umap$expr[umap$expr<q2] <- q2
    }
    umap <- umap %>% arrange(expr)
    if(raster){
      require(ggrastr)
      ggplot(umap, aes(UMAP_1, UMAP_2))+
        geom_point_rast(aes(color = expr), raster.dpi = 300, 
                        size = size, stroke = stroke)+
        scale_color_paletteer_c("grDevices::Blue-Red 2") -> p
    }else{
      ggplot(umap, aes(UMAP_1, UMAP_2))+
        geom_point(aes(color = expr), size = size, stroke = stroke)+
        scale_color_paletteer_c("grDevices::Blue-Red 2") -> p
    }
    message(paste0("DONE"))
    if(!is.null(group)){
      p <- p +
        facet_wrap(facets = . ~ g, nrow = 1)+
        theme(strip.background = element_rect(fill = "#f0f0f0"))
    }
    p <- p +
      theme_void()+
      theme(legend.key.width = unit(0.5,"cm"))+
      labs(x=NULL,y=NULL,title = x)
  })
  plist <- wrap_plots(p, nrow = nrow)
  return(plist)
}

Groupplot <- function(seu_obj, feat, reduction="umap", group = NULL,
                      raster = T, size = 1, stroke = NULL, nrow = 1){
  require(ggsci)
  require(ggplot2)
  require(paletteer)
  require(patchwork)
  p <- lapply(feat,function(x){
    message(paste0("PLOT: ",x))
    umap <- Embeddings(seu_obj,reduction = reduction) %>% as.data.frame()
    colnames(umap) <- c("UMAP_1","UMAP_2")
    umap$expr <- unlist(seu_obj[[x]])
    if(!is.null(group)){
      if(length(group)==1) umap$g <- unlist(seu_obj[[group]]) else umap$g <- group
    }
    if (is.null(stroke)) {
      if(nrow(umap)>20000) stroke = 0 else stroke = 0.25
    }
    
    if(raster){
      require(ggrastr)
      ggplot(umap, aes(UMAP_1, UMAP_2))+
        geom_point_rast(aes(color = expr), raster.dpi = 300, 
                        stroke = stroke, size = size)+
        scale_color_paletteer_d("ggsci::default_igv") -> p
    }else{
      ggplot(umap, aes(UMAP_1, UMAP_2))+
        geom_point(aes(color = expr), shape = shape, size = size)+
        scale_color_paletteer_d("ggsci::default_igv") -> p
    }
    message(paste0("DONE"))
    if(!is.null(group)){
      p <- p +
        facet_wrap(facets = . ~ g, nrow = nrow)+
        theme(strip.background = element_rect(fill = "#f0f0f0"))
    }
    p <- p +
      theme_void()+
      theme(legend.key.width = unit(0.5,"cm"))+
      labs(x=NULL,y=NULL,title = x)
  })
  plist <- wrap_plots(p, nrow = nrow)
  return(plist)
}
# ----------- Part2 end -----------------




# ----------- Part3 useful utils -----------------
#f01: char_replace
# replace characters one by one
# the length of 'from' need to equal to that of 'to'
char_replace <- function(vec,from,to){
  vec <- as.character(vec)
  if(length(from)!=length(to)){
    message("input is wrong"); return(NULL)
  }
  idy <- lapply(from,function(x){which(vec==x)})
  for(i in 1:length(from)) vec[idy[[i]]] <- to[i]
  return(vec)
}

#f02: savePDF
# save ggplot object to PDF
savePDF <- function(p, save.path, height = 4, width = 5.8){
  pdf(save.path,height = height, width = width)
  plot(p)
  dev.off()
}

#f03: staddata
# tranform data into [0-1], which 1 is depended by q
staddata <- function(x,q = 0.95){
  x <- x-min(x)
  qn <- quantile(x,q)
  if(qn==0) qn = max(x)
  x <- x/quantile(x,q)
  x[x>1] <- 1; return(x)
} 

#f04: Truncate data
# Truncate data based on q1 and q2
Trudata <- function(x, q1 = 0, q2 = 0.95){
  qx <- quantile(x, c(q1,q2))
  x[x<qx[1]] <- qx[1]
  x[x>qx[2]] <- qx[2]
  return(x)
} 

#f05: extexpr
extexpr <- function(seu.obj, gene, group = NULL, celltype = NULL, assay="RNA"){
  if(is.null(group)) group <- unique(seu.obj$group)
  if(is.null(celltype)) celltype <- unique(seu.obj$celltype)
  seu.obj@assays[[assay]]@data[gene, seu$group%in%group & seu$celltype%in%celltype]
}


#f06: which.two
which.two <- function(x){
  x[which.max(x)] <- min(x) 
  which.max(x)
} 



gene_density <- function(seu, gene, group = NULL, dims = 20, neighbors = 15){
  if(is.null(group)){
    seu[["pseudo_group"]] <- "nd"
    group <- "pseudo_group"
  }
  g <- unique(unlist(seu[[group]]))
  gene_dens <- lapply(g, function(x){
    seux <- seu[,seu[[group]]==x]
    gene_expr <- exp(seux@assays$RNA@data[gene,])-1
    pca <- Embeddings(seux, reduction = "pca")[,1:dims]
    knn <- FNN::get.knn(pca, k = neighbors)
    mtx <- sapply(1:neighbors,function(i) gene_expr[knn$nn.index[,i]])
    gene_dens <- (rowSums(mtx)+gene_expr)/(neighbors+1)
    names(gene_dens) <- colnames(seux)
    return(gene_dens)
  }) %>% unlist
  gene_dens[colnames(seu)]
}
positive_rate <- function(seu, gene, group = NULL, dims = 20, neighbors = 15){
  if(is.null(group)){
    seu[["pseudo_group"]] <- "nd"
    group <- "pseudo_group"
  }
  g <- unique(unlist(seu[[group]]))
  gene_pos <- lapply(g, function(x){
    seux <- seu[,seu[[group]]==x]
    gene_expr <- exp(seux@assays$RNA@data[gene,])-1
    gene_bi <- ifelse(gene_expr>0,1,0)
    pca <- Embeddings(seux, reduction = "pca")[,1:dims]
    knn <- FNN::get.knn(pca, k = neighbors)
    mtx <- sapply(1:neighbors,function(i) gene_bi[knn$nn.index[,i]])
    gene_pos <- (rowSums(mtx)+gene_bi)/(neighbors+1)
    names(gene_pos) <- colnames(seux)
    return(gene_pos)
  }) %>% unlist
  gene_pos[colnames(seu)]
}



findDEG_muscat <- function(seuobj, thres = 0.2, sample, group){
  require(muscat)
  DefaultAssay(seuobj) <- "RNA"
  seuobj <- seuobj[rowMeans(seuobj)>thres,]
  message(paste0("detect ", nrow(seuobj), " genes"))
  
  #SingleCellExperiment prepare
  sce <- as.SingleCellExperiment(seuobj)
  colData(sce)$sample_id <- colData(sce)[[sample]]
  colData(sce)$group_id <- colData(sce)[[group]]
  colData(sce)$cluster_id <- "1"
  
  #run muscat
  sce <- prepSCE(sce, kid = "cluster_id", gid = "group_id", sid = "sample_id")
  #differential state analysis
  ps <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id","sample_id"))
  res_DS <- muscat::pbDS(ps, method = "edgeR",filter = "none")
  
  #differential detection analysis
  pd <- muscat::aggregateData(sce, assay = "counts", fun = "num.detected", by = c("cluster_id","sample_id"))
  res_DD <- muscat::pbDS(pd, method = "DD", filter = "none")
  
  #output
  res <- muscat::stagewise_DS_DD(res_DS,res_DD,sce = sce)
  res <- data.frame(gene = res[[1]][[1]]$gene, p_adj = res[[1]][[1]]$p_adj,
                    DSp = res[[1]][[1]]$res_DS$p_val, DSpadj = res[[1]][[1]]$res_DS$p_adj.loc, DSfc = res[[1]][[1]]$res_DS$logFC,
                    DDp = res[[1]][[1]]$res_DD$p_val, DDpadj = res[[1]][[1]]$res_DD$p_adj.loc, DDfc = res[[1]][[1]]$res_DD$logFC,
                    row.names = res[[1]][[1]]$gene)
}


treat_lmer <- function(model){
  require(emmeans)
  emm <- emmeans(model, ~ sample_group | subcluster)
  treatment_effects <- pairs(emm, reverse = TRUE)  # T - NC
  contrast_effects <- pairs(treatment_effects, by = NULL, adjust = "fdr")
  effect_df <- list(
    effect_df = as.data.frame(treatment_effects) %>%
      mutate(`shifted_cells(CFD)(%)` = NA, `shifted_cells(HFD)(%)` = NA),
    contrast_effect = as.data.frame(contrast_effects) %>%
      dplyr::rename(fdr = `p.value`)
  )
  return(effect_df)
} 


shuffle_seu <- function(seuobj, n = 50){
  seuobj <- seuobj[rowSums(seuobj)>20,]
  mtx <- seuobj@assays$RNA@counts %>% as.matrix()
  seu_r <- pblapply(1:n, function(i){
    mtx_shuffle <- Rfast::rowShuffle(mtx)
    dimnames(mtx_shuffle) <- dimnames(mtx)
    seu_s <- CreateSeuratObject(mtx_shuffle,project = paste0("random",i))
    seu_s <- NormalizeData(seu_s)
  })
}



calc_shift <- function(seuobj, deg, orient = 1, thres = 0.99, group, control = "CFD"){
  seuobj@meta.data[,grep("Cluster",colnames(seuobj@meta.data))] <- NULL
  seuobj <- AddModuleScore(seuobj, list(deg))
  pca <- Embeddings(seuobj,reduction = "pca")[,1:15]
  knn <- FNN::get.knnx(data = pca[seuobj@meta.data[[group]]==control,], query = pca, k = 10)
  s <- seuobj@meta.data[["Cluster1"]]
  sr <- s[seuobj@meta.data[[group]]==control]
  sr_knn <- matrix(data = sr[knn$nn.index], ncol = ncol(knn$nn.index))
  s <- (s - rowMeans(sr_knn)) / apply(sr_knn,1,sd)
  if(orient==1){
    thres <- abs(qnorm(thres))
    shift_state <- ifelse(s > thres, "shifted", "stable")
  }else{
    thres <- -abs(qnorm(thres))
    shift_state <- ifelse(s < thres, "shifted", "stable")
  }
  return(shift_state)
}
# ----------- Part3 end -----------------



# ----------- Part4 miloR -----------------
#DOI: 10.1038/s41587-021-01033-z
#Once we have defined neighbourhoods, it’s good to take a look at how big the neighbourhoods are (i.e. how many cells form each neighbourhood). 
#This affects the power of DA testing. We can check this out using the plotNhoodSizeHist function. 
#Empirically, we found it’s best to have a distribution peaking between 50 and 100.
#Otherwise you might consider rerunning makeNhoods increasing k and/or prop
sel_kp <- function(seu.obj, k, d, p){
  require(miloR)
  require(SingleCellExperiment)
  scRNA_pre <- as.SingleCellExperiment(seu.obj)
  scRNA_milo <- Milo(scRNA_pre)
  pca_matrix <- Embeddings(seu.obj, "pca")  
  reducedDim(scRNA_milo, "PCA") <- pca_matrix
  scRNA_milo <- buildGraph(scRNA_milo, k = k, d = d, reduced.dim = "PCA")
  scRNA_milo <- makeNhoods(scRNA_milo, prop = p, k = k, d = d,
                           refined = T, reduced_dims = "PCA")
  message(ncol(scRNA_milo@nhoods))
  p <- plotNhoodSizeHist(scRNA_milo) + geom_vline(xintercept = c(0,50,100)); plot(p)
}


run_miloR <- function(seu.obj, k, d, p){
  require(miloR)
  require(SingleCellExperiment)
  #same to sel_kp
  scRNA_pre <- as.SingleCellExperiment(seu.obj)
  scRNA_milo <- Milo(scRNA_pre)
  pca_matrix <- Embeddings(seu.obj, "pca")  
  reducedDim(scRNA_milo, "PCA") <- pca_matrix
  scRNA_milo <- buildGraph(scRNA_milo, k = k, d = d, reduced.dim = "PCA")
  scRNA_milo <- makeNhoods(scRNA_milo, prop = p, k = k, d = d, refined = TRUE, reduced_dims = "PCA")
  scRNA_milo <- countCells(scRNA_milo, meta.data = as.data.frame(colData(scRNA_milo)), sample="sample")
  scRNA_milo <- calcNhoodDistance(scRNA_milo, d = d, reduced.dim = "PCA")
  
  #DA analysis
  traj_design <- data.frame(colData(scRNA_milo))[,c("sample", "sample_group")]
  traj_design <- distinct(traj_design)
  rownames(traj_design) <- traj_design$sample
  traj_design <- traj_design[colnames(nhoodCounts(scRNA_milo)), , drop=FALSE]
  traj_design$group <- factor(traj_design[["sample_group"]], level = c("CFD", "HFD"))
  da_results <- testNhoods(scRNA_milo,design = ~ group, design.df = traj_design)
  scRNA_milo <- buildNhoodGraph(scRNA_milo)
  rel <- list(scRNA_milo = scRNA_milo, da_results = da_results)
  return(rel)
}


plotnodegraph <- function (miloresults, seuobj, colour_by = "logFC", reddim = "UMAP", 
                           size_range = c(0.1, 4), node_stroke = 0.1){
  #prepare data
  require(igraph)
  require(ggraph)
  require(SingleCellExperiment)
  require(ggplot2)
  require(RColorBrewer)
  rangex <- function(x,m1,m2){
    x[x<m1] <- m1; x[x>m2] <- m2; return(x)
  }
  
  # 添加到 milo 对象
  umap_coords <- Embeddings(seuobj, "umap")
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  reducedDim(miloresults$scRNA_milo, "UMAP") <- umap_coords
  
  x <- miloresults$scRNA_milo
  signif_res <- miloresults$da_results
  if (colour_by %in% colnames(signif_res)){
    idy <- signif_res$Nhood[signif_res$PValue > 0.05]
    signif_res[signif_res$PValue > 0.05, colour_by] <- 0
    signif_res[[colour_by]] <- rangex(signif_res[[colour_by]], -2.5, 2.5)
    colData(x)[colour_by] <- NA
    colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]), colour_by] <- signif_res[, colour_by]
  }else {
    if(!colour_by %in% colnames(colData(x))) stop(colour_by, "is not a column in colData(x)")
  }
  
  #nh_graph construct
  nh_graph <- nhoodGraph(x)
  V(nh_graph)$colour_by <- colData(x)[as.numeric(vertex_attr(nh_graph)$name), colour_by]
  nh_graph <- permute(nh_graph, order(abs(vertex_attr(nh_graph)$colour_by),decreasing = T))
  layout <- reducedDim(x, reddim)[as.numeric(vertex_attr(nh_graph)$name), ]
  if (!any(class(layout) %in% c("matrix"))) {
    warning("Coercing layout to matrix format")
    layout <- as(layout, "matrix")
  }
  
  #plot data
  pl <- ggraph(simplify(nh_graph), layout = layout) + 
    geom_point_rast(data = data.frame(umap_coords), aes(UMAP1,UMAP2), size = 1, color = "#ededed")+
    geom_edge_link0(aes(width = weight), edge_colour = "#515151", edge_alpha = 0.01) +
    geom_node_point(aes(fill = colour_by,size = size), shape = 21, stroke = node_stroke, color = "#C1C1C1") + 
    scale_size(range = size_range, name = "Nhood size") + 
    scale_edge_width(range = c(0.1, 0.5), name = "overlap size") + 
    theme_void()
  if (is.numeric(V(nh_graph)$colour_by)) {
    pl <- pl + scale_fill_paletteer_c("ggthemes::Red-Black-White Diverging",direction = -1)
  }
  else {
    mycolors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(V(nh_graph)$colour_by)))
    pl <- pl + scale_color_manual(values = mycolors, name = colour_by, na.value = "white")
  }
  pl
}
#style2
plotnodegraph <- function(
    miloresults, seuobj,
    colour_by = "logFC",
    reddim = "UMAP",
    size_range = c(0.1, 4),
    node_stroke = 0.1,
    alpha_cutoff = 0.25,
    fc_limit = 2.5
) {
  require(igraph)
  require(ggraph)
  require(SingleCellExperiment)
  require(ggplot2)
  require(miloR)
  
  stopifnot(fc_limit > alpha_cutoff, alpha_cutoff >= 0)
  
  x <- miloresults$scRNA_milo
  da <- as.data.frame(miloresults$da_results)
  
  if (!all(c("Nhood", "logFC", "PValue") %in% colnames(da))) {
    stop("da_results 必须包含 Nhood、logFC 和 PValue")
  }
  
  # UMAP 按 milo 对象的细胞顺序对齐
  umap_coords <- SeuratObject::Embeddings(seuobj, "umap")
  
  if (!all(colnames(x) %in% rownames(umap_coords))) {
    stop("seuobj 的 UMAP 缺少 milo 对象中的部分细胞")
  }
  
  umap_coords <- umap_coords[colnames(x), 1:2, drop = FALSE]
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  reducedDim(x, "UMAP") <- umap_coords
  
  # 每个 neighborhood 的 index cell
  index_cells <- vapply(
    nhoodIndex(x),
    function(z) as.integer(z)[1L],
    integer(1)
  )
  
  nh_id <- as.integer(as.character(da$Nhood))
  
  if (anyNA(nh_id) ||
      any(nh_id < 1L | nh_id > length(index_cells))) {
    stop("da_results$Nhood 与 nhoodIndex(x) 不匹配")
  }
  
  # 将 DA 结果映射到 index cell
  cell_fc <- rep(NA_real_, ncol(x))
  cell_p <- rep(NA_real_, ncol(x))
  
  cell_fc[index_cells[nh_id]] <- da$logFC
  cell_p[index_cells[nh_id]] <- da$PValue
  
  if (colour_by %in% colnames(da)) {
    cell_colour <- da[[colour_by]][rep(NA_integer_, ncol(x))]
    cell_colour[index_cells[nh_id]] <- da[[colour_by]]
  } else {
    if (!colour_by %in% colnames(colData(x))) {
      stop(colour_by, " is not a column in colData(x)")
    }
    cell_colour <- colData(x)[[colour_by]]
  }
  
  # 图节点属性
  nh_graph <- nhoodGraph(x)
  node_cells <- as.integer(igraph::V(nh_graph)$name)
  
  node_fc <- cell_fc[node_cells]
  node_p <- cell_p[node_cells]
  node_colour <- cell_colour[node_cells]
  
  if (is.factor(node_colour)) {
    node_colour <- as.character(node_colour)
  }
  
  # logFC 颜色截断；不再将不显著结果设为 0
  if (colour_by == "logFC") {
    node_colour <- pmax(-fc_limit, pmin(fc_limit, node_colour))
  }
  
  # 低于阈值隐藏；达到阈值后由 0.15 渐变至 1
  node_alpha <- ifelse(
    is.na(node_fc) | abs(node_fc) < alpha_cutoff,
    0,
    0.15 + 0.85 *
      pmin((abs(node_fc) - alpha_cutoff) /
             (fc_limit - alpha_cutoff), 1)
  )
  
  igraph::V(nh_graph)$colour_by <- node_colour
  igraph::V(nh_graph)$node_alpha <- node_alpha
  igraph::V(nh_graph)$significant <- !is.na(node_p) & node_p < 0.05
  igraph::V(nh_graph)$abs_fc <- abs(node_fc)
  
  layout <- as.matrix(
    reducedDim(x, reddim)[node_cells, 1:2, drop = FALSE]
  )
  
  # 显式建立绘图数据，按显著性分层
  graph_layout <- ggraph::create_layout(
    igraph::simplify(nh_graph),
    layout = layout
  )
  
  visible_nodes <- as.data.frame(graph_layout)
  visible_nodes <- visible_nodes[
    visible_nodes$node_alpha > 0, , drop = FALSE
  ]
  
  # 同一层内，|logFC| 越大越晚绘制
  visible_nodes <- visible_nodes[
    order(visible_nodes$abs_fc), , drop = FALSE
  ]
  
  background_coords <- as.data.frame(
    reducedDim(x, reddim)[, 1:2, drop = FALSE]
  )
  colnames(background_coords) <- c("x", "y")
  
  pl <- ggraph::ggraph(graph_layout) +
    ggrastr::geom_point_rast(
      data = background_coords,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      size = 1,
      colour = "#ededed"
    ) +
    ggraph::geom_edge_link0(
      aes(width = weight),
      edge_colour = "#A1A1A1",
      edge_alpha = 0.01
    ) +
    
    # 第一层：不显著节点
    geom_point(
      data = visible_nodes[!visible_nodes$significant, , drop = FALSE],
      aes(x = x, y = y, fill = colour_by,
          size = size, alpha = node_alpha),
      inherit.aes = FALSE,
      shape = 21,
      stroke = node_stroke,
      colour = "#C1C1C1"
    ) +
    
    # 最上层：PValue < 0.05 的节点
    geom_point(
      data = visible_nodes[visible_nodes$significant, , drop = FALSE],
      aes(x = x, y = y, fill = colour_by,
          size = size, alpha = node_alpha),
      inherit.aes = FALSE,
      shape = 21,
      stroke = 0.5,
      colour = "black"
    ) +
    scale_alpha_identity() +
    scale_size(range = size_range, name = "Nhood size") +
    ggraph::scale_edge_width(
      range = c(0.1, 0.5),
      name = "overlap size"
    ) +
    theme_void()
  
  if (is.numeric(node_colour)) {
    pl <- pl +
      paletteer::scale_fill_paletteer_c(
        "ggthemes::Red-Black-White Diverging",
        direction = -1,
        name = colour_by
      )
  } else {
    categories <- unique(node_colour[!is.na(node_colour)])
    mycolors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(11, "Spectral")
    )(max(1L, length(categories)))
    
    pl <- pl +
      scale_fill_manual(
        values = mycolors,
        name = colour_by,
        na.value = "white"
      )
  }
  
  pl
}

# ----------- Part4 end -----------------





