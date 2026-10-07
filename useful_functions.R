
Geneplot <- function(seu_obj, gene, assay="RNA", reduction="umap", 
                     q = 0.95, size=0.1, reorder = T, 
                     raster = T, dpi = 300){
  require(ggsci)
  require(ggplot2)
  require(paletteer)
  require(ggrastr)
  umap <- data.frame(Embeddings(seu_obj,reduction = reduction))
  colnames(umap) <- c("UMAP_1","UMAP_2")
  if(size < 0.1) shape <- "." else shape <- 20
  p <- lapply(gene,function(x){
    message(paste(x))
    if(x %in% colnames(seu_obj)){
      umap[["expr"]] <- seu[[x]] %>% norm_scale(q2 = q)
    }else{
      umap[["expr"]] <- seu_obj@assays[[assay]]@data[x,] %>% norm_scale(q2 = q)
    }
    if(reorder) umap <- umap[order(umap$expr,decreasing = F),]
    if(raster){
      ggplot(umap, aes(UMAP_1,UMAP_2))+
        geom_point_rast(data = umap[umap$expr==0,], color = "#c0c0c0", size=size, shape = shape)+
        geom_point_rast(aes(color=expr,alpha = expr), size=size, shape = shape)+
        scale_color_gradientn(values = c(0,0.1,0.3,0.5,0.7,0.9,1),
                              colours = c( paletteer_c("grDevices::OrRd",30)[c(28,24,18,12,6,1)] ))+
        theme_minimal()+
        theme(legend.key.width = unit(0.5,"cm"),
              panel.background = element_rect(fill = "#FdFdFd", color = "#FdFdFd"))+
        labs(x=NULL,y=NULL,title = x)
    }else{
      ggplot(umap, aes(UMAP_1,UMAP_2))+
        geom_point(aes(color=expr,alpha = expr), size=size, shape = shape)+
        scale_color_gradientn(values = c(0,0.1,0.3,0.5,0.7,0.9,1),
                              colours = c( paletteer_c("grDevices::OrRd",30)[c(28,24,18,12,6,1)] ))+
        theme_minimal()+
        theme(legend.key.width = unit(0.5,"cm"),
              panel.background = element_rect(fill = "#FdFdFd", color = "#FdFdFd"))+
        labs(x=NULL,y=NULL,title = x)
    }
  })
  p
}

num_density <- function(umap,dist){
  knn <- FNN::get.knn(umap,k = 100)
  fx <- knn$nn.dist
  message(quantile(fx[,50],0.5))
  fx[fx<dist] <- 0
  fx[fx>dist] <- 1
  rowSums(1-fx)/100
}

reCreateSeuratObject <- function(seu.obj,dims=NULL){
  if("SCT"%in%Seurat::Assays(seu.obj)){
    mtx <- seu.obj@assays$SCT@counts
  }else{
    mtx <- seu.obj@assays$RNA@counts
  }
  mtx <- mtx[!grepl("^mt-", rownames(mtx)),]
  meta <- seu.obj@meta.data
  seu <- CreateSeuratObject(mtx)
  seu <- AddMetaData(seu,meta)
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  if(!is.null(dims)) seu <- RunUMAP(seu,dims=1:dims)
  return(seu)
}

create_cds <- function(seu.obj){
  require(monocle3)
  mtx <- seu.obj@assays$RNA@counts
  meta <- seu.obj@meta.data
  gene.meta <- data.frame(gene_short_name = rownames(mtx),features="none",row.names = rownames(mtx))
  cds <- new_cell_data_set(as.matrix(mtx),
                           cell_metadata = meta,
                           gene_metadata = gene.meta)
  return(cds)
}

Transfer_cds <- function(seu.obj,cds){
  cds <- preprocess_cds(cds, num_dim = 50, norm_method = "log")
  cds <- reduce_dimension(cds,max_components = 2,reduction_method="UMAP",preprocess_method = 'PCA')
  reducedDims(cds)[["PCA"]] <- Embeddings(seu.obj,reduction = "pca")
  reducedDims(cds)[["UMAP"]] <- Embeddings(seu.obj,reduction = "umap")
  cds <- cluster_cells(cds,cluster_method = 'louvain')
  s.clusters <- seu.obj$subcluster
  names(s.clusters) <- colnames(cds)
  s.clusters <- as.factor(s.clusters)
  cds@clusters$UMAP$clusters <- s.clusters
  return(cds)
}

learn_graph_cds <- function(cds, branch = 25){
  graph_control <- setNames(list(1, 0.5, branch,
                                 FALSE, TRUE, FALSE, NULL, 
                                 10, 1e-5, 0.05, 0.01), 
                            c("euclidean_distance_ratio", "geodesic_distance_ratio", "minimal_branch_len",
                              "orthogonal_proj_tip", "prune_graph", "scale", "rann.k",
                              "maxiter", "eps", "L1.gamma", "L1.sigma"))
  cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = graph_control)
}

norm_scale <- function(x, q1=0.01, q2=0.99, limq2=0.01){
  y <- x[!is.na(x)]
  if(max(y)==1 & min(y)==0){return(x)}
  y <- y - quantile(y,q1)
  y[y<0] <- 0
  q2 <- quantile(y,q2)
  if(q2 < limq2) y <- replicate(length(y),0) else y <- y/q2
  y[y>1] <- 1 
  x[!is.na(x)] <- y
  return(x)
}

floor2 <- function(x){
  y <- x[!is.na(x)]
  if(max(y)*10 %% 1 == 0) y[y==max(y)] <- max(y)-0.0001
  y <- floor(y*10)/10
  x[!is.na(x)] <- y
  return(x)
}

savePDF <- function(p,file,width,height){
  pdf(file, width = width, height = height)
  plot(p)
  dev.off()
}

binary_matrix <- function(matrix){
  matrix[matrix>1] <- 1; return(matrix)
}

compare_curves <- function(x, y1, y2) {
  calculate_derivatives <- function(x, y) {
    dy <- diff(y)
    dx <- diff(x)
    derivatives <- dy / dx
    return(derivatives)
  }
  deriv1 <- calculate_derivatives(x, y1)
  deriv2 <- calculate_derivatives(x, y2)
  deriv_cor <- cor(y1, y2, use = "complete.obs")
  #DTW distance
  library(dtw)
  dtw_deriv <- dtw(deriv1, deriv2)$distance
  results <- list(derivative_correlation = deriv_cor,
                  dtw_distance = dtw_deriv)
}




