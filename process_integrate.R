#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(glmGamPoi)
  library(qs)
})

options(future.globals.maxSize = 300 * 1024^3)
set.seed(20260907)

base_dir <- "/home/disk1/wgs/batch2"
batch2_file <- file.path(base_dir, "data.qs")
mapped_file <- "/home/disk1/wgs/NEW/s03_Integration Data mapped k5.qs"
output_file <- file.path(base_dir, "data_SCT_integrated.qs")
temporary_file <- paste0(output_file, ".tmp")

sample_order <- c(
  "AM3N378", "AM3N379", "AM3N382", "AM3N383",
  "HX327", "HX328", "HX329", "HX330"
)
reference_samples <- c("AM3N378", "AM3N382", "HX327", "HX328")
cfd_samples <- c("AM3N378", "AM3N379", "HX327", "HX329")

batch2 <- qread(batch2_file)
mapped <- qread(mapped_file)
stopifnot(inherits(batch2, "Seurat"), inherits(mapped, "Seurat"))
if (!setequal(unique(as.character(batch2$sample)), sample_order[5:8])) {
  stop("Unexpected batch2 samples")
}
if (!setequal(unique(as.character(mapped$orig.ident)), sample_order[1:4])) {
  stop("Unexpected mapped samples")
}
if (!identical(rownames(GetAssayData(batch2, assay = "RNA", slot = "counts")),
               rownames(GetAssayData(mapped, assay = "RNA", slot = "counts")))) {
  stop("RNA count feature names/order differ between inputs")
}

target_columns <- colnames(mapped@meta.data)
make_object <- function(source, sample_name, is_batch2) {
  sample_field <- if (is_batch2) "sample" else "orig.ident"
  cells <- rownames(source@meta.data)[as.character(source@meta.data[[sample_field]]) == sample_name]
  counts <- GetAssayData(source, assay = "RNA", slot = "counts")[, cells, drop = FALSE]
  obj <- CreateSeuratObject(counts = counts, project = sample_name)
  obj$percent.mt <- as.numeric(source@meta.data[cells, "percent.mt"])
  obj$sample_group <- if (sample_name %in% cfd_samples) "CFD" else "HFD"
  obj$doubletfinder <- as.character(source@meta.data[cells, "doubletfinder"])
  if (is_batch2) {
    obj$celltype <- NA_character_
    obj$cluster <- NA_character_
    obj$cluster_id <- NA_character_
  } else {
    obj$celltype <- as.character(source@meta.data[cells, "celltype"])
    obj$cluster <- as.character(source@meta.data[cells, "cluster"])
    obj$cluster_id <- as.character(source@meta.data[cells, "cluster_id"])
  }
  obj
}

objects <- setNames(vector("list", length(sample_order)), sample_order)
for (sample_name in sample_order[1:4]) objects[[sample_name]] <- make_object(mapped, sample_name, FALSE)
for (sample_name in sample_order[5:8]) objects[[sample_name]] <- make_object(batch2, sample_name, TRUE)
rm(batch2, mapped)
invisible(gc())

for (i in seq_along(objects)) {
  message("SCTransform: ", names(objects)[i], " (", ncol(objects[[i]]), " cells)")
  objects[[i]] <- SCTransform(
    objects[[i]], vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE
  )
}

features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 3000)
objects <- PrepSCTIntegration(object.list = objects, anchor.features = features, verbose = FALSE)
objects <- lapply(objects, RunPCA, features = features, npcs = 50, verbose = FALSE)

reference_indices <- match(reference_samples, names(objects))
anchors <- FindIntegrationAnchors(
  object.list = objects,
  normalization.method = "SCT",
  anchor.features = features,
  reference = reference_indices,
  reduction = "rpca",
  k.anchor = 5,
  dims = 1:30,
  verbose = TRUE
)
data <- IntegrateData(
  anchorset = anchors, normalization.method = "SCT", dims = 1:30, verbose = TRUE
)
rm(anchors)
invisible(gc())

DefaultAssay(data) <- "integrated"
data <- RunPCA(data, npcs = 50, verbose = FALSE)
data <- RunUMAP(data, dims = 1:15, spread = 2, verbose = FALSE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:15, verbose = FALSE)
data <- FindClusters(data, resolution = 1, verbose = FALSE)

missing_columns <- setdiff(target_columns, colnames(data@meta.data))
if (length(missing_columns)) stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
data@meta.data <- data@meta.data[, target_columns, drop = FALSE]
data$sample_group <- factor(as.character(data$sample_group), levels = c("CFD", "HFD"))

batch2_cells <- data$orig.ident %in% sample_order[5:8]
stopifnot(
  setequal(unique(as.character(data$orig.ident)), sample_order),
  !anyDuplicated(colnames(data)),
  !anyNA(data$sample_group),
  !anyNA(data$doubletfinder),
  all(as.character(data$sample_group) == ifelse(data$orig.ident %in% cfd_samples, "CFD", "HFD")),
  all(is.na(data$celltype[batch2_cells])),
  all(is.na(data$cluster[batch2_cells])),
  all(is.na(data$cluster_id[batch2_cells])),
  all(c("RNA", "SCT", "integrated") %in% Assays(data)),
  all(c("pca", "umap") %in% Reductions(data))
)

qsave(data, temporary_file, preset = "high")
check <- qread(temporary_file)
stopifnot(inherits(check, "Seurat"), identical(dim(check), dim(data)))
rm(check)
if (file.exists(output_file) && !file.remove(output_file)) stop("Cannot replace ", output_file)
if (!file.rename(temporary_file, output_file)) stop("Cannot rename temporary output")

message("Saved and validated: ", output_file)
message("Dimensions: ", nrow(data), " features x ", ncol(data), " cells")
print(table(sample = data$orig.ident, group = data$sample_group))
