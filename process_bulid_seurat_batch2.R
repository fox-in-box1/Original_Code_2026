#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(glmGamPoi)
  library(DoubletFinder)
  library(qs)
})

options(future.globals.maxSize = 200 * 1024^3)
set.seed(20260907)

base_dir <- "/home/disk1/wgs/batch2"
output_file <- file.path(base_dir, "data.qs")
temporary_file <- paste0(output_file, ".tmp")

sample_info <- data.frame(
  sample = c("HX327", "HX328", "HX329", "HX330"),
  sample_group = c("CFD", "HFD", "CFD", "HFD"),
  h5_file = file.path(base_dir, c(
    "AM3HX327329_CFD/AM3HX327/filtered_feature_bc_matrix.h5",
    "AM3HX328330_HFD/AM3HX328/filtered_feature_bc_matrix.h5",
    "AM3HX327329_CFD/AM3HX329/filtered_feature_bc_matrix.h5",
    "AM3HX328330_HFD/AM3HX330/filtered_feature_bc_matrix.h5"
  )),
  stringsAsFactors = FALSE
)

missing_files <- sample_info$h5_file[!file.exists(sample_info$h5_file)]
if (length(missing_files)) {
  stop("Missing count matrices: ", paste(missing_files, collapse = ", "))
}

read_gene_expression <- function(path) {
  counts <- Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("No Gene Expression matrix in ", path)
    }
    counts <- counts[["Gene Expression"]]
  }
  counts
}

process_sample <- function(sample_name, sample_group, h5_file) {
  message("[", sample_name, "] reading raw counts")
  counts <- read_gene_expression(h5_file)
  colnames(counts) <- paste(sample_name, colnames(counts), sep = "_")
  obj <- CreateSeuratObject(counts = counts, project = sample_name)
  obj$sample <- sample_name
  obj$sample_group <- sample_group
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^(p|e|s)]")

  before <- ncol(obj)
  obj <- subset(obj, subset = nFeature_RNA < 7500 & percent.mt < 10)
  message("[", sample_name, "] QC retained ", ncol(obj), "/", before, " cells")

  obj <- SCTransform(
    obj, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE
  )
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)

  message("[", sample_name, "] DoubletFinder parameter sweep")
  sweep <- paramSweep_v3(obj, PCs = 1:30, sct = TRUE, num.cores = 8)
  bcmvn <- find.pK(summarizeSweep(sweep, GT = FALSE))
  valid <- is.finite(bcmvn$BCmetric) & !is.na(bcmvn$pK)
  if (!any(valid)) stop("No valid DoubletFinder pK for ", sample_name)
  best <- which(valid)[which.max(bcmvn$BCmetric[valid])]
  best_pk <- as.numeric(as.character(bcmvn$pK[best]))
  homotypic_prop <- modelHomotypic(obj$seurat_clusters)
  n_exp_raw <- round(0.05 * ncol(obj))
  n_exp <- round(n_exp_raw * (1 - homotypic_prop))

  obj <- doubletFinder_v3(
    obj, PCs = 1:30, pN = 0.25, pK = best_pk, nExp = n_exp,
    reuse.pANN = FALSE, sct = TRUE
  )
  df_column <- tail(grep("^DF.classifications", colnames(obj@meta.data), value = TRUE), 1)
  if (length(df_column) != 1L) stop("Missing DoubletFinder classification for ", sample_name)
  obj$doubletfinder <- as.character(obj@meta.data[[df_column]])
  obj$doubletfinder_pK <- best_pk
  obj$doubletfinder_nExp_raw <- n_exp_raw
  obj$doubletfinder_nExp <- n_exp
  obj$doubletfinder_homotypic_prop <- homotypic_prop
  obj@meta.data <- obj@meta.data[, !grepl("^(pANN|DF.classifications)", colnames(obj@meta.data)), drop = FALSE]
  message("[", sample_name, "] pK=", best_pk, ", expected doublets=", n_exp)
  obj
}

objects <- lapply(seq_len(nrow(sample_info)), function(i) {
  process_sample(sample_info$sample[i], sample_info$sample_group[i], sample_info$h5_file[i])
})
names(objects) <- sample_info$sample

data <- merge(objects[[1]], y = objects[-1], project = "batch2_HX327_HX330")
data$sample <- factor(as.character(data$sample), levels = sample_info$sample)
data$orig.ident <- factor(as.character(data$orig.ident), levels = sample_info$sample)
data$sample_group <- factor(as.character(data$sample_group), levels = c("CFD", "HFD"))

stopifnot(
  !anyDuplicated(colnames(data)),
  !anyNA(data$doubletfinder),
  all(data$nFeature_RNA < 7500),
  all(data$percent.mt < 10),
  all(as.character(data$sample_group) == ifelse(
    as.character(data$sample) %in% c("HX327", "HX329"), "CFD", "HFD"
  ))
)

qsave(data, temporary_file, preset = "high")
check <- qread(temporary_file)
stopifnot(inherits(check, "Seurat"), identical(dim(check), dim(data)))
rm(check)
if (file.exists(output_file) && !file.remove(output_file)) stop("Cannot replace ", output_file)
if (!file.rename(temporary_file, output_file)) stop("Cannot rename temporary output")

message("Saved and validated: ", output_file)
print(table(sample = data$sample, group = data$sample_group))
print(table(sample = data$sample, doubletfinder = data$doubletfinder))
