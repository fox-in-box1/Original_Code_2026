library(circlize)
library(dplyr)

# ============================================
# 1. 数据准备（你的原始代码）
# ============================================

ligand <- df1 %>% group_by(ligand) %>% summarise(size = sum(value), cell = unique(source))
receptor <- df1 %>% group_by(receptor) %>% summarise(size = sum(value), cell = unique(target))
colnames(ligand) <- colnames(receptor) <- c("nodes","size", "cell")
nodes <- rbind(ligand,receptor) %>% mutate(cell = factor(cell, levels = names(cell_col))) %>% arrange(cell,desc(size))

# ============================================
# 2. 节点属性（统一字体大小，不区分）
# ============================================

node_label_text <- sub("_[^_]*$", "", nodes$nodes)
node_col <- cell_col[match(nodes$cell, names(cell_col))]

# 统一大小和颜色（关键修复：用 rep 创建与节点数相同的向量）
n_nodes <- nrow(nodes)
node_label_size <- rep(1.2, n_nodes)  # 统一大小
node_label_col <- rep("black", n_nodes)  # 统一颜色

names(node_label_text) <- names(node_label_size) <- names(node_label_col) <- names(node_col) <- nodes$nodes

# ============================================
# 3. 去重节点
# ============================================

unique_nodes <- unique(nodes$nodes)
n_sectors <- length(unique_nodes)

nodes_unique <- nodes %>%
  group_by(nodes) %>%
  summarise(size = sum(size), cell = first(cell), .groups = "drop") %>%
  mutate(cell = factor(cell, levels = names(cell_col))) %>%
  arrange(cell, desc(size))

node_col_unique <- cell_col[match(nodes_unique$cell, names(cell_col))]
names(node_col_unique) <- nodes_unique$nodes

node_label_text_unique <- sub("_[^_]*$", "", nodes_unique$nodes)
node_label_size_unique <- rep(1.2, n_sectors)  # 统一
node_label_col_unique <- rep("black", n_sectors)  # 统一

names(node_label_text_unique) <- names(node_label_size_unique) <- names(node_label_col_unique) <- names(node_col_unique) <- nodes_unique$nodes

# ============================================
# 4. 计算 gap.after（减小组间间隙）
# ============================================

cell_counts <- table(nodes_unique$cell)
cell_names <- names(cell_col)

gap_after <- c()
for (i in seq_along(cell_names)) {
  cell_name <- cell_names[i]
  count <- as.numeric(cell_counts[cell_name])
  if (!is.na(count) && count > 0) {
    if (count > 1) gap_after <- c(gap_after, rep(1.5, count - 1))
    if (i < length(cell_names)) gap_after <- c(gap_after, 6)
  }
}

if (length(gap_after) < n_sectors) {
  gap_after <- c(gap_after, rep(1.5, n_sectors - length(gap_after)))
} else if (length(gap_after) > n_sectors) {
  gap_after <- gap_after[1:n_sectors]
}

# ============================================
# 5. 绘制
# ============================================

pdf("NEW/cellchat/p06_chord.pdf", width = 10, height = 10)
par(mar = c(0, 0, 0, 0))
circos.clear()

circos.par(
  start.degree = 90,
  gap.after = gap_after,
  track.margin = c(0, 0.02),
  points.overflow.warning = FALSE
)

# 绘制弦图
chordDiagramFromDataFrame(
  df = df1,
  grid.col = node_col_unique,
  col = df1$linecol,
  transparency = 0.2,
  annotationTrack = "grid",
  annotationTrackHeight = mm_h(2),
  
  preAllocateTracks = list(
    list(track.height = mm_h(18)),  # track 1: 标签
    list(track.height = mm_h(4))    # track 2: 色块
  ),
  
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  link.arr.length = 0.05,
  link.arr.width = 0.08,
  
  order = nodes_unique$nodes,
  link.largest.ontop = TRUE
)

# ============================================
# 6. 标签轨道（track 1，最内侧靠近弦）
# ============================================

circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    xcenter <- get.cell.meta.data("xcenter")
    ylim <- get.cell.meta.data("ylim")
    
    ypos <- ylim[1] + 0.5
    
    circos.text(
      x = xcenter,
      y = ypos,
      labels = node_label_text_unique[sector.index],
      facing = "clockwise",
      niceFacing = TRUE,
      cex = node_label_size_unique[sector.index],
      col = node_label_col_unique[sector.index],
      font = 1,
      adj = c(0.5, 1)
    )
  }
)

# ============================================
# 7. 色块轨道（track 2，标签外侧）
# ============================================

circos.trackPlotRegion(
  track.index = 2,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    
    col <- node_col_unique[sector.index]
    
    circos.rect(
      xleft = get.cell.meta.data("xlim")[1],
      xright = get.cell.meta.data("xlim")[2],
      ybottom = get.cell.meta.data("ylim")[1],
      ytop = get.cell.meta.data("ylim")[2],
      col = col,
      border = NA
    )
  }
)

dev.off()
circos.clear()