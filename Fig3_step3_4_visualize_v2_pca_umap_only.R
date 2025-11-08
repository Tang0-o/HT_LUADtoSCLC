#!/usr/bin/env Rscript
# ----------------------------------------------------------
# 仅保留 PCA / UMAP 可视化
# ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# 固定 Lancet 配色
LANCET_DATASET <- c(
  "NT_bulk"= "#FF8C00",  # DarkOrange / 深橙色
  "TCGA"   = "#3F51B5",
  "OAK"    = "#ED0000",
  "POPLAR" = "#42B540",
  "ONCOSG" = "#0099B4"
)

LANCET_CLUSTER <- c(
  "cluster_high" = "#FDAF91",
  "cluster_low"  = "#00468B"
)

# ----------------------------------------------------------
# 1. 运行 PCA + UMAP
# ----------------------------------------------------------
run_pca_umap <- function(seurat_obj,
                         n_features = 2000,
                         npcs       = 50) {
  set.seed(2025)
  
  # 选高变基因
  gene_vars <- apply(seurat_obj[["RNA"]]@data, 1, var)
  var_feats <- head(names(sort(gene_vars, decreasing = TRUE)), n_features)
  VariableFeatures(seurat_obj) <- var_feats
  
  # 必须先做 ScaleData
  seurat_obj <- ScaleData(seurat_obj, features = var_feats, verbose = FALSE)
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs, features = var_feats, verbose = FALSE)
  
  # UMAP
  pca_emb <- Embeddings(seurat_obj, reduction = "pca")[, 1:npcs]
  umap_xy <- uwot::umap(pca_emb,
                        n_neighbors = 30,
                        min_dist    = 0.3,
                        metric      = "correlation",
                        n_epochs    = 500,
                        verbose     = FALSE)
  colnames(umap_xy) <- c("UMAP_1", "UMAP_2")
  seurat_obj[["umap"]] <- CreateDimReducObject(umap_xy, key = "UMAP_")
  
  return(seurat_obj)
}

# ----------------------------------------------------------
# 2. 绘制 PCA / UMAP 图
# ----------------------------------------------------------
plot_pca_umap <- function(seurat_obj,
                          out_dir = "pca_umap_plots") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 仅保留 cluster_high / cluster_low
  keep <- which(seurat_obj$cluster %in% c("cluster_high", "cluster_low"))
  if (length(keep) == 0) stop("No cluster_high or cluster_low samples found.")
  
  obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep])
  
  highlight_cells <- colnames(obj)[obj$orig.ident == "NT_bulk"]
  
  for (reduction in c("pca", "umap")) {
    if (!reduction %in% names(obj@reductions)) next
    
    # 数据集着色
    p1 <- DimPlot(obj, reduction = reduction,
                  group.by = "orig.ident",
                  pt.size  = 0.8,
                  cols     = LANCET_DATASET) +
      ggtitle(paste0(toupper(reduction), " by Dataset")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text  = element_text(size = 10),
            axis.title = element_text(face = "bold", size = 12),
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    
    # 聚类着色
    p2 <-  CellDimPlot(obj, reduction = reduction,group.by = "cluster",
                       pt.size =0.8,label.size=12, lineages_line_bg_stroke = 0,streamline_bg_stroke = 0,
                       cells.highlight = highlight_cells,  cols.highlight = "#F05C3BFF",  sizes.highlight = 1,  alpha.highlight = 1,  stroke.highlight = 1,
                       label = F , label_insitu =T)+
      ggtitle(paste0(toupper(reduction), " by Cluster")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text  = element_text(size = 10),
            axis.title = element_text(face = "bold", size = 12),
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    
    # 合并并保存
    combined <- p1 + p2 + plot_layout(ncol = 2)
    pdf(file.path(out_dir, paste0(reduction, "_visualization.pdf")),
        width = 8, height = 3)
    print(combined)
    dev.off()
  }
}

# ----------------------------------------------------------
# 3. 主函数
# ----------------------------------------------------------
main <- function(seurat_obj, out_dir = "pca_umap_plots") {
  if (!"cluster" %in% colnames(seurat_obj@meta.data))
    stop("Seurat object must contain a 'cluster' column.")
  
  message("Running PCA & UMAP ...")
  seurat_obj <- run_pca_umap(seurat_obj)
  
  message("Plotting PCA & UMAP ...")
  plot_pca_umap(seurat_obj, out_dir)
  
  message("Done! Plots saved in: ", normalizePath(out_dir))
  invisible(seurat_obj)
}

# 示例
# main(sc_4bulk, "my_output_dir")