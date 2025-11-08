#!/usr/bin/env Rscript

###
### 脚本2: 为BayesPrism解卷积结果生成可视化和报告
###

# 设置工作目录
output_dir <- "/home/data/tmh_project/SCLC/Fig4/1_HTAN_downsample_bayesprism/"
setwd(output_dir)

# 加载所需包 - 确保正确加载顺序
library(tidyverse)  # 这会加载dplyr、ggplot2等
# 如果tidyverse加载有问题，可以单独加载
if(!requireNamespace("dplyr", quietly = TRUE)) {
  library(dplyr)
}
if(!requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
}
library(ggpubr)  # 用于stat_compare_means函数
library(pheatmap) # 用于绘制热图
library(RColorBrewer) # 用于热图配色
library(patchwork)
library(scales)
library(viridis)  # 用于高质量配色方案

# 系统报错改为英文
Sys.setenv(LANGUAGE = "en")
# 禁止转化为因子
options(stringsAsFactors = FALSE)

# 创建可视化输出目录
vis_dir <- '../2_bayesprism_visualization'
dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)

# 创建各类图形的子目录
heatmap_dir <- file.path(vis_dir, "heatmaps")
line_dir <- file.path(vis_dir, "line_plots")
bar_dir <- file.path(vis_dir, "bar_plots")
area_dir <- file.path(vis_dir, "area_plots")
box_dir <- file.path(vis_dir, "box_plots")
comparison_dir <- file.path(vis_dir, "comparisons")

# 创建所有子目录
for (dir in c(heatmap_dir, line_dir, bar_dir, area_dir, box_dir, comparison_dir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# 读取脚本1生成的结果
cat("读取解卷积分析结果...\n")
results_file <- file.path(output_dir, "all_results.rds")

if (!file.exists(results_file)) {
  stop("无法找到解卷积分析结果文件。请先运行process_data.R脚本。")
}

results <- readRDS(results_file)
cat("成功读取结果数据\n")

# 定义分组和细胞类型顺序
group_orders <- list(
  NT = c("C_LUAD", "T_LUAD", "T_SCLC", "C_SCLC"),
  TCGA = c("cluster_low", "cluster_medium", "cluster_high"),
  OAK_POPLAR = c("cluster_low", "cluster_medium", "cluster_high")
)

# 更新细胞类型顺序
cell_type_order <- c("NSCLC", "SCLC-A", "SCLC-N", "SCLC-P", 
                     "B cell", "T cell", "DC", "Macrophage", "Mast", 
                     "Neutrophil", "Plasma cell", "Endothelial", "Fibroblast")

#cancer_epi <- c("NSCLC", "SCLC-A", "SCLC-N", "SCLC-P")
cancer_epi <- c( "SCLC-A", "SCLC-N", "SCLC-P")
type1 <- c( "SCLC-A", "SCLC-P")
type2 <- c(  "SCLC-N")
immu <- c("B cell", "T cell", "DC", "Macrophage", "Mast", "Neutrophil", "Plasma cell")

# CNS风格配色 - 使用提供的调色板
colsp <- c('#46732EFF','#9370DB','#FD7446FF','#197EC0FF','#FED439FF',
           '#8A9197FF','#D2AF81FF','#709AE1FF','#FF91AF','#F05C3BFF',
           '#D5E4A2FF','#F8D568','#71D0F5FF','#370335FF','#075149FF',
           '#C80813FF','#91331FFF','#1A9993FF','#FD8CC1FF','#FF6700',
           '#00AD43','#89CFF0','#BA160C','#A6A6A6','#006DB0',
           '#C154C1','#D99A6C','#96C8A2','#FBEC5D')

# 为细胞类型分配颜色 - 使用自定义颜色
main_colors <- c(
  "NSCLC" = "#ADB6CA",
  "SCLC-A" = "#ED0000",
  "SCLC-P" = "#FDAF91",
  "SCLC-N" = "#00468B",
  "B cell" = colsp[2],
  "T cell" = colsp[8],
  "DC" = colsp[9],
  "Macrophage" = colsp[10],
  "Mast" = colsp[11],
  "Neutrophil" = colsp[12],
  "Plasma cell" = colsp[13],
  "Endothelial" = colsp[18],
  "Fibroblast" = colsp[19]
)

# 组别颜色
group_colors <- list(
  NT = c(
    "C_LUAD" = "#ADB6CA",
    "T_LUAD" = "#FDAF91", 
    "C_SCLC" = "#00468B",
    "T_SCLC" = "#ED0000"
  ),
  TCGA = c(
    "cluster_low" = "#00468B",
    "cluster_medium" = "#ADB6CA",
    "cluster_high" = "#FDAF91"
  ),
  OAK_POPLAR = c(
    "cluster_low" = "#00468B",
    "cluster_medium" = "#ADB6CA",
    "cluster_high" = "#FDAF91"
  )
)

# 数据集颜色 - 使用更鲜明的对比色
dataset_colors <- c(
  "NT" = "#1B9E77",       # 绿松石色
  "TCGA" = "#D95F02",     # 橙色
  "OAK_POPLAR" = "#7570B3" # 紫色
)

# 准备数据函数
prepare_data <- function() {
  all_data <- list()
  
  for (dataset_name in names(results)) {
    result <- results[[dataset_name]]
    
    if(!is.null(result$theta) && !is.null(result$meta_df)) {
      # 合并细胞比例和元数据
      theta_df <- as.data.frame(result$theta)
      # 使用dplyr::前缀明确指定使用dplyr的函数
      data_df <- theta_df %>%
        dplyr::mutate(Sample = rownames(theta_df)) %>%
        dplyr::left_join(result$meta_df, by = "Sample")
      
      # 确保分组顺序正确
      if(all(group_orders[[dataset_name]] %in% data_df$group)) {
        data_df$group <- factor(data_df$group, levels = group_orders[[dataset_name]])
      } else {
        data_df$group <- factor(data_df$group)
      }
      
      # 添加数据集标识
      data_df$Dataset <- dataset_name
      
      all_data[[dataset_name]] <- data_df
      cat("处理了", dataset_name, "数据集，包含", nrow(data_df), "个样本\n")
    } else {
      cat("警告：", dataset_name, "数据集缺少必要的结果数据\n")
    }
  }
  
  return(all_data)
}

# 1. 绘制热图函数 - 修改以适应大数据集和新的细胞类型顺序，并增强视觉效果
plot_heatmap <- function(data_df, output_path, title, cluster_rows = TRUE, cluster_cols = TRUE, 
                         annotation_colors = NULL, show_rownames = FALSE) {
  # 提取细胞比例数据
  cell_cols <- setdiff(colnames(data_df), c("Sample", "group", "Dataset"))
  
  # 检查数据是否足够
  if(length(cell_cols) == 0 || nrow(data_df) == 0) {
    cat("警告：没有足够的数据绘制热图\n")
    return(NULL)
  }
  
  # 准备热图数据 - 使用dplyr::前缀
  heatmap_data <- data_df %>%
    dplyr::select(dplyr::all_of(c("Sample", "group", cell_cols))) %>%
    dplyr::arrange(group)
  
  # 重新排序细胞类型列
  available_cell_types <- intersect(cell_type_order, cell_cols)
  other_cell_types <- setdiff(cell_cols, cell_type_order)
  ordered_cell_cols <- c(available_cell_types, other_cell_types)
  
  # 提取细胞比例矩阵
  cell_prop_matrix <- as.matrix(heatmap_data[, ordered_cell_cols])
  rownames(cell_prop_matrix) <- heatmap_data$Sample
  
  # 创建注释数据
  annotation_df <- data.frame(
    Group = heatmap_data$group,
    row.names = heatmap_data$Sample
  )
  
  # 如果没有提供注释颜色，则创建
  if(is.null(annotation_colors)) {
    # 获取组别
    groups <- unique(annotation_df$Group)
    dataset_name <- data_df$Dataset[1]
    
    # 使用预定义的组别颜色
    if(dataset_name %in% names(group_colors)) {
      group_color_palette <- group_colors[[dataset_name]]
      # 确保所有组都有颜色
      if(!all(groups %in% names(group_color_palette))) {
        missing_groups <- setdiff(groups, names(group_color_palette))
        additional_colors <- viridis(length(missing_groups))
        names(additional_colors) <- missing_groups
        group_color_palette <- c(group_color_palette, additional_colors)
      }
    } else {
      # 如果没有预定义颜色，使用viridis调色板
      group_color_palette <- viridis(length(groups))
      names(group_color_palette) <- groups
    }
    
    annotation_colors <- list(Group = group_color_palette)
  }
  
  # 计算适当的单元格高度 - 根据样本数量动态调整
  sample_count <- nrow(cell_prop_matrix)
  cellheight <- ifelse(sample_count > 100, 0.5, 
                       ifelse(sample_count > 50, 1, 
                              ifelse(sample_count > 20, 3, 8)))
  
  # 使用更高质量的热图配色方案
  heatmap_colors <- colorRampPalette(c(
    "#053061", "#2166ac", "#4393c3", "#92c5de", 
    "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", 
    "#d6604d", "#b2182b", "#67001f"
  ))(100)
  
  # 绘制热图
  pheatmap(
    cell_prop_matrix,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    main = title,
    color = heatmap_colors,
    scale = "column",  # 按列进行标准化
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,  # 不显示样本名称
    show_colnames = TRUE,
    fontsize = 14,  # 增加字体大小
    fontsize_row = 10,
    fontsize_col = 14,
    border_color = NA,
    cellwidth = 18,  # 增加单元格宽度
    cellheight = cellheight,  # 动态调整单元格高度
    angle_col = 90,  # 倾斜列名以防止重叠
    treeheight_row = ifelse(cluster_rows, 30, 0),  # 调整树状图高度
    treeheight_col = ifelse(cluster_cols, 30, 0)
  )
}

# 2. 绘制每个数据集的热图
plot_dataset_heatmaps <- function(all_data) {
  cat("绘制各数据集细胞比例热图...\n")
  
  for (dataset_name in names(all_data)) {
    data_df <- all_data[[dataset_name]]
    
    # 带聚类的热图
    pdf_file <- file.path(heatmap_dir, paste0("heatmap_", dataset_name, "_clustered.pdf"))
    pdf(pdf_file, width = 10, height = 8)
    plot_heatmap(
      data_df, 
      pdf_file,
      title = paste(dataset_name, "- Cell Type Proportions (Clustered)"),
      cluster_rows = TRUE, 
      cluster_cols = TRUE,
      show_rownames = FALSE  # 不显示样本名称
    )
    dev.off()
    
    # 不聚类的热图
    pdf_file <- file.path(heatmap_dir, paste0("heatmap_", dataset_name, "_ordered.pdf"))
    pdf(pdf_file, width = 10, height = 8)
    plot_heatmap(
      data_df, 
      pdf_file,
      title = paste(dataset_name, "- Cell Type Proportions (Ordered)"),
      cluster_rows = FALSE, 
      cluster_cols = FALSE,
      show_rownames = FALSE  # 不显示样本名称
    )
    dev.off()
    
    # 同时保存PNG版本
    png_file <- file.path(heatmap_dir, paste0("heatmap_", dataset_name, "_ordered.png"))
    png(png_file, width = 10, height = 8, units = "in", res = 300)
    plot_heatmap(
      data_df, 
      png_file,
      title = paste(dataset_name, "- Cell Type Proportions (Ordered)"),
      cluster_rows = FALSE, 
      cluster_cols = FALSE,
      show_rownames = FALSE  # 不显示样本名称
    )
    dev.off()
  }
}

# 3. 绘制TCGA和OAK_POPLAR整合热图 - 修改为按Group分组展示不同dataset的celltypes分布
plot_integrated_heatmap <- function(all_data) {
  cat("绘制TCGA和OAK_POPLAR整合热图...\n")
  
  # 检查是否有这两个数据集
  if(!all(c("TCGA", "OAK_POPLAR") %in% names(all_data))) {
    cat("警告：缺少TCGA或OAK_POPLAR数据集，无法绘制整合热图\n")
    return()
  }
  
  # 合并数据
  integrated_data <- rbind(all_data[["TCGA"]], all_data[["OAK_POPLAR"]])
  
  # 确保Dataset是因子类型且级别正确
  integrated_data$Dataset <- factor(integrated_data$Dataset, levels = c("TCGA", "OAK_POPLAR"))
  
  # 准备热图数据 - 按分组和数据集计算平均值
  cell_cols <- setdiff(colnames(integrated_data), c("Sample", "group", "Dataset"))
  
  # 重新排序细胞类型列
  available_cell_types <- intersect(cell_type_order, cell_cols)
  other_cell_types <- setdiff(cell_cols, cell_type_order)
  ordered_cell_cols <- c(available_cell_types, other_cell_types)
  
  # 修改：按组别分组，展示不同数据集的细胞类型分布
  # 首先统一组别名称
  common_groups <- c("cluster_low", "cluster_medium", "cluster_high")
  
  # 过滤只保留共同的组别
  filtered_data <- integrated_data %>%
    dplyr::filter(group %in% common_groups)
  
  # 按组别和数据集计算平均值
  heatmap_data <- filtered_data %>%
    dplyr::select(Dataset, group, dplyr::all_of(ordered_cell_cols)) %>%
    dplyr::group_by(group, Dataset) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(ordered_cell_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::mutate(Group_Dataset = paste(group, Dataset, sep = "_"))
  
  # 创建热图矩阵
  heatmap_matrix <- as.matrix(heatmap_data[, ordered_cell_cols])
  rownames(heatmap_matrix) <- heatmap_data$Group_Dataset
  
  # 创建注释
  annotation_df <- data.frame(
    Dataset = heatmap_data$Dataset,
    Group = heatmap_data$group,
    row.names = heatmap_data$Group_Dataset
  )
  
  # 确保注释颜色与因子级别匹配
  dataset_levels <- levels(annotation_df$Dataset)
  group_levels <- unique(annotation_df$Group)
  
  # 创建注释颜色 - 确保所有级别都有对应颜色
  dataset_color_palette <- dataset_colors[dataset_levels]
  
  # 使用预定义的组别颜色
  group_color_palette <- c(
    "cluster_low" = "#00468B",
    "cluster_medium" = "#ADB6CA",
    "cluster_high" = "#FDAF91"
  )
  
  # 确保所有组别都有颜色
  if(!all(group_levels %in% names(group_color_palette))) {
    missing_groups <- setdiff(group_levels, names(group_color_palette))
    additional_colors <- viridis(length(missing_groups))
    names(additional_colors) <- missing_groups
    group_color_palette <- c(group_color_palette, additional_colors)
  }
  
  annotation_colors <- list(
    Dataset = dataset_color_palette,
    Group = group_color_palette
  )
  
  # 绘制热图
  pdf_file <- file.path(heatmap_dir, "heatmap_TCGA_OAK_POPLAR_integrated.pdf")
  pdf(pdf_file, width = 12, height = 8)
  
  pheatmap(
    heatmap_matrix,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    main = "Cell Type Proportions Across Groups and Datasets",
    color = colorRampPalette(c(
      "#053061", "#2166ac", "#4393c3", "#92c5de", 
      "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", 
      "#d6604d", "#b2182b", "#67001f"
    ))(100),
    scale = "column",  # 按列进行标准化
    cluster_rows = FALSE,  # 不对行进行聚类，保持组别顺序
    cluster_cols = FALSE,  # 不对列进行聚类，保持细胞类型顺序
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 14,
    fontsize_row = 12,
    fontsize_col = 14,
    border_color = NA,
    cellwidth = 18,
    cellheight = 18,
    angle_col = 90  # 倾斜列名以防止重叠
  )
  
  dev.off()
  
  # 保存PNG版本
  png_file <- file.path(heatmap_dir, "heatmap_TCGA_OAK_POPLAR_integrated.png")
  png(png_file, width = 12, height = 8, units = "in", res = 300)
  
  pheatmap(
    heatmap_matrix,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    main = "Cell Type Proportions Across Groups and Datasets",
    color = colorRampPalette(c(
      "#053061", "#2166ac", "#4393c3", "#92c5de", 
      "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", 
      "#d6604d", "#b2182b", "#67001f"
    ))(100),
    scale = "column",  # 按列进行标准化
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 14,
    fontsize_row = 12,
    fontsize_col = 14,
    border_color = NA,
    cellwidth = 18,
    cellheight = 18,
    angle_col = 90  # 倾斜列名以防止重叠
  )
  
  dev.off()
}

# 4. 绘制折线图函数 - 修改以使用新的细胞类型顺序和颜色，增强视觉效果
plot_line_graph <- function(data_df, cell_types, title, output_file, colors = NULL) {
  # 确保所有需要的细胞类型都存在
  available_cell_types <- intersect(cell_types, colnames(data_df))
  
  if(length(available_cell_types) == 0) {
    cat("警告：在数据中找不到任何指定的细胞类型，无法绘制图表\n")
    return(NULL)
  }
  
  # 数据处理：计算每个组和细胞类型的平均比例
  summary_data <- data_df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(available_cell_types), 
                        names_to = "Cell_Type", values_to = "Proportion") %>%
    dplyr::group_by(group, Cell_Type) %>%
    dplyr::summarise(
      mean_proportion = mean(Proportion, na.rm = TRUE),
      se_proportion = sd(Proportion, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::ungroup()
  
  # 确保细胞类型因子级别正确 - 使用新的顺序
  ordered_cell_types <- intersect(cell_type_order, available_cell_types)
  other_cell_types <- setdiff(available_cell_types, cell_type_order)
  all_ordered_types <- c(ordered_cell_types, other_cell_types)
  
  summary_data$Cell_Type <- factor(summary_data$Cell_Type, 
                                   levels = all_ordered_types)
  
  # 调整颜色调色板
  if(is.null(colors)) {
    # 使用预定义的颜色
    available_colors <- main_colors[names(main_colors) %in% levels(summary_data$Cell_Type)]
    
    # 如果有缺失的颜色，分配新颜色
    if(length(available_colors) < length(levels(summary_data$Cell_Type))) {
      missing_cell_types <- setdiff(levels(summary_data$Cell_Type), names(available_colors))
      additional_colors <- viridis(length(missing_cell_types))
      names(additional_colors) <- missing_cell_types
      available_colors <- c(available_colors, additional_colors)
    }
  } else {
    available_colors <- colors
  }
  
  # 绘制CNS风格图形
  p <- ggplot(summary_data, 
              aes(x = group, y = mean_proportion, 
                  color = Cell_Type, 
                  group = Cell_Type)) +
    # 折线 - 使用linewidth替代size
    geom_line(linewidth = 1.5) +  # 增加线宽
    
    # 数据点
    geom_point(
      size = 4.5,  # 增加点大小
      shape = 16,  # 实心圆点
      stroke = 1.2  # 增加描边宽度
    ) +
    
    # 误差线
    geom_errorbar(
      aes(ymin = mean_proportion - se_proportion, 
          ymax = mean_proportion + se_proportion),
      width = 0.2, 
      linewidth = 1.0  # 增加误差线宽度
    ) +
    
    # 颜色和主题
    scale_color_manual(values = available_colors) +
    
    # 专业的主题
    theme_minimal() +
    
    # 标签
    labs(
      title = title,
      x = "Group", 
      y = "Mean Proportion",
      caption = "BayesPrism Analysis"
    ) +
    
    # 精细主题设置
    theme(
      # 标题设置
      plot.title = element_text(
        face = "bold", 
        size = 20,
        hjust = 0.5,
        margin = margin(0, 0, 10, 0)
      ),
      plot.caption = element_text(
        size = 16,
        color = "grey50", 
        hjust = 1,
        margin = margin(10, 0, 0, 0)
      ),
      
      # 坐标轴设置
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(
        angle = 0,  # 水平显示
        hjust = 0.5, 
        size = 16,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 16,
        color = "black"
      ),
      
      # 图例设置
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1, "cm"),
      
      # 背景设置
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
    )
  
  # 保存图形
  ggsave(
    plot = p, 
    filename = output_file, 
    width = 7,
    height = 5,
    dpi = 300
  )
  
  # 同时保存PNG版本
  png_file <- gsub("\\.pdf$", ".png", output_file)
  ggsave(
    plot = p, 
    filename = png_file, 
    width = 7,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

# 5. 绘制数据集折线图
plot_dataset_lines <- function(all_data) {
  cat("绘制各数据集细胞比例折线图...\n")
  
  for (dataset_name in names(all_data)) {
    data_df <- all_data[[dataset_name]]
    
    # 获取可用的细胞类型
    all_cell_types <- setdiff(colnames(data_df), c("Sample", "group", "Dataset"))
    
    # 1. 癌症上皮细胞折线图
    cancer_cells_available <- intersect(cancer_epi, all_cell_types)
    
    if(length(cancer_cells_available) > 0) {
      output_file <- file.path(line_dir, paste0("line_cancer_", dataset_name, ".pdf"))
      plot_line_graph(
        data_df,
        cancer_cells_available,
        paste(dataset_name, "- Cancer Epithelial Cell Proportions"),
        output_file
      )
    }
    
    # 2. 免疫细胞折线图
    immune_cells_available <- intersect(immu, all_cell_types)
    
    if(length(immune_cells_available) > 0) {
      output_file <- file.path(line_dir, paste0("line_immune_", dataset_name, ".pdf"))
      plot_line_graph(
        data_df,
        immune_cells_available,
        paste(dataset_name, "- Immune Cell Proportions"),
        output_file
      )
    }
    
    # 3. 所有细胞类型折线图
    if(length(all_cell_types) > 0) {
      output_file <- file.path(line_dir, paste0("line_all_", dataset_name, ".pdf"))
      plot_line_graph(
        data_df,
        all_cell_types,
        paste(dataset_name, "- All Cell Type Proportions"),
        output_file
      )
    }
    
    # 4. 癌症上皮细胞折线图
    cancer_cells_available <- intersect(type1, all_cell_types)
    
    if(length(cancer_cells_available) > 0) {
      output_file <- file.path(line_dir, paste0("line_SCLC-A_", dataset_name, ".pdf"))
      plot_line_graph(
        data_df,
        cancer_cells_available,
        paste(dataset_name, "- Cancer Epithelial Cell Proportions"),
        output_file
      )
    }
    # 5. 癌症上皮细胞折线图
    cancer_cells_available <- intersect(type2, all_cell_types)
    
    if(length(cancer_cells_available) > 0) {
      output_file <- file.path(line_dir, paste0("line_SCLC-N_", dataset_name, ".pdf"))
      plot_line_graph(
        data_df,
        cancer_cells_available,
        paste(dataset_name, "- Cancer Epithelial Cell Proportions"),
        output_file
      )
    }
    
  }
}

# 8. 绘制堆叠柱状图函数 - 修改以使用新的细胞类型顺序和颜色，增强视觉效果
plot_stacked_bars <- function(all_data) {
  cat("绘制堆叠柱状图...\n")
  
  for (dataset_name in names(all_data)) {
    data_df <- all_data[[dataset_name]]
    
    # 1. 所有细胞类型
    all_cell_types <- setdiff(colnames(data_df), c("Sample", "group", "Dataset"))
    
    # 重新排序细胞类型
    ordered_cell_types <- intersect(cell_type_order, all_cell_types)
    other_cell_types <- setdiff(all_cell_types, cell_type_order)
    all_ordered_types <- c(ordered_cell_types, other_cell_types)
    
    # 计算每个组的平均细胞比例
    summary_data <- data_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(all_cell_types), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
      tidyr::pivot_longer(cols = -group, names_to = "Cell_Type", values_to = "Proportion")
    
    # 确保细胞类型因子级别正确
    summary_data$Cell_Type <- factor(summary_data$Cell_Type, levels = all_ordered_types)
    
    # 准备颜色
    available_colors <- main_colors[names(main_colors) %in% all_ordered_types]
    
    # 如果有缺失的颜色，分配新颜色
    if(length(available_colors) < length(all_ordered_types)) {
      missing_cell_types <- setdiff(all_ordered_types, names(available_colors))
      additional_colors <- viridis(length(missing_cell_types))
      names(additional_colors) <- missing_cell_types
      available_colors <- c(available_colors, additional_colors)
    }
    
    # 绘制堆叠柱状图 - 增强视觉效果
    p <- ggplot(summary_data, aes(x = group, y = Proportion, fill = Cell_Type)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_fill_manual(values = available_colors) +
      labs(
        title = paste(dataset_name, "- Cell Type Composition by Group"),
        x = "Group",
        y = "Mean Proportion",
        fill = "Cell Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.8, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
      )
    
    # 保存图形
    ggsave(
      plot = p,
      filename = file.path(bar_dir, paste0("stacked_bar_all_", dataset_name, ".pdf")),
      width = 10,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      plot = p,
      filename = file.path(bar_dir, paste0("stacked_bar_all_", dataset_name, ".png")),
      width = 10,
      height = 7,
      dpi = 300
    )
    
    # 2. 癌症上皮细胞
    cancer_cells_available <- intersect(cancer_epi, all_cell_types)
    
    if(length(cancer_cells_available) > 0) {
      # 按照新的顺序排列
      cancer_cells_ordered <- intersect(cell_type_order, cancer_cells_available)
      
      # 计算每个组的平均细胞比例
      summary_data <- data_df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(cancer_cells_available), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
        tidyr::pivot_longer(cols = -group, names_to = "Cell_Type", values_to = "Proportion")
      
      # 确保细胞类型因子级别正确
      summary_data$Cell_Type <- factor(summary_data$Cell_Type, levels = cancer_cells_ordered)
      
      # 绘制堆叠柱状图
      p <- ggplot(summary_data, aes(x = group, y = Proportion, fill = Cell_Type)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7) +
        scale_fill_manual(values = main_colors[cancer_cells_ordered]) +
        labs(
          title = paste(dataset_name, "- Cancer Epithelial Cell Composition"),
          x = "Group",
          y = "Mean Proportion",
          fill = "Cell Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(bar_dir, paste0("stacked_bar_cancer_", dataset_name, ".pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(bar_dir, paste0("stacked_bar_cancer_", dataset_name, ".png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
    
    # 3. 免疫细胞
    immune_cells_available <- intersect(immu, all_cell_types)
    
    if(length(immune_cells_available) > 0) {
      # 按照新的顺序排列
      immune_cells_ordered <- intersect(cell_type_order, immune_cells_available)
      
      # 计算每个组的平均细胞比例
      summary_data <- data_df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(immune_cells_available), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
        tidyr::pivot_longer(cols = -group, names_to = "Cell_Type", values_to = "Proportion")
      
      # 确保细胞类型因子级别正确
      summary_data$Cell_Type <- factor(summary_data$Cell_Type, levels = immune_cells_ordered)
      
      # 绘制堆叠柱状图
      p <- ggplot(summary_data, aes(x = group, y = Proportion, fill = Cell_Type)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7) +
        scale_fill_manual(values = main_colors[immune_cells_ordered]) +
        labs(
          title = paste(dataset_name, "- Immune Cell Composition"),
          x = "Group",
          y = "Mean Proportion",
          fill = "Cell Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(bar_dir, paste0("stacked_bar_immune_", dataset_name, ".pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(bar_dir, paste0("stacked_bar_immune_", dataset_name, ".png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }
}

# 堆叠面积图函数 - 替代堆叠柱形折线图
plot_stacked_area <- function(all_data) {
  cat("绘制堆叠面积图...\n")
  
  for (dataset_name in names(all_data)) {
    data_df <- all_data[[dataset_name]]
    
    # 获取可用的细胞类型
    all_cell_types <- setdiff(colnames(data_df), c("Sample", "group", "Dataset"))
    
    # 重新排序细胞类型
    ordered_cell_types <- intersect(cell_type_order, all_cell_types)
    other_cell_types <- setdiff(all_cell_types, cell_type_order)
    all_ordered_types <- c(ordered_cell_types, other_cell_types)
    
    # 1. 所有细胞类型
    # 计算每个组的平均细胞比例
    summary_data <- data_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(all_cell_types), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
    
    # 转换为长格式
    stacked_data <- summary_data %>%
      tidyr::pivot_longer(cols = -group, names_to = "Cell_Type", values_to = "Proportion") %>%
      dplyr::mutate(Cell_Type = factor(Cell_Type, levels = all_ordered_types))
    
    # 准备颜色
    available_colors <- main_colors[names(main_colors) %in% all_ordered_types]
    
    # 如果有缺失的颜色，分配新颜色
    if(length(available_colors) < length(all_ordered_types)) {
      missing_cell_types <- setdiff(all_ordered_types, names(available_colors))
      additional_colors <- viridis(length(missing_cell_types))
      names(additional_colors) <- missing_cell_types
      available_colors <- c(available_colors, additional_colors)
    }
    
    # 绘制堆叠面积图
    p <- ggplot(stacked_data, aes(x = group, y = Proportion, fill = Cell_Type)) +
      geom_area(position = "stack", alpha = 0.9) +
      scale_fill_manual(values = available_colors) +
      labs(
        title = paste(dataset_name, "- Cell Type Composition Trends"),
        x = "Group",
        y = "Mean Proportion",
        fill = "Cell Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.8, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
      )
    
    # 保存图形
    ggsave(
      plot = p,
      filename = file.path(area_dir, paste0("stacked_area_all_", dataset_name, ".pdf")),
      width = 10,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      plot = p,
      filename = file.path(area_dir, paste0("stacked_area_all_", dataset_name, ".png")),
      width = 10,
      height = 7,
      dpi = 300
    )
    
    # 2. 癌症上皮细胞
    cancer_cells_available <- intersect(cancer_epi, all_cell_types)
    
    if(length(cancer_cells_available) > 0) {
      # 按照新的顺序排列
      cancer_cells_ordered <- intersect(cell_type_order, cancer_cells_available)
      
      # 计算每个组的平均细胞比例
      summary_data <- data_df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(cancer_cells_available), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
      
      # 转换为长格式
      stacked_data <- summary_data %>%
        tidyr::pivot_longer(cols = -group, names_to = "Cell_Type", values_to = "Proportion") %>%
        dplyr::mutate(Cell_Type = factor(Cell_Type, levels = cancer_cells_ordered))
      
      # 绘制堆叠面积图
      p <- ggplot(stacked_data, aes(x = group, y = Proportion, fill = Cell_Type)) +
        geom_area(position = "stack", alpha = 0.9) +
        scale_fill_manual(values = main_colors[cancer_cells_ordered]) +
        labs(
          title = paste(dataset_name, "- Cancer Cell Composition Trends"),
          x = "Group",
          y = "Mean Proportion",
          fill = "Cell Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(area_dir, paste0("stacked_area_cancer_", dataset_name, ".pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(area_dir, paste0("stacked_area_cancer_", dataset_name, ".png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
    
    # 3. 免疫细胞
    immune_cells_available <- intersect(immu, all_cell_types)
    
    if(length(immune_cells_available) > 0) {
      # 按照新的顺序排列
      immune_cells_ordered <- intersect(cell_type_order, immune_cells_available)
      
      # 计算每个组的平均细胞比例
      summary_data <- data_df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(immune_cells_available), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
      
      # 转换为长格式
      stacked_data <- summary_data %>%
        tidyr::pivot_longer(cols = -group, names_to = "Cell_Type", values_to = "Proportion") %>%
        dplyr::mutate(Cell_Type = factor(Cell_Type, levels = immune_cells_ordered))
      
      # 绘制堆叠面积图
      p <- ggplot(stacked_data, aes(x = group, y = Proportion, fill = Cell_Type)) +
        geom_area(position = "stack", alpha = 0.9) +
        scale_fill_manual(values = main_colors[immune_cells_ordered]) +
        labs(
          title = paste(dataset_name, "- Immune Cell Composition Trends"),
          x = "Group",
          y = "Mean Proportion",
          fill = "Cell Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(area_dir, paste0("stacked_area_immune_", dataset_name, ".pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(area_dir, paste0("stacked_area_immune_", dataset_name, ".png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }
}
# 9. 绘制组间比较的箱线图 - 增强视觉效果以达到CNS发表级别
plot_boxplot_comparisons <- function(all_data) {
  cat("绘制组间比较箱线图...\n")
  
  for (dataset_name in names(all_data)) {
    data_df <- all_data[[dataset_name]]
    
    # 获取可用的细胞类型
    all_cell_types <- setdiff(colnames(data_df), c("Sample", "group", "Dataset"))
    
    # 1. 癌症上皮细胞比较
    cancer_cells_available <- intersect(cancer_epi, all_cell_types)
    cancer_cells_ordered <- intersect(cell_type_order, cancer_cells_available)
    
    if(length(cancer_cells_available) > 0) {
      # 转换为长格式
      long_data <- data_df %>%
        dplyr::select(Sample, group, dplyr::all_of(cancer_cells_available)) %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(cancer_cells_available),
          names_to = "Cell_Type",
          values_to = "Proportion"
        ) %>%
        dplyr::mutate(Cell_Type = factor(Cell_Type, levels = cancer_cells_ordered))
      
      # 绘制箱线图
      p <- ggplot(
        long_data,
        aes(x = group, y = Proportion, fill = group)
      ) +
        geom_boxplot(
          width = 0.7,
          alpha = 0.8,
          outlier.shape = 16,
          outlier.size = 2,
          outlier.alpha = 0.7,
          linewidth = 0.6
        ) +
        scale_fill_manual(values = group_colors[[dataset_name]]) +
        facet_wrap(~ Cell_Type, scales = "free_y", ncol = 2) +
        labs(
          title = paste(dataset_name, "- Cancer Epithelial Cell Proportions by Group"),
          x = "Group",
          y = "Cell Proportion",
          fill = "Group",
          caption = "BayesPrism Analysis"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 12, color = "grey50", hjust = 1),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 添加组间比较 - 显示所有组间比较
      if(length(levels(data_df$group)) > 1) {
        p <- p + stat_compare_means(
          aes(label = ..p.signif..),
          method = "wilcox.test",
          comparisons = combn(levels(data_df$group), 2, simplify = FALSE),
          label.y.npc = "top",
          step.increase = 0.1,
          hide.ns = TRUE,
          size = 4,
          vjust = 0.5
        )
      }
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(box_dir, paste0("boxplot_cancer_", dataset_name, ".pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(box_dir, paste0("boxplot_cancer_", dataset_name, ".png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
    
    # 2. 免疫细胞比较
    immune_cells_available <- intersect(immu, all_cell_types)
    immune_cells_ordered <- intersect(cell_type_order, immune_cells_available)
    
    if(length(immune_cells_available) > 0) {
      # 转换为长格式
      long_data <- data_df %>%
        dplyr::select(Sample, group, dplyr::all_of(immune_cells_available)) %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(immune_cells_available),
          names_to = "Cell_Type",
          values_to = "Proportion"
        ) %>%
        dplyr::mutate(Cell_Type = factor(Cell_Type, levels = immune_cells_ordered))
      
      # 绘制箱线图
      p <- ggplot(
        long_data,
        aes(x = group, y = Proportion, fill = group)
      ) +
        geom_boxplot(
          width = 0.7,
          alpha = 0.8,
          outlier.shape = 16,
          outlier.size = 2,
          outlier.alpha = 0.7,
          linewidth = 0.6
        ) +
        scale_fill_manual(values = group_colors[[dataset_name]]) +
        facet_wrap(~ Cell_Type, scales = "free_y", ncol = 2) +
        labs(
          title = paste(dataset_name, "- Immune Cell Proportions by Group"),
          x = "Group",
          y = "Cell Proportion",
          fill = "Group",
          caption = "BayesPrism Analysis"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 12, color = "grey50", hjust = 1),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 添加组间比较 - 显示所有组间比较
      if(length(levels(data_df$group)) > 1) {
        p <- p + stat_compare_means(
          aes(label = ..p.signif..),
          method = "wilcox.test",
          comparisons = combn(levels(data_df$group), 2, simplify = FALSE),
          label.y.npc = "top",
          step.increase = 0.1,
          hide.ns = TRUE,
          size = 4,
          vjust = 0.5
        )
      }
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(box_dir, paste0("boxplot_immune_", dataset_name, ".pdf")),
        width = 12,
        height = 8,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(box_dir, paste0("boxplot_immune_", dataset_name, ".png")),
        width = 12,
        height = 8,
        dpi = 300
      )
    }
  }
}

# 6. 绘制数据集间比较箱线图 - 增强视觉效果
plot_comparative_boxplots <- function(all_data) {
  cat("绘制数据集间细胞比例比较图...\n")
  
  # 合并所有数据集
  combined_data <- dplyr::bind_rows(all_data)
  
  # 确保数据集和分组顺序正确
  combined_data$Dataset <- factor(combined_data$Dataset, levels = names(all_data))
  
  # 获取所有可用的细胞类型
  all_cell_types <- setdiff(colnames(combined_data), c("Sample", "group", "Dataset"))
  
  # 1. 癌症上皮细胞比较
  cancer_cells_available <- intersect(cancer_epi, all_cell_types)
  cancer_cells_ordered <- intersect(cell_type_order, cancer_cells_available)
  
  if(length(cancer_cells_available) > 0) {
    for(cell_type in cancer_cells_ordered) {
      p <- ggplot(
        combined_data,
        aes(x = Dataset, y = .data[[cell_type]], fill = group)
      ) +
        geom_boxplot(
          alpha = 0.8, 
          outlier.shape = 16, 
          outlier.size = 2,
          outlier.alpha = 0.7,
          position = position_dodge(width = 0.8),
          width = 0.7,
          linewidth = 0.6
        ) +
        scale_fill_manual(values = unlist(lapply(names(all_data), function(ds) {
          group_colors[[ds]][levels(combined_data$group[combined_data$Dataset == ds])]
        }))) +
        labs(
          title = paste("Comparison of", cell_type, "Proportions Across Datasets"),
          x = "Dataset",
          y = "Cell Proportion",
          fill = "Group",
          caption = "BayesPrism Analysis"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 12, color = "grey50", hjust = 1),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(comparison_dir, paste0("boxplot_", cell_type, "_comparison.pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(comparison_dir, paste0("boxplot_", cell_type, "_comparison.png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }
  
  # 2. 免疫细胞比较
  immune_cells_available <- intersect(immu, all_cell_types)
  immune_cells_ordered <- intersect(cell_type_order, immune_cells_available)
  
  if(length(immune_cells_available) > 0) {
    for(cell_type in immune_cells_ordered) {
      p <- ggplot(
        combined_data,
        aes(x = Dataset, y = .data[[cell_type]], fill = group)
      ) +
        geom_boxplot(
          alpha = 0.8, 
          outlier.shape = 16, 
          outlier.size = 2,
          outlier.alpha = 0.7,
          position = position_dodge(width = 0.8),
          width = 0.7,
          linewidth = 0.6
        ) +
        scale_fill_manual(values = unlist(lapply(names(all_data), function(ds) {
          group_colors[[ds]][levels(combined_data$group[combined_data$Dataset == ds])]
        }))) +
        labs(
          title = paste("Comparison of", cell_type, "Proportions Across Datasets"),
          x = "Dataset",
          y = "Cell Proportion",
          fill = "Group",
          caption = "BayesPrism Analysis"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 12, color = "grey50", hjust = 1),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
        )
      
      # 保存图形
      ggsave(
        plot = p,
        filename = file.path(comparison_dir, paste0("boxplot_", cell_type, "_comparison.pdf")),
        width = 10,
        height = 7,
        dpi = 300
      )
      
      ggsave(
        plot = p,
        filename = file.path(comparison_dir, paste0("boxplot_", cell_type, "_comparison.png")),
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }
}

# 创建汇总报告函数 - 增强报告格式
create_summary_report <- function(all_data) {
  cat("创建汇总报告...\n")
  
  report_file <- file.path(vis_dir, "summary_report.txt")
  
  # 打开文件连接
  con <- file(report_file, "w")
  
  # 写入报告标题
  cat("=================================================\n", file = con)
  cat("       BayesPrism Deconvolution Analysis Summary\n", file = con)
  cat("=================================================\n\n", file = con)
  
  # 写入生成日期
  cat("Date Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = con)
  
  # 写入数据集信息
  cat("Analyzed Datasets:\n", file = con)
  for (dataset_name in names(all_data)) {
    data_df <- all_data[[dataset_name]]
    cat("- ", dataset_name, ": ", nrow(data_df), " samples, ", 
        length(unique(data_df$group)), " groups\n", file = con)
    cat("  Groups: ", paste(levels(data_df$group), collapse = ", "), "\n", file = con)
    
    # 获取可用的细胞类型
    all_cell_types <- setdiff(colnames(data_df), c("Sample", "group", "Dataset"))
    cat("  Detected Cell Types: ", length(all_cell_types), "\n", file = con)
    cat("  Cell Types: ", paste(all_cell_types, collapse = ", "), "\n\n", file = con)
    
    # 计算每个组的平均细胞比例
    summary_data <- data_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(all_cell_types), \(x) mean(x, na.rm = TRUE) * 100), .groups = "drop")
    
    # 写入每个组的主要细胞类型
    cat("  Major Cell Type Proportions (%):\n", file = con)
    for (g in levels(data_df$group)) {
      cat("    ", g, ":\n", file = con)
      
      # 获取该组的细胞比例
      group_data <- summary_data %>% dplyr::filter(group == g)
      
      # 按比例排序并获取前5个细胞类型
      cell_proportions <- sort(as.numeric(group_data[1, all_cell_types]), decreasing = TRUE)
      top_cells <- names(cell_proportions)[1:min(5, length(cell_proportions))]
      
      for (cell in top_cells) {
        cat("      - ", cell, ": ", sprintf("%.2f", group_data[[cell]]), "%\n", file = con)
      }
      cat("\n", file = con)
    }
  }
  
  # 写入可视化结果摘要
  cat("Generated Visualizations:\n", file = con)
  cat("- Heatmaps: Showing cell type proportions for each dataset (in ", heatmap_dir, ")\n", file = con)
  cat("- Line Plots: Showing cell type proportion trends across groups (in ", line_dir, ")\n", file = con)
  cat("- Stacked Bar Plots: Showing cell type compositions across groups (in ", bar_dir, ")\n", file = con)
  cat("- Stacked Area Plots: Showing cell type composition trends (in ", area_dir, ")\n", file = con)
  cat("- Box Plots: Comparing cell type proportions between groups (in ", box_dir, ")\n", file = con)
  cat("- Comparative Box Plots: Comparing cell types across datasets (in ", comparison_dir, ")\n", file = con)
  
  # 关闭文件连接
  close(con)
  
  cat("Summary report saved to:", report_file, "\n")
}

# 10. 主函数
main <- function() {
  # 准备数据
  all_data <- prepare_data()
  
  if(length(all_data) == 0) {
    stop("No valid deconvolution result data found.")
  }
  
  # 执行可视化
  plot_dataset_heatmaps(all_data)
  
  # 如果有TCGA和OAK_POPLAR数据，绘制整合热图
  if(all(c("TCGA", "OAK_POPLAR") %in% names(all_data))) {
    plot_integrated_heatmap(all_data)
  }
  
  # 绘制折线图
  plot_dataset_lines(all_data)
  
  # 绘制堆叠柱状图
  plot_stacked_bars(all_data)
  
  # 绘制堆叠面积图
  plot_stacked_area(all_data)
  
  # 绘制组间比较的箱线图
  plot_boxplot_comparisons(all_data)
  
  # 绘制数据集间比较图
  plot_comparative_boxplots(all_data)
  
  # 创建汇总报告
  create_summary_report(all_data)
  
  cat("\nVisualization completed! Results saved in", vis_dir, "directory.\n")
}

# 执行主程序
tryCatch({
  main()
}, error = function(e) {
  cat("Visualization error:", conditionMessage(e), "\n")
  # 记录错误到日志文件
  log_file <- file.path(vis_dir, "visualization_log.txt")
  cat("Visualization error:", conditionMessage(e), "\n", file = log_file, append = TRUE)
})



