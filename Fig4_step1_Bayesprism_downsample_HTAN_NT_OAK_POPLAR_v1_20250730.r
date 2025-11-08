#!/usr/bin/env Rscript

###
### 这个脚本修复了TCGA和OAK-POPLAR数据集的处理问题，并增加了热图可视化
###

# 创建输出文件夹和临时目录
output_dir <- "/home/data/tmh_project/SCLC/Fig4/1_HTAN_downsample_bayesprism"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# 创建临时目录
temp_dir <- file.path(output_dir, "tmp")
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = temp_dir)
options(future.globals.maxSize = 8 * 1024^3)  # 8GB

###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(BayesPrism)
library(scales)
library(ggplot2)
library(ggpubr)  # 用于stat_compare_means函数
library(pheatmap) # 用于绘制热图
library(RColorBrewer) # 用于热图配色

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)

# 定义数据集信息
bulk_data_info <- list(
  NT = list(
    path = "/home/data/tmh_project/SCLC/Fig2_NTbulk_human/1_NTbulkdata_E-MTAB-10399/GeneCounts_new.csv",
    type = "counts",
    metadata = "/home/data/tmh_project/SCLC/Fig2_NTbulk_human/1_NTbulkdata_E-MTAB-10399/metadata_all_48samples.csv"
  ),
  TCGA = list(
    path = '/home/data/tmh_project/SCLC/Fig3_cohort_bulk_human/0_Fig3_cohort_data/1_TCGA_LUAD/LUAD_fpkm_new.csv',
    type = "fpkm",
    metadata = NULL  # 将从clustering_output中获取
  ),
  OAK_POPLAR = list(
    path = '/home/data/tmh_project/SCLC/Fig3_cohort_bulk_human/0_Fig3_cohort_data/5_OAK_POPLAR/OAK_POPLAR_NOS_tpm.csv',
    type = "tpm",
    metadata = NULL  # 将从predictions_OAK中获取
  )
)

# 读取TCGA元数据
tryCatch({
  clustering_output <- readRDS('/home/data/tmh_project/SCLC/Fig3_cohort_bulk_human/1_TCGA/ssgsea/clustering_output_complete.rds')
  meta_tcga <- clustering_output$cluster_results$hierarchical
  meta_tcga <- data.frame(
    Sample = rownames(meta_tcga),
    group = meta_tcga$cluster,
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  cat("警告：读取TCGA元数据时出错:", conditionMessage(e), "\n")
  meta_tcga <- data.frame(Sample = character(0), group = character(0))
})

# 读取OAK+POPLAR元数据
tryCatch({
  predictions_OAK <- readRDS('/home/data/tmh_project/SCLC/Fig3_cohort_bulk_human/3_predictions/OAK_POPLAR/ssgsea/hierarchical/knn/predictions_knn.rds')
  meta_oak <- data.frame(
    Sample = names(predictions_OAK$class),
    group = predictions_OAK$class,
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  cat("警告：读取OAK+POPLAR元数据时出错:", conditionMessage(e), "\n")
  meta_oak <- data.frame(Sample = character(0), group = character(0))
})

# 将FPKM转换为近似计数（简化版本，不需要基因长度信息）
fpkm_to_counts <- function(fpkm_data, library_size = 20e6, avg_gene_length = 2000) {
  # 使用固定的平均基因长度
  counts <- fpkm_data * (library_size / 1e6) * (avg_gene_length / 1000)
  return(round(counts))
}

# 将TPM转换为近似计数
tpm_to_counts <- function(tpm_data, library_size = 20e6) {
  counts <- tpm_data * (library_size / 1e6)
  return(round(counts))
}

# 安全版本的plot.bulk.vs.sc函数
safe_plot_bulk_vs_sc <- function(sc.input, bulk.input, pdf.prefix) {
  tryCatch({
    # 检查数据是否有足够的非零值进行相关性计算
    sc_nonzero <- colSums(sc.input > 0) > nrow(sc.input) * 0.01
    bulk_nonzero <- colSums(bulk.input > 0) > nrow(bulk.input) * 0.01
    common_nonzero <- names(sc_nonzero[sc_nonzero & names(sc_nonzero) %in% names(bulk_nonzero[bulk_nonzero])])
    
    if(length(common_nonzero) > 100) {
      cat("使用", length(common_nonzero), "个共同表达的基因进行一致性检查\n")
      BayesPrism::plot.bulk.vs.sc(
        sc.input = sc.input[, common_nonzero],
        bulk.input = bulk.input[, common_nonzero],
        pdf.prefix = pdf.prefix
      )
    } else {
      cat("警告：找不到足够的共同表达基因进行一致性检查，跳过此步骤\n")
    }
  }, error = function(e) {
    cat("警告：执行基因表达一致性检查时出错:", conditionMessage(e), "\n")
    cat("跳过此步骤继续执行...\n")
  })
}

# 修改后的单细胞数据处理函数，限制每个分层最多1000个细胞
process_sc_data <- function(max_cells_per_stratum = 1000, seed = 123) {
  cat("正在加载和处理单细胞数据...\n")
  
  # 检查是否已存在处理好的单细胞数据
  sc_data_file <- file.path(output_dir, "sc_data_processed.rds")
  if (file.exists(sc_data_file)) {
    cat("检测到已处理的单细胞数据，直接加载...\n")
    return(readRDS(sc_data_file))
  }
  
  tryCatch({
    sce <- qs::qread("/home/data/tmh_project/SCLC/Fig4/HTAN_LUAD_SCLC_filtered1_128979.qs")
    gc()
    
    # 提取元数据
    meta <- sce@meta.data
    
    # 查看每个分层的细胞数
    cell_counts <- table(meta$cell_type_fine, meta$histo)
    cat("原始细胞数分布:\n")
    print(cell_counts)
    
    # 设置随机种子以确保可重复性
    set.seed(seed)
    
    # 分层抽样，每个cell_type_fine和histo组合最多1000个细胞
    sampled_cells <- data.frame()
    
    for (cell_type in unique(meta$cell_type_fine)) {
      for (histo_type in c('LUAD', 'SCLC')) {
        # 获取当前分层的细胞
        current_cells <- meta[meta$cell_type_fine == cell_type & meta$histo == histo_type, ]
        
        # 如果有细胞，进行抽样
        if (nrow(current_cells) > 0) {
          # 确定抽样数量，不超过max_cells_per_stratum
          n_sample <- min(max_cells_per_stratum, nrow(current_cells))
          
          # 随机抽样
          if (n_sample < nrow(current_cells)) {
            sampled_indices <- sample(1:nrow(current_cells), n_sample)
            current_sampled <- current_cells[sampled_indices, ]
          } else {
            current_sampled <- current_cells  # 如果细胞数少于限制，全部保留
          }
          
          # 添加到结果中
          sampled_cells <- rbind(sampled_cells, current_sampled)
        }
      }
    }
    
    # 子集化Seurat对象
    sce_subset <- subset(sce, cells = rownames(sampled_cells))
    
    # 查看抽样后的细胞数
    sampled_counts <- table(sce_subset$cell_type_fine, sce_subset$histo)
    cat("\n抽样后细胞数分布:\n")
    print(sampled_counts)
    cat("\n总抽样细胞数:", ncol(sce_subset), "\n")
    
    # 保存抽样后的Seurat对象
    saveRDS(sce_subset, file.path(output_dir, "HTAN_downsampled.RDS"))
    cat("抽样后的Seurat对象已保存为 HTAN_downsampled.RDS\n")
    
    # 为BayesPrism准备数据
    sampled_cells$new <- paste(sampled_cells$cell_type_fine, sampled_cells$histo)
    sc.dat <- as.data.frame(GetAssayData(sce_subset, assay = "RNA", slot = "counts"))
    cell.state.labels <- sampled_cells$new
    cell.type.labels <- sampled_cells$cell_type_fine
    colnames(sc.dat) <- cell.state.labels
    
    # 转置矩阵以适应BayesPrism格式
    sc.dat <- t(sc.dat)
    mode(sc.dat) <- "integer"
    
    sc_data <- list(
      sc.dat = sc.dat,
      cell.type.labels = cell.type.labels,
      cell.state.labels = cell.state.labels,
      sce_subset = sce_subset
    )
    
    # 保存处理好的单细胞数据
    saveRDS(sc_data, sc_data_file)
    
    return(sc_data)
  }, error = function(e) {
    cat("处理单细胞数据时出错:", conditionMessage(e), "\n")
    stop("单细胞数据处理失败，无法继续执行")
  })
}

# 热图绘制函数
plot_heatmap <- function(data_df, output_dir, dataset_name) {
  tryCatch({
    # 提取细胞比例数据
    cell_cols <- setdiff(colnames(data_df), c("Sample", "group"))
    
    # 检查数据是否足够
    if(length(cell_cols) == 0 || nrow(data_df) == 0) {
      cat("警告：没有足够的数据绘制热图\n")
      return(NULL)
    }
    
    # 准备热图数据
    heatmap_data <- data_df %>%
      select(all_of(c("Sample", "group", cell_cols))) %>%
      arrange(group)
    
    # 提取细胞比例矩阵
    cell_prop_matrix <- as.matrix(heatmap_data[, cell_cols])
    rownames(cell_prop_matrix) <- heatmap_data$Sample
    
    # 创建注释数据
    annotation_df <- data.frame(
      Group = heatmap_data$group,
      row.names = heatmap_data$Sample
    )
    
    # 创建注释颜色
    n_groups <- length(unique(annotation_df$Group))
    group_colors <- brewer.pal(min(8, max(3, n_groups)), "Set1")[1:n_groups]
    names(group_colors) <- unique(annotation_df$Group)
    annotation_colors <- list(Group = group_colors)
    
    # 绘制热图
    pdf(file.path(output_dir, paste0("cell_proportions_heatmap_", dataset_name, ".pdf")), 
        width = 10, height = 8)
    
    pheatmap(
      cell_prop_matrix,
      annotation_row = annotation_df,
      annotation_colors = annotation_colors,
      main = paste(dataset_name, "- Cell Type Proportions"),
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "column",  # 按列进行标准化
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = ifelse(nrow(cell_prop_matrix) <= 50, TRUE, FALSE),
      show_colnames = TRUE,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      border_color = NA
    )
    
    dev.off()
    
    # 同时保存PNG版本
    png(file.path(output_dir, paste0("cell_proportions_heatmap_", dataset_name, ".png")), 
        width = 10, height = 8, units = "in", res = 300)
    
    pheatmap(
      cell_prop_matrix,
      annotation_row = annotation_df,
      annotation_colors = annotation_colors,
      main = paste(dataset_name, "- Cell Type Proportions"),
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "column",  # 按列进行标准化
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = ifelse(nrow(cell_prop_matrix) <= 50, TRUE, FALSE),
      show_colnames = TRUE,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      border_color = NA
    )
    
    dev.off()
    
    # 不聚类版本的热图
    pdf(file.path(output_dir, paste0("cell_proportions_heatmap_no_clustering_", dataset_name, ".pdf")), 
        width = 10, height = 8)
    
    pheatmap(
      cell_prop_matrix,
      annotation_row = annotation_df,
      annotation_colors = annotation_colors,
      main = paste(dataset_name, "- Cell Type Proportions (No Clustering)"),
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "column",  # 按列进行标准化
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = ifelse(nrow(cell_prop_matrix) <= 50, TRUE, FALSE),
      show_colnames = TRUE,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      border_color = NA
    )
    
    dev.off()
    
    # 同时保存PNG版本
    png(file.path(output_dir, paste0("cell_proportions_heatmap_no_clustering_", dataset_name, ".png")), 
        width = 10, height = 8, units = "in", res = 300)
    
    pheatmap(
      cell_prop_matrix,
      annotation_row = annotation_df,
      annotation_colors = annotation_colors,
      main = paste(dataset_name, "- Cell Type Proportions (No Clustering)"),
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "column",  # 按列进行标准化
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = ifelse(nrow(cell_prop_matrix) <= 50, TRUE, FALSE),
      show_colnames = TRUE,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      border_color = NA
    )
    
    dev.off()
    
    cat("热图已保存到", output_dir, "\n")
    
  }, error = function(e) {
    cat("绘制热图时出错:", conditionMessage(e), "\n")
  })
}

# 绘图函数
plot_results <- function(data_df, output_dir, dataset_name) {
  tryCatch({
    # 定义分组顺序和细胞类型顺序
    if (dataset_name == "NT") {
      group_order <- c("C_LUAD", "T_LUAD", "T_SCLC", "C_SCLC")
    } else if (dataset_name == "TCGA") {
      group_order <- unique(data_df$group)
    } else if (dataset_name == "OAK_POPLAR") {
      group_order <- unique(data_df$group)
    }
    
    cell_type_order <- c("SCLC-A", "SCLC-P", "SCLC-N", "NSCLC")
    
    # CNS风格配色
    color_palette <- c(
      "SCLC-A" = "#3B4FB4",    # 深蓝色
      "SCLC-P" = "#2CA02C",    # 深绿色
      "SCLC-N" = "#D62728",    # 深红色
      "NSCLC" = "#FF7F0E"      # 橙色
    )
    
    # 检查数据中是否包含所需的肿瘤细胞类型
    available_cell_types <- intersect(cell_type_order, colnames(data_df))
    
    if(length(available_cell_types) == 0) {
      cat("警告：在数据中找不到任何指定的肿瘤细胞类型，将使用所有可用的细胞类型\n")
      available_cell_types <- setdiff(colnames(data_df), c("Sample", "group"))
    }
    
    # 数据处理：计算每个组和细胞类型的平均比例
    summary_data <- data_df %>%
      pivot_longer(cols = all_of(available_cell_types), 
                   names_to = "Cell_Type", values_to = "Proportion") %>%
      group_by(group, Cell_Type) %>%
      summarise(
        mean_proportion = mean(Proportion, na.rm = TRUE),
        se_proportion = sd(Proportion, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      ungroup()
    
    # 确保分组和细胞类型因子级别正确
    if(all(group_order %in% summary_data$group)) {
      summary_data$group <- factor(summary_data$group, levels = group_order)
    } else {
      summary_data$group <- factor(summary_data$group)
    }
    
    if(all(cell_type_order %in% summary_data$Cell_Type)) {
      summary_data$Cell_Type <- factor(summary_data$Cell_Type, levels = intersect(cell_type_order, unique(summary_data$Cell_Type)))
    } else {
      summary_data$Cell_Type <- factor(summary_data$Cell_Type)
    }
    
    # 调整颜色调色板以匹配可用的细胞类型
    available_colors <- color_palette[names(color_palette) %in% levels(summary_data$Cell_Type)]
    if(length(available_colors) < length(levels(summary_data$Cell_Type))) {
      # 为缺少颜色的细胞类型分配新颜色
      missing_cell_types <- setdiff(levels(summary_data$Cell_Type), names(available_colors))
      additional_colors <- scales::hue_pal()(length(missing_cell_types))
      names(additional_colors) <- missing_cell_types
      available_colors <- c(available_colors, additional_colors)
    }
    
    # 绘制CNS风格图形
    p_selected <- ggplot(summary_data, 
                         aes(x = group, y = mean_proportion, 
                             color = Cell_Type, 
                             group = Cell_Type)) +
      # 折线
      geom_line(size = 1.2) +
      
      # 数据点
      geom_point(
        size = 3.5, 
        shape = 16,  # 实心圆点
        stroke = 1
      ) +
      
      # 误差线
      geom_errorbar(
        aes(ymin = mean_proportion - se_proportion, 
            ymax = mean_proportion + se_proportion),
        width = 0.15, 
        size = 0.8
      ) +
      
      # 颜色和主题
      scale_color_manual(values = available_colors) +
      
      # 专业的主题
      theme_classic() +
      
      # 标签
      labs(
        title = paste(dataset_name, "- Cell Type Proportions Across Groups"),
        x = "Group", 
        y = "Mean Proportion",
        caption = "BayesPrism Analysis"
      ) +
      
      # 精细主题设置
      theme(
        # 标题设置
        plot.title = element_text(
          face = "bold", 
          size = 14, 
          hjust = 0.5,
          margin = margin(0, 0, 10, 0)
        ),
        plot.caption = element_text(
          size = 10, 
          color = "grey50", 
          hjust = 1,
          margin = margin(10, 0, 0, 0)
        ),
        
        # 坐标轴设置
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(
          angle = 45, 
          hjust = 1, 
          size = 10,
          color = "black"
        ),
        axis.text.y = element_text(
          size = 10,
          color = "black"
        ),
        
        # 图例设置
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        
        # 移除网格线
        panel.grid = element_blank()
      )
    
    # 如果组别数量合适，添加统计显著性
    if(length(unique(summary_data$group)) >= 2 && nrow(summary_data) > 0) {
      # 尝试添加统计比较，但如果失败则继续
      tryCatch({
        # 创建比较列表
        if(length(group_order) >= 2) {
          comparisons_list <- list(c(group_order[1], group_order[2]))
          
          p_selected <- p_selected + 
            stat_compare_means(
              comparisons = comparisons_list,
              method = "t.test",
              label = "p.signif",
              size = 3,
              label.y = max(summary_data$mean_proportion, na.rm = TRUE) * 1.1
            )
        }
      }, error = function(e) {
        cat("警告：添加统计显著性标记时出错:", conditionMessage(e), "\n")
      })
    }
    
    # 保存高质量图形
    tryCatch({
      ggsave(
        plot = p_selected, 
        filename = file.path(output_dir, paste0("cell_type_proportions_", dataset_name, ".pdf")), 
        width = 8, 
        height = 6, 
        dpi = 300
      )
      
      # 同时保存高质量PNG
      ggsave(
        plot = p_selected, 
        filename = file.path(output_dir, paste0("cell_type_proportions_", dataset_name, ".png")), 
        width = 8, 
        height = 6, 
        dpi = 300
      )
    }, error = function(e) {
      cat("警告：保存图形时出错:", conditionMessage(e), "\n")
      # 尝试保存为RDS格式
      saveRDS(p_selected, file.path(output_dir, paste0("plot_object_", dataset_name, ".rds")))
    })
    
    return(p_selected)
  }, error = function(e) {
    cat("绘图过程中出错:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 创建汇总报告
create_summary_report <- function(results, output_dir) {
  tryCatch({
    report_file <- file.path(output_dir, "deconvolution_summary.txt")
    cat("BayesPrism解卷积分析汇总报告\n", file = report_file)
    cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = report_file, append = TRUE)
    
    for (dataset_name in names(results)) {
      result <- results[[dataset_name]]
      
      if(!is.null(result$theta)) {
        cat("\n===", dataset_name, "数据集 ===\n", file = report_file, append = TRUE)
        
        # 样本数量
        cat("样本数量:", nrow(result$theta), "\n", file = report_file, append = TRUE)
        
        # 细胞类型
        cat("细胞类型:", paste(colnames(result$theta), collapse = ", "), "\n", file = report_file, append = TRUE)
        
        # 细胞类型平均比例
        avg_props <- colMeans(result$theta)
        cat("平均细胞比例:\n", file = report_file, append = TRUE)
        for(cell_type in names(avg_props)) {
          cat("  ", cell_type, ":", sprintf("%.2f%%", avg_props[cell_type] * 100), "\n", file = report_file, append = TRUE)
        }
        
        # 变异系数
        if(!is.null(result$theta.cv)) {
          avg_cv <- colMeans(result$theta.cv)
          cat("平均变异系数:\n", file = report_file, append = TRUE)
          for(cell_type in names(avg_cv)) {
            cat("  ", cell_type, ":", sprintf("%.4f", avg_cv[cell_type]), "\n", file = report_file, append = TRUE)
          }
        }
        
        # 组间差异
        if(!is.null(result$data_df) && "group" %in% colnames(result$data_df)) {
          cat("组间细胞比例:\n", file = report_file, append = TRUE)
          
          # 计算每个组的平均细胞比例
          cell_types <- setdiff(colnames(result$theta), c("Sample", "group"))
          groups <- unique(result$data_df$group)
          
          for(group in groups) {
            cat("  组:", group, "\n", file = report_file, append = TRUE)
            group_data <- result$data_df[result$data_df$group == group, ]
            
            for(cell_type in cell_types) {
              if(cell_type %in% colnames(group_data)) {
                mean_prop <- mean(group_data[[cell_type]], na.rm = TRUE)
                cat("    ", cell_type, ":", sprintf("%.2f%%", mean_prop * 100), "\n", file = report_file, append = TRUE)
              }
            }
          }
        }
      } else {
        cat("\n===", dataset_name, "数据集 ===\n", file = report_file, append = TRUE)
        cat("处理失败或无结果\n", file = report_file, append = TRUE)
      }
    }
    
    cat("\n汇总报告已生成:", report_file, "\n")
  }, error = function(e) {
    cat("创建汇总报告时出错:", conditionMessage(e), "\n")
  })
}

# 绘制所有数据集的综合比较图
plot_combined_results <- function(results, output_dir) {
  tryCatch({
    # 准备组合数据
    combined_data <- data.frame()
    
    for (dataset_name in names(results)) {
      result <- results[[dataset_name]]
      if(!is.null(result$data_df)) {
        # 添加数据集标识
        temp_df <- result$data_df
        temp_df$Dataset <- dataset_name
        combined_data <- rbind(combined_data, temp_df)
      }
    }
    
    # 检查是否有足够的数据
    if(nrow(combined_data) == 0) {
      cat("没有足够的数据进行综合比较\n")
      return(NULL)
    }
    
    # 获取所有细胞类型
    cell_types <- setdiff(colnames(combined_data), c("Sample", "group", "Dataset"))
    # 转换为长格式数据
    long_data <- combined_data %>%
      pivot_longer(
        cols = all_of(cell_types),
        names_to = "Cell_Type",
        values_to = "Proportion"
      )
    
    # 绘制箱线图
    for(cell_type in unique(long_data$Cell_Type)) {
      tryCatch({
        p <- ggplot(
          long_data %>% filter(Cell_Type == cell_type),
          aes(x = Dataset, y = Proportion, fill = group)
        ) +
          geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
          scale_fill_brewer(palette = "Set1") +
          labs(
            title = paste("Comparison of", cell_type, "Proportions Across Datasets"),
            x = "Dataset",
            y = "Cell Proportion",
            fill = "Group"
          ) +
          theme_classic() +
          theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "right"
          )
        
        # 保存图形
        ggsave(
          plot = p,
          filename = file.path(output_dir, paste0("combined_", cell_type, "_comparison.pdf")),
          width = 10,
          height = 6,
          dpi = 300
        )
        
        ggsave(
          plot = p,
          filename = file.path(output_dir, paste0("combined_", cell_type, "_comparison.png")),
          width = 10,
          height = 6,
          dpi = 300
        )
      }, error = function(e) {
        cat("绘制", cell_type, "比较图时出错:", conditionMessage(e), "\n")
      })
    }
    
    # 绘制热图比较
    tryCatch({
      # 准备热图数据
      heatmap_data <- combined_data %>%
        select(Dataset, group, all_of(cell_types)) %>%
        group_by(Dataset, group) %>%
        summarise(across(all_of(cell_types), mean, na.rm = TRUE), .groups = "drop") %>%
        mutate(Group_Dataset = paste(group, Dataset, sep = "_"))
      
      # 创建热图矩阵
      heatmap_matrix <- as.matrix(heatmap_data[, cell_types])
      rownames(heatmap_matrix) <- heatmap_data$Group_Dataset
      
      # 创建注释
      annotation_df <- data.frame(
        Dataset = heatmap_data$Dataset,
        Group = heatmap_data$group,
        row.names = heatmap_data$Group_Dataset
      )
      
      # 创建注释颜色
      dataset_colors <- brewer.pal(min(8, length(unique(annotation_df$Dataset))), "Set2")
      names(dataset_colors) <- unique(annotation_df$Dataset)
      
      group_colors <- brewer.pal(min(8, length(unique(annotation_df$Group))), "Set1")
      names(group_colors) <- unique(annotation_df$Group)
      
      annotation_colors <- list(
        Dataset = dataset_colors,
        Group = group_colors
      )
      
      # 绘制热图
      pdf(file.path(output_dir, "combined_cell_proportions_heatmap.pdf"), 
          width = 10, height = 8)
      
      pheatmap(
        heatmap_matrix,
        annotation_row = annotation_df,
        annotation_colors = annotation_colors,
        main = "Cell Type Proportions Across Datasets and Groups",
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        scale = "column",  # 按列进行标准化
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        fontsize = 10,
        fontsize_row = 8,
        fontsize_col = 10,
        border_color = NA
      )
      
      dev.off()
      
      # PNG版本
      png(file.path(output_dir, "combined_cell_proportions_heatmap.png"), 
          width = 10, height = 8, units = "in", res = 300)
      
      pheatmap(
        heatmap_matrix,
        annotation_row = annotation_df,
        annotation_colors = annotation_colors,
        main = "Cell Type Proportions Across Datasets and Groups",
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        scale = "column",  # 按列进行标准化
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        fontsize = 10,
        fontsize_row = 8,
        fontsize_col = 10,
        border_color = NA
      )
      
      dev.off()
      
    }, error = function(e) {
      cat("绘制综合热图时出错:", conditionMessage(e), "\n")
    })
    
    cat("综合比较图已保存到", output_dir, "\n")
  }, error = function(e) {
    cat("创建综合比较图时出错:", conditionMessage(e), "\n")
  })
}

# 主函数：处理指定数据集
process_dataset <- function(dataset_name, sc_data) {
  cat("正在处理", dataset_name, "数据集...\n")
  
  # 检查是否已经处理过此数据集
  result_file <- file.path(output_dir, paste0("results_", dataset_name, ".rds"))
  if (file.exists(result_file)) {
    cat(dataset_name, "数据集已处理，直接加载结果...\n")
    return(readRDS(result_file))
  }
  
  # 获取数据集信息
  dataset_info <- bulk_data_info[[dataset_name]]
  
  # 创建输出目录
  output_subdir <- file.path(output_dir, dataset_name)
  dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
  
  # 读取bulk数据
  tryCatch({
    bk.dat <- read.csv(dataset_info$path, header = TRUE, row.names = 1)
    
    # 根据数据类型进行预处理
    if (dataset_info$type == "fpkm") {
      cat("将FPKM转换为近似计数...\n")
      bk.dat <- fpkm_to_counts(bk.dat)
    } else if (dataset_info$type == "tpm") {
      cat("将TPM转换为近似计数...\n")
      bk.dat <- tpm_to_counts(bk.dat)
    }
    
    # 转置矩阵以适应BayesPrism格式
    bk.dat <- t(bk.dat)
    mode(bk.dat) <- "integer"
  }, error = function(e) {
    cat("读取或预处理bulk数据时出错:", conditionMessage(e), "\n")
    stop("无法继续处理此数据集")
  })
  
  # 提取单细胞数据
  sc.dat <- sc_data$sc.dat
  cell.type.labels <- sc_data$cell.type.labels
  cell.state.labels <- sc_data$cell.state.labels
  
  # 对大型数据集进行特殊处理
  if(dataset_name == "TCGA" || dataset_name == "OAK_POPLAR") {
    # 对于大型数据集，减少基因数量以节省内存
    # 只保留表达量最高的15000个基因
    gene_means <- colMeans(sc.dat)
    top_genes <- names(sort(gene_means, decreasing = TRUE)[1:min(15000, length(gene_means))])
    sc.dat <- sc.dat[, top_genes]
    bk.dat <- bk.dat[, intersect(colnames(bk.dat), top_genes)]
    cat("为减少内存使用，只保留了表达量最高的", ncol(sc.dat), "个基因\n")
  }
  
  # 确保基因一致性
  common_genes <- intersect(colnames(bk.dat), colnames(sc.dat))
  cat("共有", length(common_genes), "个基因在bulk和单细胞数据中共同存在\n")
  
  bk.dat <- bk.dat[, common_genes]
  sc.dat <- sc.dat[, common_genes]
  
  # BayesPrism质控指标
  # 1. 细胞类型和状态相关性检查
  tryCatch({
    plot.cor.phi(
      input = sc.dat,
      input.labels = cell.state.labels,
      title = paste0(dataset_name, " - Cell State Correlation"),
      pdf.prefix = file.path(output_subdir, "cell_state_correlation"),
      cexRow = 0.6, 
      cexCol = 0.6,
      margins = c(6,6)
    )
  }, error = function(e) {
    cat("警告：绘制细胞状态相关性图表时出错:", conditionMessage(e), "\n")
  })
  
  tryCatch({
    plot.cor.phi(
      input = sc.dat, 
      input.labels = cell.type.labels, 
      title = paste0(dataset_name, " - Cell Type Correlation"),
      pdf.prefix = file.path(output_subdir, "cell_type_correlation"),
      cexRow = 0.6, 
      cexCol = 0.6,
      margins = c(5,5)
    )
  }, error = function(e) {
    cat("警告：绘制细胞类型相关性图表时出错:", conditionMessage(e), "\n")
  })
  
  # 2. 基因离群值检查
  sc.stat <- NULL
  tryCatch({
    sc.stat <- plot.scRNA.outlier(
      input = sc.dat,
      cell.type.labels = cell.type.labels,
      species = "hs",
      return.raw = TRUE,
      pdf.prefix = file.path(output_subdir, "scRNA_outlier")
    )
  }, error = function(e) {
    cat("警告：检查单细胞RNA离群值时出错:", conditionMessage(e), "\n")
  })
  
  bk.stat <- NULL
  tryCatch({
    bk.stat <- plot.bulk.outlier(
      bulk.input = bk.dat,
      sc.input = sc.dat,
      cell.type.labels = cell.type.labels,
      species = "hs",
      return.raw = TRUE,
      pdf.prefix = file.path(output_subdir, "bulk_outlier")
    )
  }, error = function(e) {
    cat("警告：检查bulk RNA离群值时出错:", conditionMessage(e), "\n")
  })
  
  # 3. 基因过滤
  sc.dat.filtered <- NULL
  tryCatch({
    sc.dat.filtered <- cleanup.genes(
      input = sc.dat,
      input.type = "count.matrix",
      species = "hs", 
      gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
      exp.cells = 5
    )
  }, error = function(e) {
    cat("警告：过滤基因时出错:", conditionMessage(e), "\n")
    # 如果过滤失败，使用原始数据继续
    sc.dat.filtered <- sc.dat
    cat("使用未过滤的数据继续...\n")
  })
  
  # 基因表达一致性检查
  safe_plot_bulk_vs_sc(
    sc.input = sc.dat.filtered,
    bulk.input = bk.dat,
    pdf.prefix = file.path(output_subdir, "bulk_vs_sc_plot")
  )
  
  # 4. 仅选择蛋白编码基因
  sc.dat.filtered.pc <- NULL
  tryCatch({
    sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")
    cat("成功筛选出蛋白编码基因，共", ncol(sc.dat.filtered.pc), "个基因\n")
  }, error = function(e) {
    cat("警告：筛选蛋白编码基因时出错:", conditionMessage(e), "\n")
    # 如果筛选失败，使用过滤后的数据继续
    sc.dat.filtered.pc <- sc.dat.filtered
    cat("使用过滤后的所有基因继续...\n")
  })
  
  # 5. 选择标记基因
  diff.exp.stat <- NULL
  sc.dat.filtered.pc.sig <- NULL
  tryCatch({
    diff.exp.stat <- get.exp.stat(
      sc.dat = sc.dat[, colSums(sc.dat > 0) > 3],
      cell.type.labels = cell.type.labels,
      cell.state.labels = cell.state.labels,
      pseudo.count = 0.1,
      cell.count.cutoff = 50,
      n.cores = min(20, parallel::detectCores() - 1)
    )
    
    sc.dat.filtered.pc.sig <- select.marker(
      sc.dat = sc.dat.filtered.pc,
      stat = diff.exp.stat,
      pval.max = 0.01,
      lfc.min = 0.1
    )
    cat("成功选择标记基因，共", ncol(sc.dat.filtered.pc.sig), "个基因\n")
  }, error = function(e) {
    cat("警告：选择标记基因时出错:", conditionMessage(e), "\n")
    # 如果选择标记基因失败，使用蛋白编码基因继续
    sc.dat.filtered.pc.sig <- sc.dat.filtered.pc
    cat("使用所有蛋白编码基因继续...\n")
  })
  
  # 确定要使用的单细胞数据
  sc.dat.final <- sc.dat.filtered.pc.sig
  if(is.null(sc.dat.final)) {
    sc.dat.final <- sc.dat.filtered.pc
  }
  if(is.null(sc.dat.final)) {
    sc.dat.final <- sc.dat.filtered
  }
  if(is.null(sc.dat.final)) {
    sc.dat.final <- sc.dat
  }
  
  # 再次确保基因一致性
  common_genes <- intersect(colnames(bk.dat), colnames(sc.dat.final))
  cat("最终使用", length(common_genes), "个基因进行解卷积分析\n")
  
  bk.dat <- bk.dat[, common_genes]
  sc.dat.final <- sc.dat.final[, common_genes]
  
  # BayesPrism参数调整
  gibbs_control <- list(
    chain.length = 5000,   # 增加链长度以更好地探索后验分布
    burn.in = 1000,        # 增加烧入期
    thinning = 5,          # 增加薄化间隔以减少自相关
    n.cores = min(20, parallel::detectCores() - 1),
    seed = 123
  )
  
  # 对大型数据集减少计算量
  if(dataset_name == "TCGA" || dataset_name == "OAK_POPLAR") {
    gibbs_control$chain.length <- 3000
    gibbs_control$burn.in <- 500
    cat("为加快处理速度，减少了Gibbs采样参数\n")
  }
  
  # 构建Prism对象
  myPrism <- NULL
  tryCatch({
    myPrism <- new.prism(
      reference = sc.dat.final, 
      mixture = bk.dat,
      input.type = "count.matrix", 
      cell.type.labels = cell.type.labels, 
      cell.state.labels = cell.state.labels,
      key = NULL,
      outlier.cut = 0.01,
      outlier.fraction = 0.1
    )
  }, error = function(e) {
    cat("构建Prism对象时出错:", conditionMessage(e), "\n")
    stop("无法继续解卷积分析")
  })
  
  # 运行BayesPrism
  bp.res <- NULL
  tryCatch({
    cat("运行BayesPrism解卷积...\n")
    bp.res <- run.prism(prism = myPrism, gibbs.control = gibbs_control)
    
    # 保存结果
    tryCatch({
      qs::qsave(bp.res, file.path(output_subdir, paste0("runprism_", dataset_name, ".qs")))
    }, error = function(e) {
      cat("警告：保存qs格式结果时出错:", conditionMessage(e), "\n")
    })
    
    saveRDS(bp.res, file.path(output_subdir, paste0("runprism_", dataset_name, ".rds")))
  }, error = function(e) {
    cat("运行BayesPrism时出错:", conditionMessage(e), "\n")
    stop("解卷积分析失败")
  })
  
  # 提取结果
  theta <- NULL
  theta.cv <- NULL
  tryCatch({
    theta <- get.fraction(
      bp = bp.res,
      which.theta = "final",
      state.or.type = "type"
    )
    
    theta.cv <- bp.res@posterior.theta_f@theta.cv
    
    # 保存结果
    write.csv(theta, file = file.path(output_subdir, "cell_type_proportions.csv"))
    write.csv(theta.cv, file = file.path(output_subdir, "cell_type_proportions_cv.csv"))
  }, error = function(e) {
    cat("提取或保存解卷积结果时出错:", conditionMessage(e), "\n")
    return(NULL)  # 返回NULL表示处理失败
  })
  
  # 读取元数据
  meta_df <- NULL
  tryCatch({
    if (dataset_name == "NT") {
      meta <- read.csv(dataset_info$metadata, header = TRUE, sep = ",", quote = "\"")
      meta_df <- meta[, c(1,2)]
      colnames(meta_df) <- c('Sample', 'group')
    } else if (dataset_name == "TCGA") {
      meta_df <- meta_tcga
    } else if (dataset_name == "OAK_POPLAR") {
      meta_df <- meta_oak
    }
  }, error = function(e) {
    cat("读取元数据时出错:", conditionMessage(e), "\n")
    # 创建一个基本的元数据框架
    if(!is.null(theta)) {
      meta_df <- data.frame(
        Sample = rownames(theta),
        group = "Unknown"
      )
    }
  })
  
  # 转换数据框
  data_df <- NULL
  tryCatch({
    if(!is.null(theta) && !is.null(meta_df)) {
      theta_df <- as.data.frame(theta)
      data_df <- theta_df %>%
        rownames_to_column(var = "Sample") %>%
        left_join(meta_df, by = "Sample")
      
      # 绘图和可视化
      plot_results(data_df, output_subdir, dataset_name)
      
      # 绘制热图
      plot_heatmap(data_df, output_subdir, dataset_name)
    }
  }, error = function(e) {
    cat("准备绘图数据时出错:", conditionMessage(e), "\n")
  })
  
  result <- list(theta = theta, theta.cv = theta.cv, data_df = data_df)
  saveRDS(result, result_file)
  
  return(result)
}

# 主程序
main <- function() {
  # 设置全局随机种子以确保可重复性
  set.seed(123)
  
  # 创建日志文件
  log_file <- file.path(output_dir, "deconvolution_log.txt")
  cat("开始BayesPrism解卷积分析 -", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file)
  
  # 处理单细胞数据（每个分层最多1000个细胞）
  sc_data <- NULL
  tryCatch({
    sc_data <- process_sc_data(max_cells_per_stratum = 1000, seed = 123)
    cat("单细胞数据处理成功\n", file = log_file, append = TRUE)
  }, error = function(e) {
    cat("处理单细胞数据时出错:", conditionMessage(e), "\n", file = log_file, append = TRUE)
    stop("单细胞数据处理失败，无法继续执行")
  })
  
  # 处理每个数据集
  results <- list()
  for (dataset in names(bulk_data_info)) {
    cat("\n开始处理", dataset, "数据集...\n")
    cat("\n开始处理", dataset, "数据集 -", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)
    
    tryCatch({
      results[[dataset]] <- process_dataset(dataset, sc_data)
      cat(dataset, "数据集处理完成\n")
      cat(dataset, "数据集处理完成 -", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)
    }, error = function(e) {
      cat("处理", dataset, "数据集时出错:", conditionMessage(e), "\n")
      cat("处理", dataset, "数据集时出错:", conditionMessage(e), "\n", file = log_file, append = TRUE)
      cat("跳过此数据集继续处理下一个...\n")
    })
  }
  
  # 保存所有结果
  tryCatch({
    saveRDS(results, file.path(output_dir, "all_results.rds"))
  }, error = function(e) {
    cat("保存最终结果时出错:", conditionMessage(e), "\n", file = log_file, append = TRUE)
  })
  
  # 创建汇总报告
  create_summary_report(results, output_dir)
  
  # 绘制综合比较图
  plot_combined_results(results, output_dir)
  
  cat("\n所有数据集处理完成！结果已保存在", output_dir, "目录下。\n")
  cat("\n所有数据集处理完成 -", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)
}

# 执行主程序
tryCatch({
  main()
}, error = function(e) {
  cat("主程序执行出错:", conditionMessage(e), "\n")
  # 记录错误到日志文件
  log_file <- file.path(output_dir, "deconvolution_log.txt")
  cat("主程序执行出错:", conditionMessage(e), "\n", file = log_file, append = TRUE)
})

    