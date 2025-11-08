
# 创建输出目录
output_dir <- "/home/data/tmh_project/SCLC/Fig3_multicohort_v1/6_TCGA_others"
dir.create(output_dir)

# 加载必要的库
library(ggplot2)
library(readxl)
library(dplyr)
library(stats)
library(rstatix)
library(ggpubr)

# 定义统计分析函数
perform_statistical_test <- function(data, value_column) {
  # 提取cluster_low和cluster_high的数据
  low_data <- data[data$group == "cluster_low", ][[value_column]]
  high_data <- data[data$group == "cluster_high", ][[value_column]]
  
  # 进行Shapiro-Wilk正态性检验
  sw_low <- shapiro.test(low_data)
  sw_high <- shapiro.test(high_data)
  
  # 如果两组都符合正态分布（p > 0.05），使用t检验
  if (sw_low$p.value > 0.05 && sw_high$p.value > 0.05) {
    test_result <- t.test(low_data, high_data)
    test_method <- "t-test"
  } else {
    # 否则使用Mann-Whitney U检验
    test_result <- wilcox.test(low_data, high_data)
    test_method <- "Mann-Whitney U test"
  }
  
  return(list(
    p_value = test_result$p.value,
    method = test_method
  ))
}

# 格式化p值函数
format_pvalue <- function(p_value) {
  if (p_value < 0.001) return("p < 0.001")
  if (p_value < 0.01) return(sprintf("p = %.3f", p_value))
  return(sprintf("p = %.2f", p_value))
}

# 定义TCGA颜色
tcga_colors <- c(
  "cluster_low" = "#00468B",
  "cluster_medium" = "#ADB6CA",
  "cluster_high" = "#FDAF91"
)

# 读取clustering输出数据
clustering_output <- readRDS('/home/data/tmh_project/SCLC/Fig3_multicohort_v1/1_TCGA/ssgsea/clustering_output_complete.rds')
meta_tcga <- clustering_output$cluster_results$hierarchical
meta_data <- data.frame(
  Sample = rownames(meta_tcga),
  group = meta_tcga$cluster,
  stringsAsFactors = FALSE
)
meta_data$Sample <- gsub("\\.", "-", meta_data$Sample)

# 创建子目录
violin_dir <- file.path(output_dir, "violin_plots_stats")
boxplot_dir <- file.path(output_dir, "box_plots_stats")
scatter_dir <- file.path(output_dir, "scatter_plots_stats")
density_dir <- file.path(output_dir, "density_plots_stats")

dir.create(violin_dir, showWarnings = FALSE)
dir.create(boxplot_dir, showWarnings = FALSE)
dir.create(scatter_dir, showWarnings = FALSE)
dir.create(density_dir, showWarnings = FALSE)

# 获取Excel文件中的所有sheet
sheets <- excel_sheets("0_Fig3_cohort_data/1_TCGA_LUAD/TCGA附加数据汇总.xlsx")

# 循环绘制每个sheet中的数据
for (sheet in sheets) {
  # 读取当前sheet的数据
  data <- read_excel("0_Fig3_cohort_data/1_TCGA_LUAD/TCGA附加数据汇总.xlsx", sheet = sheet)
  
  # 清理数据
  colnames(data)[colnames(data) == "SampleName"] <- "Sample"
  data <- data[, !colnames(data) %in% c("...1")]
  
  # 合并数据
  common_samples <- intersect(meta_data$Sample, data$Sample)
  
  if (length(common_samples) > 0) {
    merged_data <- merge(data, meta_data, by = "Sample", all.x = TRUE)
    merged_data <- merged_data[merged_data$Sample %in% common_samples, ]
    merged_data$group <- factor(merged_data$group, levels = c("cluster_low", "cluster_medium", "cluster_high"))
    
    # 获取数字类型的列
    numeric_columns <- sapply(merged_data, is.numeric)
    
    if (any(numeric_columns)) {
      for (col in names(merged_data)[numeric_columns]) {
        if (col %in% c("Sample", "...1")) next
        
        # 进行统计检验
        stat_test <- perform_statistical_test(
          merged_data[merged_data$group %in% c("cluster_low", "cluster_high"), ],
          col
        )
        
        # 计算y轴位置
        y_max <- max(merged_data[[col]], na.rm = TRUE)
        y_min <- min(merged_data[[col]], na.rm = TRUE)
        y_range <- y_max - y_min
        y_pos <- y_max + y_range * 0.1
        
        # 小提琴图
        p1 <- ggplot(merged_data, aes(x = group, y = .data[[col]], fill = group)) +
          geom_violin(trim = FALSE) +
          geom_boxplot(width = 0.2, fill = "white", alpha = 0.5) +
          scale_fill_manual(values = tcga_colors) +
          labs(title = paste(sheet, "-", col),
               subtitle = paste("Statistical test:", stat_test$method, "\n", format_pvalue(stat_test$p_value)),
               x = "Group",
               y = col) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 30, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 18),
            plot.subtitle = element_text(size = 12)
          )
        
        # 如果p值显著，添加显著性标记
        if (stat_test$p_value < 0.05) {
          p1 <- p1 + geom_signif(
            comparisons = list(c("cluster_low", "cluster_high")),
            y_position = y_pos,
            annotations = format_pvalue(stat_test$p_value),
            tip_length = 0.02
          )
        }
        
        # 保存图形
        violin_file <- file.path(violin_dir, paste(sheet, col, "violin_plot_with_stats.pdf", sep = "_"))
        ggsave(violin_file, plot = p1, width = 8, height = 6)
        
        # 密度图
        p4 <- ggplot(merged_data, aes(x = .data[[col]], fill = group)) +
          geom_density(alpha = 0.6) +
          scale_fill_manual(values = tcga_colors) +
          labs(title = paste(sheet, "-", col),
               subtitle = paste("Statistical test:", stat_test$method, "\n", format_pvalue(stat_test$p_value)),
               x = col,
               y = "Density") +
          theme_minimal() +
          theme(
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 18),
            plot.subtitle = element_text(size = 12)
          )
        
        density_file <- file.path(density_dir, paste(sheet, col, "density_plot_with_stats.pdf", sep = "_"))
        ggsave(density_file, plot = p4, width = 8, height = 6)
        
        # 将统计结果保存到文件
        stats_df <- data.frame(
          Sheet = sheet,
          Variable = col,
          Test_Method = stat_test$method,
          P_Value = stat_test$p_value,
          Significance = ifelse(stat_test$p_value < 0.05, "*", "ns")
        )
        
        # 追加写入统计结果
        write.table(
          stats_df,
          file = file.path(output_dir, "statistical_results.csv"),
          append = TRUE,
          sep = ",",
          row.names = FALSE,
          col.names = !file.exists(file.path(output_dir, "statistical_results.csv"))
        )
      }
    }
  }
}

# 修改Neoantigen密度图函数
neoantigen_plot <- function() {
  data <- read_excel("0_Fig3_cohort_data/1_TCGA_LUAD/TCGA附加数据汇总.xlsx", sheet = "Neoantigen")
  
  colnames(data)[colnames(data) == "SampleName"] <- "Sample"
  data <- data[, !colnames(data) %in% c("...1")]
  
  merged_data <- merge(data, meta_data, by = "Sample")
  merged_data$group <- factor(merged_data$group, levels = c("cluster_low", "cluster_medium", "cluster_high"))
  
  # 进行统计检验
  stat_test <- perform_statistical_test(
    merged_data[merged_data$group %in% c("cluster_low", "cluster_high"), ],
    "Neoantigen"
  )
  
  p <- ggplot(merged_data, aes(x = Neoantigen, fill = group)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = tcga_colors) +
    xlim(0, 180) +
    labs(title = "Density Plot of Neoantigen by Group",
         subtitle = paste("Statistical test:", stat_test$method, "\n", format_pvalue(stat_test$p_value)),
         x = "Neoantigen",
         y = "Density") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  ggsave(file.path(output_dir, "Neoantigen_density_plot_180_with_stats.pdf"), 
         plot = p, 
         width = 8, 
         height = 6)
  
  # 保存Neoantigen的统计结果
  stats_df <- data.frame(
    Sheet = "Neoantigen",
    Variable = "Neoantigen",
    Test_Method = stat_test$method,
    P_Value = stat_test$p_value,
    Significance = ifelse(stat_test$p_value < 0.05, "*", "ns")
  )
  
  write.table(
    stats_df,
    file = file.path(output_dir, "statistical_results.csv"),
    append = TRUE,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(file.path(output_dir, "statistical_results.csv"))
  )
}

# 执行Neoantigen图的绘制
neoantigen_plot()




