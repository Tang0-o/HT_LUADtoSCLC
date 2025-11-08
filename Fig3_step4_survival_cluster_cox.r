# ================================================================================
# 完整的多数据集生存分析森林图绘制代码（改进版）
# Complete Multi-Dataset Survival Analysis with Forest Plot (Improved Version)
# ================================================================================

# 加载必需的R包
suppressMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(RColorBrewer)
  library(cowplot)
  library(patchwork)
})

# 设置工作目录和路径
base_dir <- "/home/data/tmh_project/SCLC"
output_dir <- file.path(base_dir, "Fig5_Risk_prediction_model/3_survival_analysis")

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# ================================================================================
# 步骤1: 数据加载和预处理（修复版）
# ================================================================================
cat("=== 步骤1: 数据加载和预处理 ===\n")

load_and_prepare_data <- function() {
  
  cat("开始加载数据...\n")
  combined_data <- NULL
  
  # 检查文件是否存在的辅助函数
  check_file <- function(file_path, description) {
    if (!file.exists(file_path)) {
      cat("警告: 文件不存在 -", description, ":", file_path, "\n")
      return(FALSE)
    } else {
      cat("✓ 找到文件 -", description, "\n")
      return(TRUE)
    }
  }
  
  # 定义所有需要的文件路径
  files_to_check <- list(
    tcga_risk = list(
      path = file.path(base_dir, "Fig5_Risk_prediction_model/1_risk_prediction/TCGA_risk_scores.csv"),
      desc = "TCGA风险评分"
    ),
    oak_risk = list(
      path = file.path(base_dir, "Fig5_Risk_prediction_model/1_risk_prediction/OAK_POPLAR_risk_scores.csv"),
      desc = "OAK风险评分"
    ),
    tcga_surv = list(
      path = file.path(base_dir, "Fig3_multicohort_v1/0_Fig3_cohort_data/1_TCGA_LUAD/tcga_all_LUAD_os_new.csv"),
      desc = "TCGA生存数据"
    ),
    oak_mpdl_surv = list(
      path = file.path(base_dir, "Fig3_multicohort_v1/0_Fig3_cohort_data/5_OAK_POPLAR/oak_poplar_nosquamous_MPDL3280A_os_new.csv"),
      desc = "OAK MPDL3280A生存数据"
    ),
    oak_doce_surv = list(
      path = file.path(base_dir, "Fig3_multicohort_v1/0_Fig3_cohort_data/5_OAK_POPLAR/oak_poplar_nosquamous_Docetaxel_os_new.csv"),
      desc = "OAK Docetaxel生存数据"
    ),
    predictions_oak = list(
      path = file.path(base_dir, "Fig3_multicohort_v1/3_predictions/OAK_POPLAR/ssgsea/hierarchical/knn/predictions_knn.rds"),
      desc = "OAK聚类预测"
    ),
    clustering_output = list(
      path = file.path(base_dir, "Fig3_multicohort_v1/1_TCGA/ssgsea/clustering_output_complete.rds"),
      desc = "TCGA聚类结果"
    )
  )
  
  # 检查所有文件
  all_files_exist <- TRUE
  for (file_info in files_to_check) {
    if (!check_file(file_info$path, file_info$desc)) {
      all_files_exist <- FALSE
    }
  }
  
  if (!all_files_exist) {
    cat("❌ 部分必需文件缺失，无法继续分析\n")
    return(NULL)
  }
  
  cat("\n开始加载各个数据文件...\n")
  
  # 加载数据并进行错误处理
  tryCatch({
    # 加载风险评分数据
    cat("加载风险评分数据...\n")
    tcga_risk_scores <- read.csv(files_to_check$tcga_risk$path, stringsAsFactors = FALSE)
    oak_risk_scores <- read.csv(files_to_check$oak_risk$path, stringsAsFactors = FALSE)
    cat("✓ 风险评分数据加载完成\n")
    
    # 加载聚类数据
    cat("加载聚类数据...\n")
    predictions_OAK <- readRDS(files_to_check$predictions_oak$path)
    clustering_output <- readRDS(files_to_check$clustering_output$path)
    cat("✓ 聚类数据加载完成\n")
    
    # 加载生存数据
    cat("加载生存数据...\n")
    tcga_surv <- read.csv(files_to_check$tcga_surv$path, stringsAsFactors = FALSE)
    oak_mpdl_surv <- read.csv(files_to_check$oak_mpdl_surv$path, stringsAsFactors = FALSE)
    oak_doce_surv <- read.csv(files_to_check$oak_doce_surv$path, stringsAsFactors = FALSE)
    cat("✓ 生存数据加载完成\n")
    
  }, error = function(e) {
    cat("❌ 数据加载错误:", e$message, "\n")
    return(NULL)
  })
  
  # 检查加载的数据
  if (!exists("tcga_risk_scores") || !exists("oak_risk_scores") || 
      !exists("tcga_surv") || !exists("oak_mpdl_surv") || !exists("oak_doce_surv") ||
      !exists("predictions_OAK") || !exists("clustering_output")) {
    cat("❌ 数据加载失败\n")
    return(NULL)
  }
  
  cat("\n数据预处理...\n")
  
  # 标准化样本ID格式
  fmt <- function(x) {
    if (is.null(x) || length(x) == 0) return(character(0))
    return(trimws(gsub("-", "\\.", as.character(x))))
  }
  
  # 处理TCGA数据
  cat("处理TCGA数据...\n")
  
  # 检查必需的列
  if (!"Sample" %in% colnames(tcga_risk_scores)) {
    cat("❌ TCGA风险评分缺少Sample列\n")
    return(NULL)
  }
  
  tcga_risk_scores$sample_id <- fmt(tcga_risk_scores$Sample)
  
  # 检查tcga_surv的列名并映射到标准名称
  tcga_surv_cols <- colnames(tcga_surv)
  cat("TCGA生存数据列名:", paste(tcga_surv_cols, collapse = ", "), "\n")
  
  # 智能识别生存时间和状态列
  sample_col <- NULL
  time_col <- NULL
  status_col <- NULL
  
  # 识别样本ID列
  if ("X" %in% tcga_surv_cols) {
    sample_col <- "X"
  } else if ("sample" %in% tcga_surv_cols) {
    sample_col <- "sample"
  } else if ("SampleID" %in% tcga_surv_cols) {
    sample_col <- "SampleID"
  } else {
    sample_col <- tcga_surv_cols[1]
  }
  
  # 识别生存时间列
  if ("OS.time" %in% tcga_surv_cols) {
    time_col <- "OS.time"
  } else if ("futime" %in% tcga_surv_cols) {
    time_col <- "futime"
  } else if ("time" %in% tcga_surv_cols) {
    time_col <- "time"
  } else if ("days_to_death" %in% tcga_surv_cols) {
    time_col <- "days_to_death"
  }
  
  # 识别生存状态列
  if ("OS" %in% tcga_surv_cols) {
    status_col <- "OS"
  } else if ("fustat" %in% tcga_surv_cols) {
    status_col <- "fustat"
  } else if ("status" %in% tcga_surv_cols) {
    status_col <- "status"
  } else if ("vital_status" %in% tcga_surv_cols) {
    status_col <- "vital_status"
  }
  
  if (is.null(time_col) || is.null(status_col)) {
    cat("❌ 无法识别TCGA生存数据的时间或状态列\n")
    cat("可用列名:", paste(tcga_surv_cols, collapse = ", "), "\n")
    return(NULL)
  }
  
  cat("使用列: 样本ID =", sample_col, ", 时间 =", time_col, ", 状态 =", status_col, "\n")
  
  # 创建标准化的生存数据
  tcga_surv$sample_id <- fmt(tcga_surv[[sample_col]])
  tcga_surv$OS.time <- as.numeric(tcga_surv[[time_col]])
  tcga_surv$OS <- as.numeric(tcga_surv[[status_col]])
  
  # 检查并修正生存状态编码（确保1=死亡事件，0=删失）
  unique_status <- unique(tcga_surv$OS[!is.na(tcga_surv$OS)])
  cat("生存状态唯一值:", paste(unique_status, collapse = ", "), "\n")
  
  # 提取TCGA聚类信息
  if (is.null(clustering_output$cluster_results) || 
      is.null(clustering_output$cluster_results$hierarchical)) {
    cat("❌ TCGA聚类结果格式错误\n")
    return(NULL)
  }
  
  tcga_cluster <- data.frame(
    sample_id = fmt(rownames(clustering_output$cluster_results$hierarchical)),
    cluster = clustering_output$cluster_results$hierarchical$cluster,
    stringsAsFactors = FALSE
  )
  
  # 合并TCGA数据
  tcga_combined <- tcga_surv %>%
    select(sample_id, OS, OS.time) %>%
    inner_join(tcga_risk_scores[, c("sample_id", "Normalized_Risk_Score")], 
               by = "sample_id") %>%
    inner_join(tcga_cluster, by = "sample_id") %>%
    mutate(
      Dataset = "TCGA",
      treatment_arm = "Standard_Care"
    ) %>%
    filter(!is.na(OS), !is.na(OS.time), !is.na(cluster), OS.time > 0)
  
  cat("TCGA数据匹配完成，样本数:", nrow(tcga_combined), "\n")
  
  # 处理OAK数据 - 检查OAK数据的列名
  cat("处理OAK数据...\n")
  
  oak_mpdl_cols <- colnames(oak_mpdl_surv)
  oak_doce_cols <- colnames(oak_doce_surv)
  
  cat("OAK MPDL数据列名:", paste(head(oak_mpdl_cols, 10), collapse = ", "), "\n")
  cat("OAK Docetaxel数据列名:", paste(head(oak_doce_cols, 10), collapse = ", "), "\n")
  
  # 智能识别OAK数据的列名
  identify_oak_columns <- function(df, dataset_name) {
    cols <- colnames(df)
    
    # 样本ID列
    sample_col <- if ("SampleID" %in% cols) "SampleID" else 
      if ("Sample" %in% cols) "Sample" else cols[1]
    
    # 生存时间列
    time_col <- if ("OS.time" %in% cols) "OS.time" else
      if ("futime" %in% cols) "futime" else
        if ("time" %in% cols) "time" else
          if ("survival_time" %in% cols) "survival_time" else NULL
    
    # 生存状态列  
    status_col <- if ("OS" %in% cols) "OS" else
      if ("fustat" %in% cols) "fustat" else
        if ("status" %in% cols) "status" else
          if ("survival_status" %in% cols) "survival_status" else NULL
    
    cat(dataset_name, "- 使用列: 样本ID =", sample_col, ", 时间 =", time_col, ", 状态 =", status_col, "\n")
    
    return(list(sample = sample_col, time = time_col, status = status_col))
  }
  
  # 识别OAK MPDL数据列
  oak_mpdl_cols_map <- identify_oak_columns(oak_mpdl_surv, "OAK MPDL")
  oak_doce_cols_map <- identify_oak_columns(oak_doce_surv, "OAK Docetaxel")
  
  # 处理OAK MPDL3280A数据
  if (is.null(oak_mpdl_cols_map$time) || is.null(oak_mpdl_cols_map$status)) {
    cat("❌ OAK MPDL数据列识别失败\n")
    return(NULL)
  }
  
  oak_mpdl_surv$sample_id <- fmt(oak_mpdl_surv[[oak_mpdl_cols_map$sample]])
  oak_mpdl_surv$OS.time <- as.numeric(oak_mpdl_surv[[oak_mpdl_cols_map$time]])
  oak_mpdl_surv$OS <- as.numeric(oak_mpdl_surv[[oak_mpdl_cols_map$status]])
  
  oak_mpdl_combined <- oak_mpdl_surv %>%
    select(sample_id, OS, OS.time) %>%
    inner_join(oak_risk_scores[, c("Sample", "Normalized_Risk_Score")], 
               by = c("sample_id" = "Sample")) %>%
    mutate(
      Dataset = "OAK_MPDL3280A",
      treatment_arm = "MPDL3280A"
    ) %>%
    filter(!is.na(OS), !is.na(OS.time), OS.time > 0)
  
  # 处理OAK Docetaxel数据
  if (is.null(oak_doce_cols_map$time) || is.null(oak_doce_cols_map$status)) {
    cat("❌ OAK Docetaxel数据列识别失败\n")
    return(NULL)
  }
  
  oak_doce_surv$sample_id <- fmt(oak_doce_surv[[oak_doce_cols_map$sample]])
  oak_doce_surv$OS.time <- as.numeric(oak_doce_surv[[oak_doce_cols_map$time]])
  oak_doce_surv$OS <- as.numeric(oak_doce_surv[[oak_doce_cols_map$status]])
  
  oak_doce_combined <- oak_doce_surv %>%
    select(sample_id, OS, OS.time) %>%
    inner_join(oak_risk_scores[, c("Sample", "Normalized_Risk_Score")], 
               by = c("sample_id" = "Sample")) %>%
    mutate(
      Dataset = "OAK_Docetaxel", 
      treatment_arm = "Docetaxel"
    ) %>%
    filter(!is.na(OS), !is.na(OS.time), OS.time > 0)
  
  # 提取OAK聚类信息
  oak_cluster <- data.frame(
    sample_id = fmt(names(predictions_OAK$class)),
    cluster = predictions_OAK$class,
    stringsAsFactors = FALSE
  )
  
  # 为OAK数据添加聚类信息
  oak_mpdl_combined <- oak_mpdl_combined %>%
    inner_join(oak_cluster, by = "sample_id") %>%
    filter(!is.na(cluster))
  
  oak_doce_combined <- oak_doce_combined %>%
    inner_join(oak_cluster, by = "sample_id") %>%
    filter(!is.na(cluster))
  
  cat("OAK MPDL3280A数据匹配完成，样本数:", nrow(oak_mpdl_combined), "\n")
  cat("OAK Docetaxel数据匹配完成，样本数:", nrow(oak_doce_combined), "\n")
  
  # 合并所有数据
  combined_data <- rbind(
    tcga_combined,
    oak_mpdl_combined, 
    oak_doce_combined
  )
  
  # 标准化聚类名称
  combined_data$cluster <- paste0("cluster_", gsub("cluster_", "", combined_data$cluster))
  
  cat("\n数据加载完成总结:\n")
  cat("  - TCGA样本数:", sum(combined_data$Dataset == "TCGA"), "\n")
  cat("  - OAK MPDL3280A样本数:", sum(combined_data$Dataset == "OAK_MPDL3280A"), "\n")
  cat("  - OAK Docetaxel样本数:", sum(combined_data$Dataset == "OAK_Docetaxel"), "\n")
  cat("  - 总样本数:", nrow(combined_data), "\n")
  cat("  - 聚类分布:\n")
  print(table(combined_data$Dataset, combined_data$cluster))
  
  return(combined_data)
}

# ================================================================================
# 步骤2: Cox回归分析函数
# ================================================================================

perform_comprehensive_cox_analysis <- function(data, min_samples = 15, min_events = 5) {
  
  cat("开始综合Cox回归分析...\n")
  
  cox_results <- list()
  
  # 分析1: 聚类间比较（在每个数据集内部）
  cat("\n--- 聚类间生存差异分析 ---\n")
  
  datasets <- unique(data$Dataset)
  available_clusters <- sort(unique(data$cluster))
  ref_cluster <- "cluster_low"  # 设定参考聚类
  
  if (!ref_cluster %in% available_clusters) {
    ref_cluster <- available_clusters[1]
    cat("参考聚类设为:", ref_cluster, "\n")
  }
  
  for (dataset in datasets) {
    cat("\n处理数据集:", dataset, "\n")
    
    dataset_data <- data[data$Dataset == dataset, ]
    dataset_clusters <- unique(dataset_data$cluster)
    
    for (target_cluster in available_clusters) {
      if (target_cluster == ref_cluster || !target_cluster %in% dataset_clusters) {
        next
      }
      
      # 筛选比较数据
      comparison_data <- dataset_data[dataset_data$cluster %in% c(ref_cluster, target_cluster), ]
      
      # 检查样本量
      group_counts <- table(comparison_data$cluster)
      group_events <- tapply(comparison_data$OS, comparison_data$cluster, sum)
      
      if (any(group_counts < min_samples) || any(group_events < min_events)) {
        cat("  ", target_cluster, "vs", ref_cluster, ": 样本量不足\n")
        next
      }
      
      # 设置因子水平（参考组在前）
      comparison_data$cluster_factor <- factor(comparison_data$cluster, 
                                               levels = c(ref_cluster, target_cluster))
      
      # Cox回归
      tryCatch({
        cox_model <- coxph(Surv(OS.time, OS) ~ cluster_factor, data = comparison_data)
        cox_summary <- summary(cox_model)
        
        if (length(coef(cox_model)) > 0 && !any(is.na(coef(cox_model)))) {
          # 提取结果
          hr <- exp(coef(cox_model)[1])
          ci <- confint(cox_model)
          ci_lower <- exp(ci[1, 1])
          ci_upper <- exp(ci[1, 2])
          p_value <- cox_summary$coefficients[1, "Pr(>|z|)"]
          
          # 保存结果
          result_key <- paste("cluster", dataset, target_cluster, "vs", ref_cluster, sep = "_")
          cox_results[[result_key]] <- list(
            analysis_type = "cluster_comparison",
            dataset = dataset,
            target_cluster = target_cluster,
            reference_cluster = ref_cluster,
            Group_Label = paste0(dataset, ": ", target_cluster, " vs ", ref_cluster),
            HR = hr,
            CI_lower = ci_lower,
            CI_upper = ci_upper,
            P_value = p_value,
            N_total = nrow(comparison_data),
            N_events = sum(comparison_data$OS),
            N_target = group_counts[target_cluster],
            N_reference = group_counts[ref_cluster],
            comparison_type = "cluster"
          )
          
          cat("  ", target_cluster, "vs", ref_cluster, ": HR =", round(hr, 3), 
              ", p =", format.pval(p_value, digits = 3), "\n")
        }
        
      }, error = function(e) {
        cat("  ", target_cluster, "vs", ref_cluster, ": 分析失败 -", e$message, "\n")
      })
    }
  }
  
  # 分析2: 治疗组间比较（在每个聚类内部）
  cat("\n--- 治疗组间生存差异分析 ---\n")
  
  treatment_datasets <- c("OAK_MPDL3280A", "OAK_Docetaxel")
  ref_treatment <- "OAK_MPDL3280A"  # 免疫治疗作为参考
  target_treatment <- "OAK_Docetaxel"  # 化疗作为目标
  
  for (cluster in available_clusters) {
    cat("\n在", cluster, "中比较治疗组:\n")
    
    # 筛选该聚类的治疗数据
    cluster_data <- data[data$cluster == cluster & data$Dataset %in% treatment_datasets, ]
    
    available_treatments <- unique(cluster_data$Dataset)
    if (length(available_treatments) < 2) {
      cat("  治疗组不足\n")
      next
    }
    
    # 检查样本量
    treatment_counts <- table(cluster_data$Dataset)
    treatment_events <- tapply(cluster_data$OS, cluster_data$Dataset, sum)
    
    if (any(treatment_counts < min_samples) || any(treatment_events < min_events)) {
      cat("  样本量不足\n")
      next
    }
    
    # 设置因子水平（参考组在前）
    cluster_data$treatment_factor <- factor(cluster_data$Dataset, 
                                            levels = c(ref_treatment, target_treatment))
    
    # Cox回归
    tryCatch({
      cox_model <- coxph(Surv(OS.time, OS) ~ treatment_factor, data = cluster_data)
      cox_summary <- summary(cox_model)
      
      if (length(coef(cox_model)) > 0 && !any(is.na(coef(cox_model)))) {
        # 提取结果
        hr <- exp(coef(cox_model)[1])
        ci <- confint(cox_model)
        ci_lower <- exp(ci[1, 1])
        ci_upper <- exp(ci[1, 2])
        p_value <- cox_summary$coefficients[1, "Pr(>|z|)"]
        
        # 保存结果
        result_key <- paste("treatment", cluster, target_treatment, "vs", ref_treatment, sep = "_")
        cox_results[[result_key]] <- list(
          analysis_type = "treatment_comparison",
          cluster_context = cluster,
          target_cluster = target_treatment,
          reference_cluster = ref_treatment,
          Group_Label = paste0(cluster, ": Docetaxel vs MPDL3280A"),
          HR = hr,
          CI_lower = ci_lower,
          CI_upper = ci_upper,
          P_value = p_value,
          N_total = nrow(cluster_data),
          N_events = sum(cluster_data$OS),
          N_target = treatment_counts[target_treatment],
          N_reference = treatment_counts[ref_treatment],
          comparison_type = "treatment"
        )
        
        cat("  Docetaxel vs MPDL3280A: HR =", round(hr, 3), 
            ", p =", format.pval(p_value, digits = 3), "\n")
      }
      
    }, error = function(e) {
      cat("  治疗比较失败 -", e$message, "\n")
    })
  }
  
  cat("\nCox分析完成，共生成", length(cox_results), "个比较结果\n")
  return(cox_results)
}

# ================================================================================
# 步骤3: 生存曲线创建函数
# ================================================================================

create_survival_curves <- function(data, output_dir) {
  
  cat("创建生存曲线图...\n")
  
  # 为每个数据集创建聚类生存曲线
  datasets <- unique(data$Dataset)
  
  for (dataset in datasets) {
    dataset_data <- data[data$Dataset == dataset, ]
    
    if (nrow(dataset_data) < 20) {
      cat("跳过", dataset, "- 样本数不足\n")
      next
    }
    
    tryCatch({
      # 创建生存对象
      surv_obj <- Surv(dataset_data$OS.time, dataset_data$OS)
      
      # 拟合生存模型
      surv_fit <- survfit(surv_obj ~ cluster, data = dataset_data)
      
      # 创建生存曲线图
      surv_plot <- ggsurvplot(
        surv_fit,
        data = dataset_data,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        linetype = "strata",
        surv.median.line = "hv",
        ggtheme = theme_bw(),
        palette = c("#E7B800", "#2E9FDF", "#FC4E07"),
        title = paste("Survival Analysis by Cluster -", dataset),
        xlab = "Time (days)",
        ylab = "Survival probability"
      )
      
      # 保存生存曲线
      filename <- paste0("survival_curve_", gsub("[^A-Za-z0-9]", "_", dataset), ".png")
      ggsave(file.path(output_dir, filename), 
             plot = surv_plot$plot, width = 10, height = 8, dpi = 300)
      
      cat("保存生存曲线:", filename, "\n")
      
    }, error = function(e) {
      cat("生存曲线创建失败 -", dataset, ":", e$message, "\n")
    })
  }
  
  # 为OAK数据创建治疗组比较的生存曲线
  oak_data <- data[grepl("OAK_", data$Dataset), ]
  
  if (nrow(oak_data) > 50) {
    clusters <- unique(oak_data$cluster)
    
    for (cluster in clusters) {
      cluster_data <- oak_data[oak_data$cluster == cluster, ]
      
      if (length(unique(cluster_data$Dataset)) < 2) {
        next
      }
      
      tryCatch({
        # 创建生存对象
        surv_obj <- Surv(cluster_data$OS.time, cluster_data$OS)
        
        # 拟合生存模型
        surv_fit <- survfit(surv_obj ~ Dataset, data = cluster_data)
        
        # 创建生存曲线图
        surv_plot <- ggsurvplot(
          surv_fit,
          data = cluster_data,
          pval = TRUE,
          conf.int = TRUE,
          risk.table = TRUE,
          risk.table.col = "strata",
          linetype = "strata",
          surv.median.line = "hv",
          ggtheme = theme_bw(),
          palette = c("#FF7F0E", "#2CA02C"),
          title = paste("Treatment Comparison in", cluster),
          xlab = "Time (days)",
          ylab = "Survival probability"
        )
        
        # 保存生存曲线
        filename <- paste0("treatment_survival_", gsub("[^A-Za-z0-9]", "_", cluster), ".png")
        ggsave(file.path(output_dir, filename), 
               plot = surv_plot$plot, width = 10, height = 8, dpi = 300)
        
        cat("保存治疗组生存曲线:", filename, "\n")
        
      }, error = function(e) {
        cat("治疗组生存曲线创建失败 -", cluster, ":", e$message, "\n")
      })
    }
  }
}

# ================================================================================
# 步骤4: 改进版森林图创建函数
# ================================================================================

create_improved_forest_plot <- function(cox_results, output_dir) {
  
  if (length(cox_results) == 0) {
    cat("❌ 没有Cox回归结果可用于创建森林图\n")
    return(NULL)
  }
  
  cat("创建改进版森林图，包含", length(cox_results), "个比较结果\n")
  
  # 转换结果为数据框
  forest_df <- data.frame()
  
  for (result in cox_results) {
    new_row <- data.frame(
      Dataset = if (!is.null(result$dataset)) result$dataset else "OAK",
      Comparison_Type = result$comparison_type,
      Target_Cluster = result$target_cluster,
      Reference_Cluster = result$reference_cluster,
      Group_Label = result$Group_Label,
      HR = result$HR,
      CI_lower = result$CI_lower,
      CI_upper = result$CI_upper,
      P_value = result$P_value,
      N_total = result$N_total,
      N_events = result$N_events,
      N_target = result$N_target,
      N_reference = result$N_reference,
      Cluster_Context = if (!is.null(result$cluster_context)) result$cluster_context else "",
      stringsAsFactors = FALSE
    )
    forest_df <- rbind(forest_df, new_row)
  }
  
  # 添加格式化列
  forest_df$Significance <- ifelse(forest_df$P_value < 0.001, "***",
                                   ifelse(forest_df$P_value < 0.01, "**", 
                                          ifelse(forest_df$P_value < 0.05, "*", "ns")))
  
  forest_df$HR_95CI <- paste0(sprintf("%.2f", forest_df$HR), 
                              " (", sprintf("%.2f", forest_df$CI_lower), 
                              "-", sprintf("%.2f", forest_df$CI_upper), ")")
  
  forest_df$P_formatted <- ifelse(forest_df$P_value < 0.001, "<0.001", 
                                  sprintf("%.3f", forest_df$P_value))
  
  # 改进的排序逻辑
  cluster_comparisons <- forest_df[forest_df$Comparison_Type == "cluster", ]
  treatment_comparisons <- forest_df[forest_df$Comparison_Type == "treatment", ]
  
  # 为聚类比较添加排序
  if (nrow(cluster_comparisons) > 0) {
    # 数据集排序：TCGA, OAK_MPDL3280A, OAK_Docetaxel
    cluster_comparisons$Dataset_Order <- match(cluster_comparisons$Dataset, 
                                               c("TCGA", "OAK_MPDL3280A", "OAK_Docetaxel"))
    # 聚类排序：high, medium
    cluster_comparisons$Cluster_Order <- match(cluster_comparisons$Target_Cluster, 
                                               c("cluster_high", "cluster_medium"))
    cluster_comparisons <- cluster_comparisons[order(cluster_comparisons$Cluster_Order, 
                                                     cluster_comparisons$Dataset_Order, na.last = TRUE), ]
  }
  
  # 为治疗比较添加排序
  if (nrow(treatment_comparisons) > 0) {
    treatment_comparisons$Cluster_Order <- match(treatment_comparisons$Cluster_Context, 
                                                 c("cluster_high", "cluster_medium", "cluster_low"))
    treatment_comparisons <- treatment_comparisons[order(treatment_comparisons$Cluster_Order, na.last = TRUE), ]
  }
  
  # 合并数据
  all_columns <- unique(c(colnames(cluster_comparisons), colnames(treatment_comparisons)))
  
  for (col in all_columns) {
    if (!col %in% colnames(cluster_comparisons) && nrow(cluster_comparisons) > 0) {
      cluster_comparisons[[col]] <- NA
    }
    if (!col %in% colnames(treatment_comparisons) && nrow(treatment_comparisons) > 0) {
      treatment_comparisons[[col]] <- NA
    }
  }
  
  if (nrow(cluster_comparisons) > 0 && nrow(treatment_comparisons) > 0) {
    cluster_comparisons <- cluster_comparisons[, all_columns]
    treatment_comparisons <- treatment_comparisons[, all_columns]
    final_forest_df <- rbind(cluster_comparisons, treatment_comparisons)
  } else if (nrow(cluster_comparisons) > 0) {
    final_forest_df <- cluster_comparisons
  } else {
    final_forest_df <- treatment_comparisons
  }
  
  final_forest_df$Plot_Order <- nrow(final_forest_df):1
  rownames(final_forest_df) <- NULL
  
  # 改进的颜色和形状方案
  final_forest_df$Color_Group <- NA
  final_forest_df$Shape_Group <- final_forest_df$Comparison_Type
  
  # 聚类比较的颜色：按数据集分组
  cluster_idx <- final_forest_df$Comparison_Type == "cluster"
  if (sum(cluster_idx) > 0) {
    final_forest_df$Color_Group[cluster_idx] <- paste0("Cluster_", final_forest_df$Dataset[cluster_idx])
  }
  
  # 治疗比较的颜色：按聚类背景分组
  treatment_idx <- final_forest_df$Comparison_Type == "treatment"
  if (sum(treatment_idx) > 0) {
    final_forest_df$Color_Group[treatment_idx] <- paste0("Treatment_", final_forest_df$Cluster_Context[treatment_idx])
  }
  
  # 定义改进的配色方案
  color_palette <- c(
    # 聚类比较 - 按数据集
    "Cluster_TCGA" = "#1f77b4",              # 蓝色 - TCGA
    "Cluster_OAK_MPDL3280A" = "#ff7f0e",     # 橙色 - OAK免疫治疗
    "Cluster_OAK_Docetaxel" = "#2ca02c",     # 绿色 - OAK化疗
    
    # 治疗比较 - 按聚类背景
    "Treatment_cluster_high" = "#d62728",     # 红色 - 高风险聚类中的治疗比较
    "Treatment_cluster_medium" = "#9467bd",   # 紫色 - 中风险聚类中的治疗比较
    "Treatment_cluster_low" = "#8c564b"       # 棕色 - 低风险聚类中的治疗比较
  )
  
  # 形状方案
  shape_palette <- c("cluster" = 16, "treatment" = 17)  # 圆形 vs 三角形
  
  cat("森林图数据处理完成，颜色分组:", paste(unique(final_forest_df$Color_Group), collapse = ", "), "\n")
  
  # 设置显示范围
  final_forest_df$HR_plot <- pmax(pmin(final_forest_df$HR, 5), 0.2)
  final_forest_df$CI_lower_plot <- pmax(final_forest_df$CI_lower, 0.2)
  final_forest_df$CI_upper_plot <- pmin(final_forest_df$CI_upper, 5)
  
  # 计算治疗组比较的数量用于背景
  n_treatment <- sum(final_forest_df$Comparison_Type == "treatment", na.rm = TRUE)
  
  # 创建改进版森林图
  library(ggplot2)
  
  improved_forest_plot <- ggplot(final_forest_df, aes(x = HR_plot, y = Plot_Order)) +
    
    # 背景区域（仅当有治疗组比较时）
    {if (n_treatment > 0) {
      annotate("rect", xmin = 0.2, xmax = 5, 
               ymin = 0.5, 
               ymax = n_treatment + 0.5,
               fill = "#f0f0f0", alpha = 0.3)
    }} +
    
    # 参考线 HR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "#d62728", 
               linewidth = 1, alpha = 0.8) +
    
    # 置信区间
    geom_errorbarh(aes(xmin = CI_lower_plot, xmax = CI_upper_plot, 
                       color = Color_Group), 
                   height = 0.15, linewidth = 1, alpha = 0.8) +
    
    # 点估计
    geom_point(aes(color = Color_Group, shape = Shape_Group), 
               size = 4, alpha = 0.9, stroke = 1.2) +
    
    # 应用配色和形状
    scale_color_manual(values = color_palette, 
                       name = "Data Group",
                       labels = function(x) gsub("_", " ", gsub("Cluster_|Treatment_", "", x))) +
    scale_shape_manual(values = shape_palette, 
                       name = "Comparison",
                       labels = c("cluster" = "Risk Groups", "treatment" = "Treatments")) +
    
    # X轴（对数尺度，范围优化）
    scale_x_log10(
      limits = c(0.2, 5),
      breaks = c(0.2, 0.5, 1, 2, 5),
      labels = c("0.2", "0.5", "1.0", "2.0", "5.0")
    ) +
    
    # Y轴
    scale_y_continuous(
      breaks = final_forest_df$Plot_Order,
      labels = final_forest_df$Group_Label,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    # 标题和标签
    labs(
      title = "Multi-Dataset Survival Analysis Forest Plot",
      subtitle = paste("Risk Stratification and Treatment Effects |", 
                       sum(final_forest_df$Significance != "ns", na.rm = TRUE), "of", 
                       nrow(final_forest_df), "comparisons significant"),
      x = "Hazard Ratio (95% CI)",
      y = "",
      caption = paste("Reference: cluster_low for risk comparisons, MPDL3280A for treatment comparisons |", Sys.Date())
    ) +
    
    # 精致主题
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = 9, hjust = 0.5, color = "gray50", 
                                  margin = margin(t = 8)),
      axis.text.y = element_text(size = 10, hjust = 1),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.margin = margin(t = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
      plot.margin = margin(15, 120, 15, 15),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    
    # 数值标注（位置优化）
    geom_text(aes(x = 5.2, y = Plot_Order, 
                  label = paste0(HR_95CI, "\np=", P_formatted, " ", Significance)), 
              hjust = 0, vjust = 0.5, size = 2.8, color = "black") +
    
    # 显著性标记
    geom_text(aes(x = 0.18, y = Plot_Order,
                  label = ifelse(Significance != "ns", "●", "")),
              color = "#d62728", size = 3, hjust = 0.5, vjust = 0.5)
  
  # 添加分组标签（优化位置）
  if (n_treatment > 0 && sum(final_forest_df$Comparison_Type == "cluster", na.rm = TRUE) > 0) {
    n_cluster <- sum(final_forest_df$Comparison_Type == "cluster", na.rm = TRUE)
    improved_forest_plot <- improved_forest_plot +
      annotate("text", x = 0.16, y = n_treatment + n_cluster/2, 
               label = "Risk Group\nComparisons", 
               angle = 90, vjust = 0.5, hjust = 0.5, 
               size = 3.5, fontface = "bold", color = "#1f77b4") +
      annotate("text", x = 0.16, y = n_treatment/2, 
               label = "Treatment\nComparisons", 
               angle = 90, vjust = 0.5, hjust = 0.5, 
               size = 3.5, fontface = "bold", color = "#d62728")
  }
  
  # 保存改进版森林图（精致尺寸）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 动态计算图片高度，但保持精致
  plot_height <- max(6, min(12, nrow(final_forest_df) * 0.8))
  plot_width <- 12
  
  # 保存PDF
  ggsave(file.path(output_dir, "improved_forest_plot.pdf"), 
         improved_forest_plot, width = plot_width, height = plot_height, device = "pdf")
  
  # 保存PNG（高分辨率）
  ggsave(file.path(output_dir, "improved_forest_plot.png"), 
         improved_forest_plot, width = plot_width, height = plot_height, dpi = 300, device = "png")
  
  cat("✅ 改进版森林图已创建并保存!\n")
  cat("   尺寸:", plot_width, "x", plot_height, "英寸\n")
  cat("   PDF: improved_forest_plot.pdf\n")
  cat("   PNG: improved_forest_plot.png\n")
  
  # 保存改进结果数据
  write.csv(final_forest_df, file.path(output_dir, "improved_forest_results.csv"), 
            row.names = FALSE)
  cat("✅ 改进结果数据已保存: improved_forest_results.csv\n")
  
  # 创建简化的结果表格用于展示
  summary_table <- final_forest_df[, c("Group_Label", "Comparison_Type", "HR", "CI_lower", "CI_upper", 
                                       "P_value", "Significance", "N_total", "N_events")]
  colnames(summary_table) <- c("Comparison", "Type", "HR", "CI_Lower", "CI_Upper", 
                               "P_Value", "Significance", "N_Total", "N_Events")
  
  summary_table$HR <- round(summary_table$HR, 3)
  summary_table$CI_Lower <- round(summary_table$CI_Lower, 3)
  summary_table$CI_Upper <- round(summary_table$CI_Upper, 3)
  summary_table$P_Value <- ifelse(summary_table$P_Value < 0.001, "<0.001", 
                                  sprintf("%.3f", summary_table$P_Value))
  
  write.csv(summary_table, file.path(output_dir, "improved_cox_analysis_summary.csv"), 
            row.names = FALSE)
  cat("✅ 简化摘要表格已保存: improved_cox_analysis_summary.csv\n")
  
  return(list(
    plot = improved_forest_plot,
    data = final_forest_df,
    summary = summary_table
  ))
}

# ================================================================================
# 主分析函数
# ================================================================================

main_analysis <- function() {
  
  cat("\n🚀 开始多数据集森林图分析...\n")
  cat("===========================================\n")
  
  analysis_results <- list()
  
  # 步骤1: 数据加载和预处理
  cat("\n执行步骤1: 数据加载和预处理\n")
  combined_data <- load_and_prepare_data()
  
  if (is.null(combined_data)) {
    cat("❌ 数据加载失败，分析终止\n")
    return(NULL)
  }
  
  analysis_results$data <- combined_data
  
  # 步骤2: Cox回归分析
  cat("\n执行步骤2: 综合Cox回归分析\n")
  cox_results <- perform_comprehensive_cox_analysis(combined_data)
  
  if (length(cox_results) == 0) {
    cat("❌ Cox分析未产生结果，分析终止\n")
    return(NULL)
  }
  
  analysis_results$cox_results <- cox_results
  
  # 步骤3: 创建改进森林图
  cat("\n执行步骤3: 创建改进版森林图\n")
  forest_results <- create_improved_forest_plot(cox_results, output_dir)
  
  if (!is.null(forest_results)) {
    analysis_results$forest_plot <- forest_results$plot
    analysis_results$forest_data <- forest_results$data
    analysis_results$forest_summary <- forest_results$summary
  }
  
  # 步骤4: 创建生存曲线
  cat("\n执行步骤4: 创建生存曲线\n")
  create_survival_curves(combined_data, output_dir)
  
  # 最终结果摘要
  cat("\n=== 🎊 分析结果摘要 ===\n")
  
  if (!is.null(forest_results)) {
    # 按类型统计结果
    cluster_results <- forest_results$data[forest_results$data$Comparison_Type == "cluster", ]
    treatment_results <- forest_results$data[forest_results$data$Comparison_Type == "treatment", ]
    
    # 聚类比较结果
    if (nrow(cluster_results) > 0) {
      cat("\n📊 聚类间比较结果:\n")
      for (i in 1:nrow(cluster_results)) {
        row <- cluster_results[i, ]
        significance_symbol <- if (row$Significance != "ns") "🔥" else "  "
        
        cat(paste0(significance_symbol, " ", row$Group_Label, "\n"))
        cat(paste0("     HR: ", row$HR_95CI, ", p = ", row$P_formatted, " ", row$Significance, "\n"))
        cat(paste0("     样本: ", row$N_total, " (target: ", row$N_target, 
                   " vs ref: ", row$N_reference, "), 事件: ", row$N_events, "\n\n"))
      }
    }
    
    # 治疗组比较结果
    if (nrow(treatment_results) > 0) {
      cat("💊 治疗组间比较结果 (Docetaxel vs MPDL3280A):\n")
      for (i in 1:nrow(treatment_results)) {
        row <- treatment_results[i, ]
        significance_symbol <- if (row$Significance != "ns") "🔥" else "  "
        
        cat(paste0(significance_symbol, " ", row$Group_Label, "\n"))
        cat(paste0("     HR: ", row$HR_95CI, ", p = ", row$P_formatted, " ", row$Significance, "\n"))
        cat(paste0("     样本: ", row$N_total, " (Docetaxel: ", row$N_target, 
                   " vs MPDL3280A: ", row$N_reference, "), 事件: ", row$N_events, "\n\n"))
      }
    }
    
    # 总体统计
    total_significant <- sum(forest_results$data$Significance != "ns", na.rm = TRUE)
    cluster_significant <- if (nrow(cluster_results) > 0) sum(cluster_results$Significance != "ns", na.rm = TRUE) else 0
    treatment_significant <- if (nrow(treatment_results) > 0) sum(treatment_results$Significance != "ns", na.rm = TRUE) else 0
    
    cat("📈 总体统计:\n")
    cat("- 总比较数:", nrow(forest_results$data), "\n")
    cat("- 总显著结果:", total_significant, "/", nrow(forest_results$data), "\n")
    cat("- 聚类比较显著:", cluster_significant, "/", nrow(cluster_results), "\n")
    cat("- 治疗比较显著:", treatment_significant, "/", nrow(treatment_results), "\n")
  }
  
  cat("\n🎊 分析完全完成!\n")
  
  return(analysis_results)
}

# ================================================================================
# 执行分析
# ================================================================================

cat("开始执行完整分析...\n")
results <- main_analysis()

# 如果分析成功，显示最终摘要
if (!is.null(results)) {
  cat("\n📊 最终分析摘要:\n")
  cat("- 总样本数:", nrow(results$data), "\n")
  cat("- Cox比较数:", length(results$cox_results), "\n")
  
  if (!is.null(results$forest_data)) {
    significant_results <- sum(results$forest_data$Significance != "ns", na.rm = TRUE)
    cat("- 显著结果数:", significant_results, "\n")
  }
  
  cat("\n✅ 所有结果已保存至:", output_dir, "\n")
  
  # 列出生成的文件
  cat("\n📁 生成的文件:\n")
  if (dir.exists(output_dir)) {
    files <- list.files(output_dir, pattern = "\\.(pdf|png|csv)$")
    for (file in files) {
      cat("  -", file, "\n")
    }
  }
  
  # 显示数据集分布
  cat("\n📈 数据集分布:\n")
  dataset_summary <- results$data %>%
    group_by(Dataset, cluster) %>%
    summarise(
      n_samples = n(),
      n_events = sum(OS),
      event_rate = round(mean(OS) * 100, 1),
      .groups = 'drop'
    ) %>%
    arrange(Dataset, cluster)
  
  print(dataset_summary)
  
} else {
  cat("❌ 分析失败，请检查数据和路径设置\n")
  cat("可能的问题:\n")
  cat("1. 检查文件路径是否正确\n")
  cat("2. 检查数据文件是否存在\n") 
  cat("3. 检查数据文件格式是否正确\n")
  cat("4. 检查R包是否正确安装\n")
}

cat("\n🎊 多数据集森林图分析程序执行完毕！\n")

