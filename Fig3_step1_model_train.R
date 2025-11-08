#################################################################################
# 肿瘤转移概率预测模型 - 模型构建及预测部分
# 统一使用Lancet配色方案
#################################################################################

# Load required packages
required_packages <- c(
  "ggplot2", "reshape2", "scales", "cluster", "fpc", 
  "boot", "caret", "clValid", "viridis", "uwot", 
  "randomForest", "e1071", "pheatmap", "RColorBrewer",
  "survival", "survminer", "dplyr", "tibble", "ggpp",
  "GSVA", "GSEABase", "class", "RANN", "igraph"
)

# 改进包的加载部分
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#################################################################################
# 全局配色方案设置 (Lancet配色)
#################################################################################
LANCET_COLORS <- list(
  # 三类别聚类配色（保留原有命名兼容性）
  cluster_colors = c(
    "trans_prob_low" = "#ADB6CA",
    "trans_prob_median" = "#00468B", 
    "trans_prob_high" = "#FDAF91",
    # 新增通用命名配色
    "cluster_low" = "#ADB6CA",
    "cluster_medium" = "#00468B",
    "cluster_high" = "#FDAF91",
    "cluster_1" = "#00468B",
    "cluster_2" = "#FDAF91", 
    "cluster_3" = "#ADB6CA",
    "cluster_4" = "#2E8B57",
    "cluster_5" = "#DC143C"
  ),
  # 生存曲线配色
  survival_colors = c("#FDAF91", "#ADB6CA","#00468B", "#2E8B57", "#DC143C"),
  # 热图配色
  heatmap_colors = colorRampPalette(c("#2E8B57", "#FFFFFF", "#DC143C"))(100),
  # 扩展配色
  extended_colors = c("#00468B", "#FDAF91", "#ADB6CA", "#2E8B57", "#DC143C", 
                      "#42B883", "#E17C05", "#6495ED", "#FF6347", "#32CD32")
)

#################################################################################
# 1. 训练和保存模型 (优化版)
#################################################################################
#' Train and save models for multiple clustering methods
#' @param clustering_results 聚类结果列表，包含merged_data和cluster_results
#' @param save_path 模型保存路径，默认为"02_models"
#' @param clustering_methods 要处理的聚类方法，默认全部处理
#' @param gsva_method GSVA分析方法，可以是"ssgsea", "gsva", "zscore", "plage"之一
#' 
train_and_save_models <- function(
    clustering_results, 
    save_path = "02_models",
    clustering_methods = c("hierarchical", "kmeans", "pam"),
    gsva_method = "ssgsea"  # 默认使用ssgsea方法
) {
  # 验证输入参数
  if(!is.list(clustering_results) || 
     is.null(clustering_results$merged_data) || 
     is.null(clustering_results$cluster_results)) {
    stop("无效的聚类结果对象")
  }
  
  # 验证输入的聚类方法
  available_methods <- names(clustering_results$merged_data)
  invalid_methods <- setdiff(clustering_methods, available_methods)
  if(length(invalid_methods) > 0) {
    stop("无效的聚类方法: ", paste(invalid_methods, collapse = ", "))
  }
  
  # 验证GSVA方法
  valid_gsva_methods <- c("ssgsea", "gsva", "zscore", "plage")
  if(!gsva_method %in% valid_gsva_methods) {
    stop("无效的GSVA方法: ", gsva_method, "。必须是以下之一: ", paste(valid_gsva_methods, collapse = ", "))
  }
  
  # 创建保存目录，包含GSVA方法
  save_path <- file.path(save_path, gsva_method)
  if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  # 为每种聚类方法训练模型
  all_models <- list()
  
  for(method in clustering_methods) {
    message("\n========================================")
    message(sprintf("处理聚类方法: %s [GSVA方法: %s]", method, gsva_method))
    message("========================================")
    
    # 创建方法特定的保存目录
    method_path <- file.path(save_path, method)
    if(!dir.exists(method_path)) dir.create(method_path, recursive = TRUE)
    
    # 准备训练数据
    merged_data <- clustering_results$merged_data[[method]]
    
    # 确保行名一致性
    common_samples <- intersect(
      rownames(merged_data),
      rownames(clustering_results$cluster_results[[method]])
    )
    
    if(length(common_samples) == 0) {
      warning("方法 ", method, " 没有共同样本，跳过")
      next
    }
    
    # 只选择数值列，并排除cluster列
    X_train <- merged_data[common_samples, 1:min(13, ncol(merged_data))]
    
    # 获取标签
    y_train <- clustering_results$cluster_results[[method]][common_samples, "cluster"]
    
    # 移除包含NA或Inf的行
    valid_rows <- complete.cases(X_train) & 
      !apply(X_train, 1, function(x) any(is.infinite(x)))
    X_train <- X_train[valid_rows, ]
    y_train <- y_train[valid_rows]
    
    # 打印数据维度
    message("清理后的训练数据维度: ", paste(dim(X_train), collapse = " x "))
    message("清理后的标签数量: ", length(y_train))
    message("类别分布: ", paste(names(table(y_train)), table(y_train), sep = "=", collapse = ", "))
    
    if(nrow(X_train) < 10 || length(unique(y_train)) < 2) {
      warning("数据量不足，跳过方法: ", method)
      next
    }
    
    # 初始化当前方法的模型列表
    method_models <- list()
    
    # 随机森林模型
    message("\n训练随机森林模型...")
    rf_model <- try({
      randomForest(x = X_train, 
                   y = factor(y_train), 
                   ntree = 500,
                   importance = TRUE)
    }, silent = TRUE)
    
    if(!inherits(rf_model, "try-error")) {
      saveRDS(rf_model, file.path(method_path, "rf_model.rds"))
      method_models$rf <- rf_model
      message("✓ 随机森林模型保存成功")
    } else {
      message("✗ 随机森林训练失败: ", rf_model[1])
    }
    
    # SVM模型
    message("训练SVM模型...")
    svm_model <- try({
      svm(x = as.matrix(X_train), 
          y = factor(y_train), 
          probability = TRUE, 
          kernel = "radial")
    }, silent = TRUE)
    
    if(!inherits(svm_model, "try-error")) {
      saveRDS(svm_model, file.path(method_path, "svm_model.rds"))
      method_models$svm <- svm_model
      message("✓ SVM模型保存成功")
    } else {
      message("✗ SVM训练失败: ", svm_model[1])
    }
    
    # PCA模型
    message("计算PCA并保存训练数据...")
    pca_model <- try({
      # 标准化数据
      scaled_data <- scale(X_train)
      
      # 执行PCA
      pca_result <- prcomp(scaled_data, scale. = FALSE)
      
      # 计算聚类中心
      pca_data <- data.frame(pca_result$x[, 1:min(2, ncol(pca_result$x))])
      centers <- aggregate(pca_data, 
                           by = list(cluster = y_train), 
                           FUN = mean)
      
      # 创建完整的PCA模型对象
      list(
        pca = pca_result,
        centers = centers,
        training_data = X_train,
        training_labels = y_train,
        scaling_params = list(
          center = attr(scaled_data, "scaled:center"),
          scale = attr(scaled_data, "scaled:scale")
        )
      )
    }, silent = TRUE)
    
    if(!inherits(pca_model, "try-error")) {
      saveRDS(pca_model, file.path(method_path, "pca_model.rds"))
      method_models$pca <- pca_model
      message("✓ PCA模型保存成功")
    } else {
      message("✗ PCA计算失败: ", pca_model[1])
    }
    
    # 保存所有模型的性能指标
    performance_metrics <- list()
    
    if(!inherits(rf_model, "try-error")) {
      rf_pred <- predict(rf_model, X_train)
      performance_metrics$rf <- list(
        accuracy = mean(rf_pred == y_train),
        confusion_matrix = table(Predicted = rf_pred, Actual = y_train)
      )
    }
    
    if(!inherits(svm_model, "try-error")) {
      svm_pred <- predict(svm_model, X_train)
      performance_metrics$svm <- list(
        accuracy = mean(svm_pred == y_train),
        confusion_matrix = table(Predicted = svm_pred, Actual = y_train)
      )
    }
    
    # 保存性能指标
    saveRDS(performance_metrics, 
            file.path(method_path, "performance_metrics.rds"))
    
    # 将当前方法的模型添加到总列表中
    all_models[[method]] <- list(
      models = method_models,
      performance = performance_metrics,
      data_info = list(
        n_samples = nrow(X_train),
        n_features = ncol(X_train),
        n_clusters = length(unique(y_train))
      )
    )
    
    message("✓ ", method, " 聚类方法的模型训练完成")
  }
  
  # 保存汇总信息
  summary_info <- lapply(all_models, function(x) {
    list(
      available_models = names(x$models),
      performance_summary = lapply(x$performance, function(p) p$accuracy),
      data_info = x$data_info
    )
  })
  
  # 在文件名中包含GSVA方法
  saveRDS(summary_info, file.path(save_path, paste0("model_summary_", gsva_method, ".rds")))
  
  message("\n========================================")
  message(sprintf("所有模型训练完成！[GSVA方法: %s]", gsva_method))
  message("保存位置: ", save_path)
  message("========================================")
  
  return(all_models)
}