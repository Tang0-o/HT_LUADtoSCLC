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
    "trans_prob_low" = "#00468B",
    "trans_prob_median" =  "#ADB6CA",
    "trans_prob_high" = "#FDAF91",
    # 新增通用命名配色
    "cluster_low" = "#00468B",
    "cluster_medium" = "#ADB6CA",
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
# 2. 预测函数 (优化版)
#################################################################################
predict_trans_cluster <- function(exp_data, 
                                  gmt_path,
                                  method = c("rf", "svm", "knn", "pca", "graph"),
                                  models_path = "02_models",
                                  clustering_method = c("hierarchical", "kmeans", "pam"),
                                  scale_data = FALSE,
                                  gsva_method = "ssgsea",
                                  output_path = "02_models/OAK") {
  
  # 参数检查
  method <- match.arg(method)
  clustering_method <- match.arg(clustering_method)
  
  if (is.null(gsva_method) || gsva_method == "") {
    stop("`gsva_method` must be specified to locate the correct model directory.")
  }
  
  # 根据GSVA方法和聚类方法构建模型路径
  current_models_path <- file.path(models_path, gsva_method, clustering_method)
  
  # 创建输出目录
  if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  # 检查输入数据
  if(!is.matrix(exp_data) && !is.data.frame(exp_data)) {
    stop("Input data must be a matrix or data frame")
  }
  exp_data <- na.omit(exp_data) ##删除含NA的行
  
  # 首先计算regulon得分
  message(sprintf("Step 1: Calculating regulon scores... [Clustering method: %s]", clustering_method))
  # 确保表达矩阵格式正确
  if(max(exp_data) > 100) {
    message("Performing log2 transformation...")
    exp_data <- log2(exp_data + 1)
  }
  
  # 读取regulon基因集
  tryCatch({
    geneset <- getGmt(gmt_path)
  }, error = function(e) {
    stop("无法读取GMT文件: ", gmt_path, "\n错误信息: ", e$message)
  })
  
  # 运行ssGSEA
  tryCatch({
    ssgsea_par <- ssgseaParam(as.matrix(exp_data), geneset)
    new_data <- gsva(ssgsea_par)
  }, error = function(e) {
    stop("ssGSEA计算失败: ", e$message)
  })
  
  # 如果需要标准化
  if(scale_data) {
    new_data <- t(scale(t(new_data)))
  }
  
  message(sprintf("Step 2: Performing cluster prediction... [Method: %s, Clustering: %s]", 
                  method, clustering_method))
  
  # 保存原始样本名
  original_samples <- colnames(new_data)
  if(is.null(original_samples)) {
    original_samples <- paste0("Sample_", 1:ncol(new_data))
    colnames(new_data) <- original_samples
  }
  
  # 确保输入数据格式正确（样本为行）
  new_data <- as.data.frame(t(new_data))
  
  # 检查模型目录
  if(!dir.exists(current_models_path)) {
    stop(sprintf("Models directory does not exist: %s\nPlease run train_and_save_models() first with clustering_method = '%s'",
                 current_models_path, clustering_method))
  }
  
  tryCatch({
    predictions <- switch(method,
                          "rf" = {
                            # Random Forest prediction
                            model_file <- file.path(current_models_path, "rf_model.rds")
                            if(!file.exists(model_file)) {
                              stop("Random Forest model file not found: ", model_file)
                            }
                            model <- readRDS(model_file)
                            list(
                              class = predict(model, new_data),
                              probabilities = predict(model, new_data, type = "prob")
                            )
                          },
                          "svm" = {
                            # SVM prediction
                            model_file <- file.path(current_models_path, "svm_model.rds")
                            if(!file.exists(model_file)) {
                              stop("SVM model file not found: ", model_file)
                            }
                            model <- readRDS(model_file)
                            list(
                              class = predict(model, new_data),
                              probabilities = attr(predict(model, new_data, probability = TRUE), 
                                                   "probabilities")
                            )
                          },
                          "knn" = {
                            # KNN prediction
                            model_file <- file.path(current_models_path, "pca_model.rds")
                            if(!file.exists(model_file)) {
                              stop("PCA model file not found: ", model_file)
                            }
                            pca_model <- readRDS(model_file)
                            
                            # 数据标准化
                            train_scaled <- scale(pca_model$training_data)
                            test_scaled <- scale(new_data)
                            
                            # KNN预测
                            pred <- class::knn(train = train_scaled,
                                               test = test_scaled,
                                               cl = factor(pca_model$training_labels),
                                               k = 5,
                                               prob = TRUE)
                            
                            # 创建概率矩阵
                            unique_labels <- unique(pca_model$training_labels)
                            prob_matrix <- matrix(0, nrow = nrow(new_data), 
                                                  ncol = length(unique_labels))
                            colnames(prob_matrix) <- unique_labels
                            
                            # 设置预测概率
                            for(i in seq_along(pred)) {
                              prob_matrix[i, as.character(pred[i])] <- attr(pred, "prob")[i]
                            }
                            
                            list(
                              class = pred,
                              probabilities = prob_matrix
                            )
                          },
                          "pca" = {
                            # PCA-based prediction
                            model_file <- file.path(current_models_path, "pca_model.rds")
                            if(!file.exists(model_file)) {
                              stop("PCA model file not found: ", model_file)
                            }
                            pca_model <- readRDS(model_file)
                            
                            # 确保列名匹配
                            if(!all(colnames(pca_model$training_data) %in% colnames(new_data))) {
                              colnames(new_data) <- colnames(pca_model$training_data)
                            }
                            
                            test_pca <- predict(pca_model$pca, newdata = new_data)[, 1:2]
                            
                            distances <- sapply(1:nrow(pca_model$centers), function(i) {
                              apply(test_pca, 1, function(x) 
                                dist(rbind(x, as.numeric(pca_model$centers[i, 2:3]))))
                            })
                            
                            prob_matrix <- t(apply(-distances, 1, function(x) exp(x)/sum(exp(x))))
                            colnames(prob_matrix) <- pca_model$centers$cluster
                            
                            list(
                              class = pca_model$centers$cluster[max.col(-distances)],
                              probabilities = prob_matrix
                            )
                          },
                          "graph" = {
                            # Graph-based prediction
                            model_file <- file.path(current_models_path, "pca_model.rds")
                            if(!file.exists(model_file)) {
                              stop("PCA model file not found: ", model_file)
                            }
                            pca_model <- readRDS(model_file)
                            
                            # 使用PCA方法作为fallback
                            if(nrow(new_data) < 5) {
                              warning("样本数量过少，使用PCA方法替代")
                              return(predict_trans_cluster(exp_data, gmt_path, "pca", models_path, scale_data, gsva_method, output_path))
                            }
                            
                            # 确保列名匹配
                            if(!all(colnames(pca_model$training_data) %in% colnames(new_data))) {
                              colnames(new_data) <- colnames(pca_model$training_data)
                            }
                            
                            # PCA转换
                            test_pca <- predict(pca_model$pca, newdata = new_data)[, 1:2]
                            
                            # 计算距离和概率
                            distances <- sapply(1:nrow(pca_model$centers), function(i) {
                              apply(test_pca, 1, function(x) 
                                dist(rbind(x, as.numeric(pca_model$centers[i, 2:3]))))
                            })
                            
                            prob_matrix <- t(apply(-distances, 1, function(x) exp(x)/sum(exp(x))))
                            colnames(prob_matrix) <- pca_model$centers$cluster
                            
                            list(
                              class = pca_model$centers$cluster[max.col(-distances)],
                              probabilities = prob_matrix
                            )
                          }
    )
    
    # 添加样本名到结果中
    names(predictions$class) <- original_samples
    rownames(predictions$probabilities) <- original_samples
    
    # 添加预测方法信息
    attr(predictions, "method") <- method
    
    # 添加预测摘要
    predictions$summary <- data.frame(
      Sample = original_samples,
      Predicted_Class = predictions$class,
      row.names = original_samples
    )
    
    # 添加regulon得分到结果中
    predictions$regulon_scores <- new_data
    
    # 添加原始表达矩阵到结果中
    predictions$expression_matrix <- exp_data
    
    # 添加聚类方法信息到预测结果中
    attr(predictions, "clustering_method") <- clustering_method
    
    # 添加GSVA方法信息到预测结果中
    attr(predictions, "gsva_method") <- gsva_method
    
    # 保存预测结果
    saveRDS(predictions, file.path(output_path, paste0("predictions_", method, ".rds")))
    
    message(sprintf("✓ 预测完成 [聚类方法: %s, 预测方法: %s]", clustering_method, method))
    message(sprintf("✓ 结果已保存到: %s", output_path))
    return(predictions)
    
  }, error = function(e) {
    warning(paste("预测错误:", e$message))
    if(method != "pca") {
      warning("切换到PCA方法")
      result <- predict_trans_cluster(exp_data, gmt_path, "pca", models_path, scale_data, gsva_method, output_path)
      return(result)
    } else {
      stop("PCA方法也失败了: ", e$message)
    }
  })
}

#################################################################################
# 3. 绘制预测结果热图 (使用Lancet配色)
#################################################################################
plot_prediction_heatmap <- function(predictions, output_path = "02_models/OAK") {
  # 获取预测方法和聚类方法名称
  method_name <- attr(predictions, "method")
  clustering_method <- attr(predictions, "clustering_method")
  
  # 参数检查
  if(!is.list(predictions) || is.null(predictions$regulon_scores)) {
    stop("Invalid predictions object. Must be output from predict_trans_cluster()")
  }
  
  # 创建输出目录
  if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  # 确保Predicted_Class是有序因子，并按指定顺序排列
  predictions$class <- factor(
    predictions$class, 
    levels = c("cluster_low", "cluster_medium", "cluster_high"),
    ordered = TRUE
  )
  
  # 根据排序后的聚类水平获取排序索引
  sort_index <- order(predictions$class)
  
  # 准备热图数据 - 转置数据使得Regulon为行，并按cluster排序样本
  data_h <- t(as.matrix(predictions$regulon_scores[sort_index, ]))
  
  # 准备注释数据 - 同样按照排序后的顺序
  annotation_df <- data.frame(
    Predicted_Class = predictions$class[sort_index],
    row.names = names(predictions$class)[sort_index]
  )
  
  # 定义颜色映射，确保顺序和水平一致
  color_map <- c(
    "cluster_low" = "#00468B", 
    "cluster_medium" =  "#ADB6CA",
    "cluster_high" = "#FDAF91"
  )
  
  # 确保注释颜色与因子水平完全匹配
  ann_colors <- list(
    Predicted_Class = setNames(
      color_map[levels(predictions$class)], 
      levels(predictions$class)
    )
  )
  
  # 打印配色信息用于调试
  message("热图聚类配色映射:")
  for(cluster_name in levels(predictions$class)) {
    color_val <- ann_colors$Predicted_Class[cluster_name]
    message("  ", cluster_name, " -> ", color_val)
  }
  
  # 创建热图文件名，包含方法名称和聚类方法
  heatmap_filename <- paste0("prediction_heatmap_", clustering_method, "_", method_name, ".pdf")

  # 创建热图
  pdf(file.path(output_path, heatmap_filename), width = 7, height = 7)
  
  # p1 <- pheatmap(
  #   data_h,
  #   annotation_col = annotation_df,
  #   annotation_colors = ann_colors,
  #   color = LANCET_COLORS$heatmap_colors,
  #   cluster_cols = FALSE,
  #   cluster_rows = TRUE,
  #   scale = "row",
  #   show_colnames = FALSE,
  #   show_rownames = TRUE,
  #   border_color = NA,
  #   fontsize = 10,
  #   fontsize_row = 10,
  #   fontsize_col = 8,
  #   legend = TRUE,
  #   main = paste("Regulon Activity Pattern -", toupper(clustering_method), "-", toupper(method_name)),
  #   annotation_names_col = FALSE,
  #   cellwidth = NA,
  #   cellheight = 15,
  #   treeheight_row = 30,
  #   treeheight_col = 0,
  #   annotation_legend = TRUE,
  #   name = "Activity"
  # )
  data_h_scaled <- t(scale(t(data_h)))   # 行 z-score
  breaks <- seq(-3, 3, length.out = 101) # 下限-3，上限3，共101格
  
  p1 <- pheatmap(
    data_h,
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    color = LANCET_COLORS$heatmap_colors,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    scale = "row",
    show_colnames = FALSE,
    show_rownames = TRUE,
    border_color = NA,
    fontsize = 12,
    fontsize_row = max(10, min(14, 14 - floor(nrow(data_h)/15))),  # 保持行标签字体大小
    fontsize_col = 12,  # 保持列标签字体大小
    legend = TRUE,
    main = paste("Regulon Activity Pattern -", toupper(clustering_method), "-", toupper(method_name)),
    annotation_names_col = FALSE,
    cellwidth = NA,
    cellheight = 15,
    treeheight_row = 30,
    treeheight_col = 0,
    annotation_legend = TRUE,
    name = "Activity",
    breaks = seq(-3, 3, length.out = 101)  # Set range from -4 to 4 with 101 breaks
  )
  
  
  dev.off()
  
  message("✓ 热图已保存到: ", file.path(output_path, heatmap_filename))
  
  return(p1)
}

#################################################################################
# 4. 绘制KM生存曲线 (使用Lancet配色)
#################################################################################
plot_km_curves <- function(predictions,
                           surv_data,
                           time_col = "time",
                           status_col = "status",
                           sample_col = "sample",
                           title = "Kaplan-Meier Survival Curves",
                           output_path = "02_models/OAK") {
  # 获取预测方法和聚类方法名称
  method_name <- attr(predictions, "method")
  clustering_method <- attr(predictions, "clustering_method")
  
  # 创建输出目录
  if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  # 准备数据
  pred_df <- data.frame(
    sample = predictions$summary$Sample,
    cluster = predictions$summary$Predicted_Class,
    stringsAsFactors = TRUE
  )
  
  # 合并数据
  combined_data <- merge(surv_data, pred_df, by.x = sample_col, by.y = "sample", all = FALSE)
  
  if(nrow(combined_data) == 0) {
    stop("生存数据与预测结果之间没有匹配的样本")
  }
  
  # 确保数据完整性并创建所需的向量
  complete_cases <- complete.cases(combined_data[c(time_col, status_col, "cluster")])
  combined_data <- combined_data[complete_cases, ]
  
  # 提取需要的列并创建新的数据框
  surv_df <- data.frame(
    time = combined_data[[time_col]],
    status = combined_data[[status_col]],
    cluster = combined_data$cluster
  )
  
  # 创建生存对象
  fit <- survfit(Surv(time, status) ~ cluster, data = surv_df)
  
  # 计算两两比较的p值
  ps <- pairwise_survdiff(Surv(time, status) ~ cluster, data = surv_df)
  
  # 创建p值表格
  addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                           round(ps$p.value, 3))))
  addTab[is.na(addTab)] <- "-"
  
  # 定义固定的颜色映射
  color_map <- c(
    "cluster_low" = "#00468B",      # 深蓝色
    "cluster_medium" = "#ADB6CA",   # 浅蓝灰色
    "cluster_high" = "#FDAF91"      # 粉橙色
  )
  
  # 获取实际的聚类水平
  cluster_levels <- levels(factor(surv_df$cluster))
  
  # 创建颜色向量，确保与实际聚类水平对应
  custom_colors <- color_map[cluster_levels]
  
  # 如果某些聚类没有对应的颜色，使用默认颜色
  missing_colors <- is.na(custom_colors)
  if(any(missing_colors)) {
    custom_colors[missing_colors] <- LANCET_COLORS$survival_colors[1:sum(missing_colors)]
  }
  
  # 创建生存曲线文件名
  survival_filename <- paste0("survival_curves_", clustering_method, "_", method_name, ".pdf")
  
  # 绘制生存曲线
  km_plot <- ggsurvplot(
    fit,
    data = surv_df,
    pval = FALSE,  # 移除总体p值显示
    conf.int = FALSE,
    risk.table = TRUE,
    palette = custom_colors,
    title = paste(title, "-", toupper(clustering_method), "-", toupper(method_name)),
    xlab = "Time (months)",
    ylab = "Survival Probability",
    legend.title = "Predicted Cluster",
    legend.labs = cluster_levels,
    risk.table.height = 0.25,
    tables.theme = theme_cleantable(),
    ncensor.plot = FALSE,
    surv.plot.height = 0.7,
    break.time.by = 12,
    ggtheme = theme_classic(),
    font.main = c(16, "bold", "black"),
    font.x = c(14, "plain", "black"),
    font.y = c(14, "plain", "black"),
    font.legend = c(12, "plain", "black"),
    font.tickslab = c(12, "plain", "black"),
    tables.col = "strata",
    risk.table.fontsize = 4,
    pval.method = FALSE  # 不显示p值计算方法
  )
  
  # 添加Log-rank标签和两两比较的p值表格
  km_plot$plot <- km_plot$plot +
    annotate("text", x = 24, y = 0.25, label = "Log-rank", hjust = 0, size = 4)
  
  # 添加两两比较的p值表格
  if(!is.null(ps$p.value) && any(!is.na(ps$p.value))) {
    x_pos <- max(surv_df$time) * 0.2  # 调整位置
    df <- tibble(x = x_pos,
                 y = 0.1,
                 tb = list(addTab))
    km_plot$plot <- km_plot$plot +
      ggpp::geom_table(data = df,
                       aes(x = x, y = y, label = tb),
                       table.rownames = TRUE,
                       size = 3.5,  # 减小表格大小
                       hjust = 0,
                       vjust = 0)
  }
  
  # 保存图片，减小尺寸
  pdf(file.path(output_path, survival_filename), width = 6, height = 6)  # 减小尺寸
  print(km_plot)
  dev.off()
  
  message("✓ 生存曲线图已保存到: ", file.path(output_path, survival_filename))
  
  return(km_plot)
}



#################################################################################
# 5. 模型评估 (使用Lancet配色)
#################################################################################
evaluate_predictions_simple <- function(predictions, 
                                        output_dir = "02_models/OAK/evaluation/", 
                                        prefix = "model_eval") {
  
  # 获取预测方法和聚类方法名称
  method_name <- attr(predictions, "method")
  clustering_method <- attr(predictions, "clustering_method")
  
  # 创建输出目录
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 初始化评估结果列表
  evaluation_results <- list()
  
  # 1. 概率分布分析
  prob_melted <- melt(predictions$probabilities)
  colnames(prob_melted) <- c("Sample", "Cluster", "Probability")
  
  # 创建概率分布图文件名，包含方法名称和聚类方法
  prob_dist_filename <- paste0(prefix, "_", clustering_method, "_", method_name, "_probability_distribution.pdf")
  
  pdf(file.path(output_dir, prob_dist_filename), 
      width = 10, height = 8)
  
  p1 <- ggplot(prob_melted, aes(x = Probability, fill = Cluster)) +
    geom_density(alpha = 0.7) +
    theme_minimal() +
    scale_fill_manual(values = LANCET_COLORS$extended_colors) +
    labs(title = paste("Prediction Probability Distribution -", toupper(clustering_method), "-", toupper(method_name)),
         x = "Prediction Probability",
         y = "Density") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 12))
  print(p1)
  dev.off()
  
  # 2. 计算基本评估指标
  metrics <- list(
    n_samples = nrow(predictions$probabilities),
    n_clusters = ncol(predictions$probabilities),
    avg_confidence = mean(apply(predictions$probabilities, 1, max)),
    min_confidence = min(apply(predictions$probabilities, 1, max)),
    cluster_sizes = table(predictions$class)
  )
  
  # 3. 生成评估报告
  report <- sprintf('
========================================
肿瘤转移概率预测模型评估报告 - %s方法
========================================
生成时间: %s
预测方法: %s

1. 基本信息:
   - 总样本数: %d
   - 聚类数目: %d
   - 预测方法: %s

2. 预测置信度分析:
   a) 平均置信度: %.3f
      说明: 模型预测的平均置信水平
      - 范围: 0 到 1
      - 评估标准:
        * > 0.9: 高置信度
        * 0.7-0.9: 良好置信度
        * 0.5-0.7: 中等置信度
        * < 0.5: 低置信度
      当前评估: %s

   b) 最低置信度: %.3f
      说明: 最不确定预测的置信水平
      - 标准:
        * > 0.5: 可接受
        * < 0.5: 需要注意
      当前评估: %s

3. 聚类分布:
%s
   说明: 显示各聚类的样本分布情况
   标准:
   - 理想: 均匀分布 (1:%d比例)
   - 警告: 单个聚类 > 50%%
   当前评估: %s

4. 总体评估结论:
', 
                    # 时间和方法信息
                    toupper(method_name),
                    Sys.time(),
                    attr(predictions, "method"),
                    
                    # 基本信息
                    metrics$n_samples,
                    metrics$n_clusters,
                    attr(predictions, "method"),
                    
                    # 置信度评估
                    metrics$avg_confidence,
                    if(metrics$avg_confidence > 0.9) "高置信度预测"
                    else if(metrics$avg_confidence > 0.7) "良好置信度水平"
                    else if(metrics$avg_confidence > 0.5) "中等置信度水平"
                    else "低置信度水平",
                    
                    metrics$min_confidence,
                    if(metrics$min_confidence > 0.5) "所有预测都有可接受的置信度"
                    else "部分预测置信度较低",
                    
                    # 聚类分布
                    paste(capture.output(metrics$cluster_sizes), collapse = "\n"),
                    metrics$n_clusters,
                    if(max(prop.table(metrics$cluster_sizes)) < 0.3) "分布均衡"
                    else if(max(prop.table(metrics$cluster_sizes)) < 0.5) "分布可接受"
                    else "分布不均衡"
  )
  
  # 4. 生成建议
  suggestions <- character()
  
  # 基于预测置信度的建议
  if(metrics$avg_confidence < 0.7) {
    suggestions <- c(suggestions, "⚠ 警告: 预测置信度较低，建议:
      * 检查模型参数设置
      * 考虑增加训练数据量
      * 优化特征选择策略")
  }
  
  # 基于聚类平衡性的建议
  cluster_props <- prop.table(metrics$cluster_sizes)
  if(max(cluster_props) > 0.5) {
    suggestions <- c(suggestions, "⚠ 警告: 检测到聚类不平衡，建议:
      * 考虑调整聚类数目
      * 检查数据分布情况  
      * 优化聚类参数设置")
  }
  
  if(length(suggestions) == 0) {
    suggestions <- c("✓ 模型整体表现良好:
      * 预测置信度高
      * 聚类分布合理
    建议:
      * 继续使用当前模型
      * 定期评估模型性能")
  }
  
  report <- paste(report, paste(suggestions, collapse = "\n\n"), 
                  "\n========================================", sep = "\n")
  
  # 创建评估报告文件名，包含方法名称和聚类方法
  report_filename <- paste0(prefix, "_", clustering_method, "_", method_name, "_evaluation_report.txt")
  data_filename <- paste0(prefix, "_", clustering_method, "_", method_name, "_evaluation_data.rds")
  
  # 5. 保存评估结果
  writeLines(report, file.path(output_dir, report_filename))
  evaluation_results$report <- report
  evaluation_results$metrics <- metrics
  
  # 保存评估数据
  saveRDS(evaluation_results, file.path(output_dir, data_filename))
  
  # 在控制台输出关键建议
  message("\n========================================")
  message(sprintf("模型评估关键建议 [%s - %s]", toupper(clustering_method), toupper(method_name)))
  message("========================================")
  message(paste(suggestions, collapse = "\n"))
  message("详细报告已保存到: ", file.path(output_dir, report_filename))
  
  return(evaluation_results)
}

