#' SCLC数据分析包 - LANCET期刊风格版本
#' 包含两个主要函数：
#' 1. perform_regulon_clustering: 进行regulon分析和聚类
#' 2. perform_survival_analysis: 进行生存分析

# ====================================================================
# LANCET期刊配色方案
# ====================================================================

LANCET_COLORS <- list(
  # 三类别聚类配色（保留原有命名兼容性）
  cluster_colors = c(
    "trans_prob_low" =  "#00468B",
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
  survival_colors = c("#00468B", "#FDAF91", "#ADB6CA", "#2E8B57", "#DC143C"),
  # 热图配色
  heatmap_colors = colorRampPalette(c("#2E8B57", "#FFFFFF", "#DC143C"))(100),
  # 扩展配色
  extended_colors = c("#00468B", "#FDAF91", "#ADB6CA", "#2E8B57", "#DC143C", 
                      "#42B883", "#E17C05", "#6495ED", "#FF6347", "#32CD32")
)

# Helper function for sample ID normalization
normalize_sample_ids <- function(ids) {
  if (is.null(ids) || length(ids) == 0) return(ids)
  
  ids_normalized <- as.character(ids)
  
  # Replace hyphens and underscores with periods
  ids_normalized <- gsub("[_\\-]", ".", ids_normalized)
  
  # Convert to UPPERCASE
  ids_normalized <- toupper(ids_normalized)
  
  # Standardize TCGA-like IDs: aim for TCGA.XX.YYYY.NN or TCGA.XX.YYYY.NNA
  for (i in seq_along(ids_normalized)) {
    if (startsWith(ids_normalized[i], "TCGA.")) {
      # Try to match up to the sample vial letter (e.g., TCGA.XX.YYYY.01A) - 16 chars
      match_16_char <- regexpr("^TCGA\\.[A-Z0-9]{2}\\.[A-Z0-9]{4}\\.\\d{2}[A-Z]", ids_normalized[i])
      if (match_16_char != -1 && attr(match_16_char, "match.length") > 0) {
        ids_normalized[i] <- substr(ids_normalized[i], 1, attr(match_16_char, "match.length"))
        next # Move to next ID
      }
      
      # Try to match up to the sample type (e.g., TCGA.XX.YYYY.01) - 15 chars
      match_15_char <- regexpr("^TCGA\\.[A-Z0-9]{2}\\.[A-Z0-9]{4}\\.\\d{2}", ids_normalized[i])
      if (match_15_char != -1 && attr(match_15_char, "match.length") > 0) {
        ids_normalized[i] <- substr(ids_normalized[i], 1, attr(match_15_char, "match.length"))
        next # Move to next ID
      }
    }
  }
  
  return(ids_normalized)
}

# GSVA分析函数
run_gsva_analysis <- function(exp_file, 
                              geneset_file,
                              output_prefix = "gsva_result",
                              filter_quantile = 0.1,
                              log_transform = TRUE,
                              log_threshold = 50,
                              ssgsea_alpha = 0.25,
                              method = "ssgsea",
                              kcdf_method = "Gaussian",
                              verbose = TRUE) {
  
  # 加载必要的包
  if(!requireNamespace("GSVA", quietly = TRUE)) install.packages("GSVA")
  library(GSVA)
  if(!requireNamespace("GSEABase", quietly = TRUE)) install.packages("GSEABase")
  library(GSEABase)
  
  # 验证方法参数
  valid_methods <- c("ssgsea", "gsva", "zscore", "plage")
  method <- tolower(method)
  if(!method %in% valid_methods) {
    warning(paste("无效的方法:", method, ". 将回退到ssgsea。"))
    method <- "ssgsea"
  }
  
  if(verbose) cat("=== 开始GSVA分析 ===\n")
  if(verbose) cat("使用方法:", method, "\n")
  
  # 读取表达数据
  if(verbose) cat("读取表达数据:", exp_file, "\n")
  exp_data <- read.csv(exp_file, row.names = 1, header = TRUE, sep = ",", quote = "\"")
  
  if(verbose && ncol(exp_data) > 0) {
    cat("读取后表达数据样本ID (前5个):", paste(head(colnames(exp_data), 5), collapse=", "), "\n")
  }
  colnames(exp_data) <- normalize_sample_ids(colnames(exp_data))
  if(verbose && ncol(exp_data) > 0) {
    cat("标准化后表达数据样本ID (前5个):", paste(head(colnames(exp_data), 5), collapse=", "), "\n")
  }
  
  # 读取基因集文件
  if(verbose) cat("读取基因集文件:", geneset_file, "\n")
  geneset <- GSEABase::getGmt(geneset_file)
  
  if(verbose) {
    cat("表达数据维度:", paste(dim(exp_data), collapse="x"), "\n")
    cat("基因集数量:", length(geneset), "\n")
  }
  
  # 数据预处理
  if(verbose) cat("\n=== 数据预处理 ===\n")
  
  exp_matrix <- as.matrix(exp_data)
  if(verbose) {
    cat("矩阵维度:", paste(dim(exp_matrix), collapse="x"), "\n")
    cat("数据类型:", class(exp_matrix), "\n")
  }
  
  # 数据清理
  if(verbose) cat("\n=== 数据清理 ===\n")
  
  na_count <- sum(is.na(exp_matrix))
  nan_count <- sum(is.nan(exp_matrix))
  inf_count <- sum(is.infinite(exp_matrix))
  
  if(verbose) {
    cat("NA值数量:", na_count, "\n")
    cat("NaN值数量:", nan_count, "\n")
    cat("无穷大值数量:", inf_count, "\n")
  }
  
  # 处理异常值
  if(na_count > 0) {
    exp_matrix[is.na(exp_matrix)] <- 0
    if(verbose) cat("已将", na_count, "个NA值替换为0\n")
  }
  if(nan_count > 0) {
    exp_matrix[is.nan(exp_matrix)] <- 0
    if(verbose) cat("已将", nan_count, "个NaN值替换为0\n")
  }
  if(inf_count > 0) {
    exp_matrix[is.infinite(exp_matrix)] <- 0
    if(verbose) cat("已将", inf_count, "个无穷大值替换为0\n")
  }
  
  # 数据转换
  if(verbose) cat("\n=== 数据转换 ===\n")
  
  if(verbose && length(exp_matrix) > 0) cat("原始数据范围: [", min(exp_matrix, na.rm = TRUE), ",", max(exp_matrix, na.rm = TRUE), "]\n") 
  
  # 根据参数决定是否进行log转换
  if(length(exp_matrix) > 0) {
    max_val <- max(exp_matrix, na.rm = TRUE)
    if(log_transform && max_val > log_threshold) {
      if(verbose) cat("应用log2转换（阈值:", log_threshold, "）\n")
      exp_matrix_final <- log2(exp_matrix + 1)
      if(verbose) cat("log2转换后数据范围: [", min(exp_matrix_final, na.rm = TRUE), ",", max(exp_matrix_final, na.rm = TRUE), "]\n") 
    } else {
      exp_matrix_final <- exp_matrix
      if(verbose) cat("使用原始数据或未达到log转换阈值\n")
    }
  } else {
    exp_matrix_final <- exp_matrix
    if(verbose) cat("表达矩阵为空，跳过log转换\n")
  }
  
  # 基因过滤
  if(verbose) cat("\n=== 基因过滤 ===\n")
  
  if(nrow(exp_matrix_final) > 0 && ncol(exp_matrix_final) > 0) {
    gene_means <- rowMeans(exp_matrix_final, na.rm = TRUE)
    
    valid_genes <- !is.na(gene_means) & !is.nan(gene_means) & !is.infinite(gene_means)
    gene_means_clean <- gene_means[valid_genes]
    exp_matrix_clean <- exp_matrix_final[valid_genes, , drop = FALSE] 
    
    if(verbose) {
      cat("清理前基因数:", length(gene_means), "\n")
      cat("清理后基因数:", length(gene_means_clean), "\n")
    }
    
    if(filter_quantile > 0 && filter_quantile < 1 && length(gene_means_clean) > 0) {
      min_expr_threshold <- quantile(gene_means_clean, filter_quantile, na.rm = TRUE)
      keep_genes <- gene_means_clean >= min_expr_threshold
      exp_matrix_filtered <- exp_matrix_clean[keep_genes, , drop = FALSE] 
      if(verbose) {
        cat("过滤阈值（保留表达量 >=", filter_quantile*100, "%分位数）:", min_expr_threshold, "\n")
        cat("最终保留基因数:", nrow(exp_matrix_filtered), "\n")
      }
    } else {
      exp_matrix_filtered <- exp_matrix_clean
      if(verbose) cat("跳过基因表达量过滤或无有效基因进行过滤\n")
    }
  } else {
    exp_matrix_filtered <- exp_matrix_final
    if(verbose) cat("表达矩阵过滤步骤跳过 (矩阵为空或维度不足)\n")
  }
  
  # 基因重叠检查
  if(verbose) cat("\n=== 基因重叠检查 ===\n")
  
  geneset_list <- NULL
  if (is(geneset, "GeneSetCollection")) {
    geneset_list <- geneIds(geneset)
  } else if (is.list(geneset)) {
    geneset_list <- geneset 
  } else {
    warning("基因集格式无法识别。期望是GeneSetCollection或列表。跳过GSVA。")
    return(list(enrichment_result = NULL, method = method, processed_matrix = exp_matrix_filtered, overlap_genes = character(0), geneset = NULL))
  }
  
  geneset_genes <- unique(unlist(geneset_list))
  overlap_genes <- character(0)
  if(nrow(exp_matrix_filtered) > 0) {
    overlap_genes <- intersect(rownames(exp_matrix_filtered), geneset_genes)
  }
  
  if(verbose) {
    cat("基因集中总基因数:", length(geneset_genes), "\n")
    cat("表达矩阵中用于分析的基因数:", nrow(exp_matrix_filtered), "\n")
    cat("重叠基因数:", length(overlap_genes), "\n")
    if (length(geneset_genes) > 0 && nrow(exp_matrix_filtered) > 0) { 
      cat("重叠比例 (基于表达矩阵基因):", round(length(overlap_genes)/length(intersect(rownames(exp_matrix_filtered),unique(rownames(exp_matrix_filtered))))*100, 2), "%\n")
      cat("重叠比例 (基于基因集基因):", round(length(overlap_genes)/length(geneset_genes)*100, 2), "%\n")
    } else {
      cat("重叠比例: N/A (基因集或过滤后表达矩阵为空)\n")
    }
    
    if(length(overlap_genes) == 0 && length(geneset_genes) > 0 && nrow(exp_matrix_filtered) > 0) {
      cat("⚠️ 警告：重叠基因数为0。请检查基因ID匹配（例如，ENSEMBL vs Symbol）和物种。\n")
    } else if (length(geneset_genes) > 0 && nrow(exp_matrix_filtered) > 0 && length(overlap_genes) < length(geneset_genes) * 0.1) {
      cat("⚠️ 警告：重叠基因比例非常低 (<10%)，可能严重影响分析结果\n")
    } else if (length(geneset_genes) == 0) {
      cat("⚠️ 警告：基因集为空或未能正确加载。\n")
    } else if (nrow(exp_matrix_filtered) == 0) {
      cat("⚠️ 警告：过滤后表达矩阵不含任何基因。\n")
    }
    else if (length(overlap_genes) > 0) { 
      cat("✅ 基因重叠比例良好\n")
    }
  }
  
  # 执行富集分析
  if(verbose) cat(paste("\n=== 执行", toupper(method), "分析 ===\n"))
  
  enrichment_result <- NULL
  if(nrow(exp_matrix_filtered) > 0 && length(geneset_list) > 0 && length(overlap_genes) > 0) {
    tryCatch({
      valid_geneSets <- Filter(function(x) length(x) > 0 && is.character(x), geneset_list)
      if(length(valid_geneSets) == 0) {
        stop("没有有效的基因集（所有基因集为空或非字符型）。")
      }
      if(length(valid_geneSets) < length(geneset_list)) {
        warning(paste(length(geneset_list) - length(valid_geneSets), "个基因集为空或格式无效，已被移除。"))
      }
      
      # 根据方法选择不同的参数
      gsva_params <- switch(method,
                            "ssgsea" = ssgseaParam(
                              exprData = exp_matrix_filtered,
                              geneSets = valid_geneSets, 
                              alpha = ssgsea_alpha,
                              normalize = TRUE
                            ),
                            "gsva" = gsvaParam(
                              exprData = exp_matrix_filtered,
                              geneSets = valid_geneSets,
                              kcdf = kcdf_method
                            ),
                            "zscore" = zscoreParam(
                              exprData = exp_matrix_filtered,
                              geneSets = valid_geneSets
                            ),
                            "plage" = plageParam(
                              exprData = exp_matrix_filtered,
                              geneSets = valid_geneSets
                            )
      )
      
      if(verbose) cat(paste("开始", toupper(method), "计算...\n"))
      enrichment_result <- gsva(gsva_params, verbose = verbose)
      
      if(verbose) {
        cat(paste("✅", toupper(method), "分析成功完成！\n"))
        cat("结果矩阵维度:", paste(dim(enrichment_result), collapse="x"), "\n")
      }
      
    }, error = function(e) {
      if(verbose) cat(paste("❌", toupper(method), "分析失败:", e$message, "\n"))
    })
  } else {
    if(verbose) cat(paste("❌ 跳过", toupper(method), "分析：过滤后表达矩阵为空，或基因集列表为空/无效，或无重叠基因。\n"))
  }
  
  # 保存结果
  if(verbose) cat("\n=== 保存结果 ===\n")
  
  if(!is.null(enrichment_result) && nrow(enrichment_result) > 0 && ncol(enrichment_result) > 0) {
    result_file <- paste0(output_prefix, "_", method, ".csv")
    tryCatch(write.csv(enrichment_result, result_file, row.names = TRUE), 
             error = function(e) warning(paste("无法保存", toupper(method), "结果到", result_file,":", e$message)))
    if(verbose) cat(paste(toupper(method), "结果已尝试保存到:", result_file, "\n"))
  } else {
    if(verbose) cat(paste(toupper(method), "结果为空或未生成，不保存。\n"))
  }
  
  # 分析总结
  if(verbose) {
    cat("\n=== 分析总结 ===\n")
    cat("原始表达数据:", nrow(exp_data), "基因 ×", ncol(exp_data), "样本\n")
    cat("用于", toupper(method), "的矩阵:", nrow(exp_matrix_filtered), "基因 ×", ncol(exp_matrix_filtered), "样本\n")
    cat("基因集:", length(geneset_list), "个 (有效用于分析的可能更少)\n")
    cat("重叠基因:", length(overlap_genes), "个\n")
    if(!is.null(enrichment_result) && nrow(enrichment_result) > 0) {
      cat(paste("✅", toupper(method), "分析成功\n"))
    } else {
      cat(paste("❌", toupper(method), "分析未成功或无结果\n"))
    }
    
    cat("🎉 分析完成！\n")
  }
  
  # 内存清理
  gc()
  
  # 返回主要结果
  return(list(
    enrichment_result = enrichment_result,
    method = method,
    processed_matrix = exp_matrix_filtered,
    overlap_genes = overlap_genes,
    geneset = geneset_list 
  ))
} 

#' Regulon分析和聚类函数
perform_regulon_clustering <- function(
    exp_data_path,                    
    gmt_file_path,                    
    output_dir,                       
    clustering_methods = c("hierarchical"), 
    n_clusters = 3,                   
    gsva_method = "ssgsea",
    verbose_clustering = TRUE
) {
  # 验证聚类方法
  valid_methods <- c("hierarchical", "kmeans", "pam")
  invalid_methods <- setdiff(clustering_methods, valid_methods)
  if(length(invalid_methods) > 0) {
    stop(paste("无效的聚类方法:", paste(invalid_methods, collapse = ", ")))
  }
  
  # 加载必要的包 
  required_packages <- c(
    "limma", "data.table", "ggplot2", "dplyr", "tibble",
    "pheatmap", "cluster" 
  )
  
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 运行富集分析
  if(verbose_clustering) print(paste("正在使用 run_gsva_analysis 进行富集分析 (请求的方法:", gsva_method, ")"))
  
  analysis_results <- run_gsva_analysis(
    exp_file = exp_data_path,
    geneset_file = gmt_file_path,
    output_prefix = file.path(output_dir, "gsva_run_output"), 
    method = tolower(gsva_method),
    verbose = verbose_clustering
  )
  
  gsva_data <- NULL
  effective_gsva_method <- analysis_results$method
  
  # 获取富集分析结果
  gsva_data <- analysis_results$enrichment_result
  if(is.null(gsva_data) || nrow(gsva_data) == 0 || ncol(gsva_data) == 0) {
    stop(paste("方法 '", gsva_method, "' 失败或未产生有效结果。无法继续。", sep=""))
  } else {
    if(verbose_clustering) print(paste("使用", toupper(effective_gsva_method), "结果进行后续分析。"))
  }
  
  if(is.null(gsva_data) || nrow(gsva_data) == 0 || ncol(gsva_data) == 0){
    stop("GSVA/ssGSEA富集分数矩阵为空或无效，无法进行聚类。")
  }
  
  # 准备聚类数据
  rt <- as.matrix(gsva_data)
  if(!is.numeric(rt)) {
    warning("GSVA富集分数数据在转换为矩阵后非数值型，尝试强制转换。")
    rt_dimnames <- dimnames(rt)
    rt_temp <- matrix(as.numeric(rt), nrow=nrow(rt), ncol=ncol(rt))
    if(any(is.na(rt_temp)) && !any(is.na(rt))) warning("强制转换为数值型时引入了NA值。")
    rt <- rt_temp
    dimnames(rt) <- rt_dimnames
  }
  data <- rt 
  
  if(ncol(data) < n_clusters && n_clusters > 1) { 
    stop(paste("样本数量 (", ncol(data), ") 小于请求的聚类数 (", n_clusters, ")，无法进行聚类。"))
  }
  if(ncol(data) == 0) stop("GSVA/ssGSEA结果不包含任何样本 (列)，无法聚类。")
  if(nrow(data) == 0) stop("GSVA/ssGSEA结果不包含任何基因集 (行)，无法聚类。")
  
  # 执行聚类和绘图
  all_results <- list()
  all_merged_data <- list()
  
  for(method in clustering_methods) {
    if(verbose_clustering) print(paste("正在执行", method, "聚类..."))
    
    cluster_results_df <- NULL
    if (n_clusters == 1) { 
      if(verbose_clustering) print("n_clusters is 1, assigning all samples to a single cluster.")
      cluster_labels <- rep(1, ncol(data))
      cluster_results_df <- data.frame(cluster = factor(cluster_labels), row.names = colnames(data))
    } else { 
      cluster_results_df <- switch(
        method,
        "hierarchical" = {
          if(ncol(data) < 2) stop("层次聚类至少需要2个样本。")
          dist_matrix <- dist(t(data)) 
          if(any(!is.finite(dist_matrix))) stop("距离矩阵中存在NA/NaN/Inf。检查GSVA分数是否存在恒定值或全零样本。")
          hc <- hclust(dist_matrix, method = "ward.D2")
          setNames(as.data.frame(cutree(hc, k = n_clusters)), 'cluster')
        },
        "kmeans" = {
          if(ncol(data) < n_clusters) stop(paste("K-means：样本数 (", ncol(data), ") 小于聚类数 (", n_clusters, ").", sep=""))
          set.seed(123) 
          km_data <- t(data)
          if(nrow(km_data) <= n_clusters) stop("K-means: number of samples must be greater than number of clusters.")
          col_vars <- apply(km_data, 2, var, na.rm = TRUE)
          if(any(col_vars == 0, na.rm = TRUE)) {
            warning("Kmeans：部分regulon在样本间无变异，可能影响聚类或对结果无贡献。")
          }
          km <- kmeans(km_data, centers = n_clusters, nstart = 25, iter.max = 50) 
          setNames(as.data.frame(km$cluster), 'cluster')
        },
        "pam" = {
          if(ncol(data) < n_clusters) stop(paste("PAM：样本数 (", ncol(data), ") 小于聚类数 (", n_clusters, ").", sep=""))
          pam_data <- t(data)
          if(nrow(pam_data) <= n_clusters) stop("PAM: number of samples must be greater than number of clusters.")
          pam_result <- cluster::pam(pam_data, k = n_clusters)
          setNames(as.data.frame(pam_result$clustering), 'cluster')
        },
        stop(paste("未知的聚类方法:", method))
      )
    } 
    
    if(is.null(cluster_results_df)) {
      warning(paste("聚类方法 '", method, "' 未能生成结果。跳过后续步骤。"))
      next 
    }
    rownames(cluster_results_df) <- colnames(data) 
    
    # 保存聚类结果
    cluster_file <- file.path(output_dir, paste0("clustering_results_", method, "_", effective_gsva_method, ".csv"))
    tryCatch({
      write.csv(cluster_results_df, cluster_file, row.names = TRUE)
      if(verbose_clustering) cat("📁 聚类结果已保存:", cluster_file, "\n")
    }, error = function(e) {
      warning(paste("保存聚类结果失败:", e$message))
    })
    
    gene_regulon <- as.data.frame(t(data)) 
    cluster_means <- calculate_cluster_means(gene_regulon, cluster_results_df, verbose = verbose_clustering)
    cluster_results_df_typed <- assign_cluster_types(cluster_results_df, cluster_means, n_clusters, verbose = verbose_clustering) 
    
    merged_data <- create_merged_data(gene_regulon, cluster_results_df_typed, verbose = verbose_clustering)
    
    if (!is.null(merged_data) && nrow(merged_data) > 0 && ncol(merged_data) > 1) { 
      plot_heatmap(merged_data, paste0(method, "_", effective_gsva_method), output_dir, verbose = verbose_clustering) 
    } else {
      if(verbose_clustering) warning(paste("合并后的数据为空或不足以绘制热图 for method", method))
    }
    
    all_results[[method]] <- cluster_results_df_typed
    all_merged_data[[method]] <- merged_data
  }
  
  # 保存完整聚类输出
  clustering_output_complete <- list(
    cluster_results = all_results,
    merged_data = all_merged_data,
    gsva_raw_output = analysis_results,
    clustering_params = list(
      methods = clustering_methods,
      n_clusters = n_clusters,
      gsva_method = gsva_method,
      effective_method = effective_gsva_method
    )
  )
  
  output_file <- file.path(output_dir, "clustering_output_complete.rds")
  tryCatch({
    saveRDS(clustering_output_complete, output_file)
    if(verbose_clustering) cat("💾 完整聚类输出已保存:", output_file, "\n")
  }, error = function(e) {
    warning(paste("保存完整聚类输出失败:", e$message))
  })
  
  # 新增：保存 hierarchical 聚类结果为 CSV
  if (!is.null(clustering_output_complete$cluster_results$hierarchical)) {
    hierarchical_clusters_file <- file.path(output_dir, "hierarchical_clusters.csv")
    write.csv(clustering_output_complete$cluster_results$hierarchical, 
              hierarchical_clusters_file, 
              row.names = TRUE)
    if(verbose_clustering) cat("💾 层次聚类结果已保存:", hierarchical_clusters_file, "\n")
  }
  
  return(clustering_output_complete)
} 

#' 生存分析函数 - LANCET期刊风格版本
perform_survival_analysis <- function(
    cluster_results,
    survival_data_path,
    output_dir,
    verbose = TRUE
) {
  # 加载必要的包
  required_packages <- c("survival", "survminer", "dplyr", "ggplot2", "grid", "tibble", "ggpp") 
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # 读取生存数据
  if(verbose) print(paste("📖 读取生存数据从:", survival_data_path))
  survival_data_raw <- read.csv(survival_data_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  if(verbose) print(paste("📊 读取后，生存数据原始维度:", paste(dim(survival_data_raw), collapse="x")))
  
  # 处理重复行名
  if(any(duplicated(rownames(survival_data_raw)))) {
    if(verbose) warning("生存数据文件第一列 (用作行名) 包含重复ID。将移除重复项，保留第一个出现的。")
    survival_data_raw <- survival_data_raw[!duplicated(rownames(survival_data_raw)), , drop = FALSE]
    if(verbose) print(paste("去重后，生存数据维度:", paste(dim(survival_data_raw), collapse="x")))
  }
  
  original_survival_ids_as_read <- rownames(survival_data_raw) 
  if(verbose && nrow(survival_data_raw) > 0) {
    print(paste("原始生存数据样本ID (前5个):", paste(head(original_survival_ids_as_read,5), collapse=", ")))
  }
  
  # 替换连字符为点号
  if(nrow(survival_data_raw) > 0) {
    rnames_hyphen_to_dot <- gsub("-", ".", rownames(survival_data_raw))
    rownames(survival_data_raw) <- rnames_hyphen_to_dot
    if(verbose) print(paste("替换连字符后生存ID (前5个):", paste(head(rownames(survival_data_raw),5), collapse=", ")))
  }
  
  # 完整样本ID标准化
  if(nrow(survival_data_raw) > 0) {
    rownames(survival_data_raw) <- normalize_sample_ids(rownames(survival_data_raw)) 
    if(verbose) print(paste("标准化后生存数据样本ID (前5个):", paste(head(rownames(survival_data_raw),5), collapse=", ")))
  }
  
  cli_processed <- survival_data_raw
  
  # 检查生存数据中是否包含必要的列
  required_cols <- c("futime", "fustat") 
  missing_cols <- setdiff(required_cols, colnames(cli_processed))
  if(length(missing_cols) > 0) {
    stop(paste("❌ 生存数据缺少必要的列:", paste(missing_cols, collapse=", ")))
  }
  
  # 确保生存数据列为数值型
  for(col in required_cols) {
    if(!is.numeric(cli_processed[[col]])) {
      if(verbose) warning(paste("生存数据列 '", col, "' 非数值型，尝试转换。", sep=""))
      original_na_count <- sum(is.na(cli_processed[[col]]))
      cli_processed[[col]] <- as.numeric(gsub(",", ".", as.character(cli_processed[[col]])))
      new_na_count <- sum(is.na(cli_processed[[col]]))
      if (new_na_count > original_na_count && verbose) {
        warning(paste("转换 '", col, "' 为数值型时引入了NA值。", sep=""))
      }
    }
  }
  
  # 过滤无效生存数据
  cli_processed <- cli_processed[!is.na(cli_processed$futime) & !is.na(cli_processed$fustat) & cli_processed$futime > 0, ]
  if(nrow(cli_processed) > 0) {
    if(verbose) print(paste("🧹 移除NA和futime<=0后生存数据行数:", nrow(cli_processed)))
  } else {
    stop("❌ 在futime/fustat中移除NA或futime<=0后，无有效生存数据。")
  }
  
  # TCGA样本ID兼容性处理
  if (length(cluster_results) > 0 && !is.null(names(cluster_results)) && length(names(cluster_results)) > 0 && nrow(cli_processed) > 0) {
    first_method_name <- names(cluster_results)[1] 
    first_method_results <- cluster_results[[first_method_name]]
    
    if (!is.null(first_method_results) && nrow(first_method_results) > 0) {
      example_cluster_id_normalized <- rownames(first_method_results)[1] 
      
      is_cluster_id_sample_level_tcga <- (
        (nchar(example_cluster_id_normalized) == 15 && grepl("^TCGA\\.[A-Z0-9]{2}\\.[A-Z0-9]{4}\\.\\d{2}$", example_cluster_id_normalized)) ||
          (nchar(example_cluster_id_normalized) == 16 && grepl("^TCGA\\.[A-Z0-9]{2}\\.[A-Z0-9]{4}\\.\\d{2}[A-Z]$", example_cluster_id_normalized))
      )
      
      if (is_cluster_id_sample_level_tcga) {
        current_survival_ids_normalized <- rownames(cli_processed)
        new_survival_ids_normalized <- current_survival_ids_normalized
        
        patient_level_ids_idx <- which(
          nchar(current_survival_ids_normalized) == 12 & 
            grepl("^TCGA\\.[A-Z0-9]{2}\\.[A-Z0-9]{4}$", current_survival_ids_normalized)
        )
        
        if (length(patient_level_ids_idx) > 0) {
          if(verbose) print(paste("🔗 检测到", length(patient_level_ids_idx), "个12字符TCGA格式的生存ID。由于聚类ID似乎是样本级别，将尝试为这些生存ID附加'.01'。"))
          
          new_survival_ids_normalized[patient_level_ids_idx] <- paste0(current_survival_ids_normalized[patient_level_ids_idx], ".01")
          rownames(cli_processed) <- new_survival_ids_normalized
          
          if(verbose && length(patient_level_ids_idx) > 0) print(paste("附加'.01'后生存数据样本ID (前5个):", paste(head(rownames(cli_processed),5), collapse=", ")))
        }
      }
    }
  }
  
  survival_results_list <- list() 
  
  if(length(cluster_results) == 0 || all(sapply(cluster_results, function(x) is.null(x) || nrow(x) == 0 ))) {
    if(verbose) warning("❌ 聚类结果为空或null或不包含任何样本，跳过生存分析。")
    return(list(error = "聚类结果为空或null或不包含样本"))
  }
  
  for(method_iter_name in names(cluster_results)) { 
    if(verbose) print(paste("🧬 正在进行", method_iter_name, "的生存分析..."))
    current_clusters_df <- cluster_results[[method_iter_name]] 
    
    if(is.null(current_clusters_df) || nrow(current_clusters_df) == 0) {
      if(verbose) warning(paste("方法 '", method_iter_name, "' 无聚类数据，跳过。", sep=""))
      survival_results_list[[method_iter_name]] <- list(error = paste("方法 '", method_iter_name, "' 无聚类数据", sep=""))
      next
    }
    
    if(verbose) {
      print(paste("聚类结果维度 ('", method_iter_name, "'): ", paste(dim(current_clusters_df), collapse="x"), sep=""))
      if("cluster" %in% colnames(current_clusters_df)){
        print(paste("聚类类别: ", paste(sort(unique(as.character(current_clusters_df$cluster))), collapse=", "), sep="")) 
      } else {
        warning(paste("列 'cluster' 在方法 '", method_iter_name, "' 的聚类结果中未找到。", sep=""))
        survival_results_list[[method_iter_name]] <- list(error = "列 'cluster' 在聚类数据中未找到")
        next
      }
    }
    
    if(nrow(cli_processed) == 0) {
      if(verbose) warning("生存数据在过滤和标准化后为空，无法查找共同样本。")
      survival_results_list[[method_iter_name]] <- list(error = "生存数据为空")
      next
    }
    
    sameSample <- intersect(rownames(current_clusters_df), rownames(cli_processed))
    if(verbose) print(paste("🔍 共同样本数量:", length(sameSample)))
    
    if(length(sameSample) == 0) {
      warning(paste("方法 '", method_iter_name, "' 的聚类结果与生存数据无共同样本。", sep=""))
      if(verbose){
        cat("--- 调试信息: 无共同样本 ---\n")
        cat("聚类样本ID (前5个):"); print(head(rownames(current_clusters_df), 5))
        cat("生存数据样本ID (前5个):"); print(head(rownames(cli_processed), 5))
        cat("---------------------------------\n")
      }
      
      survival_results_list[[method_iter_name]] <- list(
        error = "无共同样本",
        cluster_samples_head_norm = head(rownames(current_clusters_df), 20),
        survival_samples_head_fully_processed = head(rownames(cli_processed), 20),
        survival_samples_head_as_read_dedup = head(original_survival_ids_as_read, 20)
      )
      next
    }
    
    if (!"cluster" %in% colnames(current_clusters_df)) {
      if(verbose) warning(paste0("方法 '", method_iter_name, "' 的聚类结果中缺少 'cluster' 列。", sep=""))
      survival_results_list[[method_iter_name]] <- list(error = "缺少 'cluster' 列")
      next
    }
    
    rt_surv <- cbind(cli_processed[sameSample, c("futime", "fustat"), drop = FALSE], 
                     cluster_assignment = current_clusters_df[sameSample, "cluster", drop = FALSE]) 
    colnames(rt_surv)[ncol(rt_surv)] <- "regulon_cluster" 
    
    if(verbose) print(paste("📊 合并后数据维度:", paste(dim(rt_surv), collapse="x")))
    
    rt_surv_clean <- na.omit(rt_surv) 
    if(verbose) print(paste("🧹 清除NA后数据行数:", nrow(rt_surv_clean)))
    
    if(nrow(rt_surv_clean) == 0) {
      if(verbose) warning("清除NA后没有剩余数据用于生存分析。")
      survival_results_list[[method_iter_name]] <- list(
        error = "清除NA后无有效数据",
        raw_merged_data_head = head(rt_surv)
      )
      next 
    }
    
    rt_surv_clean$regulon_cluster <- factor(rt_surv_clean$regulon_cluster)
    cluster_counts <- table(rt_surv_clean$regulon_cluster)
    if(verbose) print(paste("📈 各类别样本数:", paste(names(cluster_counts), cluster_counts, sep="=", collapse=", ")))
    
    min_samples_per_group <- 1 
    if(any(cluster_counts < min_samples_per_group) || length(unique(rt_surv_clean$regulon_cluster)) < 2) { 
      if(verbose) warning(paste("存在样本数过少 (少于", min_samples_per_group, "个) 的类别，或少于两个有效类别用于生存分析。", sep=""))
      survival_results_list[[method_iter_name]] <- list(
        error = paste("样本数过少或类别不足", sep=""),
        cluster_counts = cluster_counts,
        data_head = head(rt_surv_clean)
      )
      next 
    }
    
    if(verbose) {
      print(paste("⏰ 生存时间范围:", min(rt_surv_clean$futime, na.rm=TRUE), "-", max(rt_surv_clean$futime, na.rm=TRUE)))
      print(paste("📊 生存状态计数:", paste(names(table(rt_surv_clean$fustat)), table(rt_surv_clean$fustat), sep="=", collapse=", ")))
    }
    
    tryCatch({
      if (nlevels(rt_surv_clean$regulon_cluster) < 2) { 
        stop("需要至少两个聚类组进行生存差异分析。")
      }
      
      # 创建生存对象
      fit <- survfit(Surv(futime, fustat) ~ regulon_cluster, data = rt_surv_clean)
      diff_result <- NULL
      pValue_val <- NA
      pValue_str <- "NA"
      
      if (nlevels(rt_surv_clean$regulon_cluster) > 1) { 
        diff_result <- tryCatch(survdiff(Surv(futime, fustat) ~ regulon_cluster, data = rt_surv_clean), error = function(e) NULL)
        if(!is.null(diff_result) && !is.null(diff_result$chisq)){
          pValue_val <- 1 - pchisq(diff_result$chisq, df = nlevels(rt_surv_clean$regulon_cluster) - 1)
          pValue_str <- ifelse(pValue_val < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", pValue_val)))
        } else {
          if(verbose) warning("survdiff未能成功计算，P值将为NA。")
        }
      } 
      
      # 计算两两比较的p值（与model.r保持一致）
      ps <- pairwise_survdiff(Surv(futime, fustat) ~ regulon_cluster, data = rt_surv_clean)
      
      # 创建p值表格（与model.r保持一致）
      addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                               round(ps$p.value, 3))))
      addTab[is.na(addTab)] <- "-"
      
      # 智能配色匹配 - 与model.r完全一致的逻辑
      cluster_levels <- levels(factor(rt_surv_clean$regulon_cluster))
      n_clusters <- length(cluster_levels)
      
      # 智能配色匹配 - 根据聚类名称匹配对应颜色
      cluster_colors <- rep(NA, n_clusters)
      names(cluster_colors) <- cluster_levels
      
      # 遍历每个聚类，根据名称匹配颜色
      for(cluster_name in cluster_levels) {
        cluster_str <- as.character(cluster_name)
        
        # 检查是否匹配特定的聚类名称模式
        if(grepl("low|Low|LOW", cluster_str)) {
          cluster_colors[cluster_name] <- LANCET_COLORS$cluster_colors["cluster_low"]
        } else if(grepl("medium|Medium|MEDIUM|median|Median|MEDIAN", cluster_str)) {
          cluster_colors[cluster_name] <- LANCET_COLORS$cluster_colors["cluster_medium"]
        } else if(grepl("high|High|HIGH", cluster_str)) {
          cluster_colors[cluster_name] <- LANCET_COLORS$cluster_colors["cluster_high"]
        } else if(cluster_str %in% names(LANCET_COLORS$cluster_colors)) {
          # 如果聚类名称直接在配色方案中存在
          cluster_colors[cluster_name] <- LANCET_COLORS$cluster_colors[cluster_str]
        } else {
          # 如果没有匹配到，按顺序分配
          cluster_index <- which(cluster_levels == cluster_name)
          if(cluster_index <= length(LANCET_COLORS$survival_colors)) {
            cluster_colors[cluster_name] <- LANCET_COLORS$survival_colors[cluster_index]
          } else if(cluster_index <= length(LANCET_COLORS$extended_colors)) {
            cluster_colors[cluster_name] <- LANCET_COLORS$extended_colors[cluster_index]
          } else {
            cluster_colors[cluster_name] <- rainbow(n_clusters)[cluster_index]
          }
        }
      }
      
      # 移除名称避免冲突
      names(cluster_colors) <- NULL
      
      # 创建配色映射
      sorted_cluster_levels <- sort(cluster_levels)
      custom_colors <- setNames(cluster_colors[match(sorted_cluster_levels, cluster_levels)], sorted_cluster_levels)
      
      # 打印配色信息用于调试
      if(verbose) {
        message("生存分析配色映射:")
        for(cluster_name in cluster_levels) {
          color_val <- cluster_colors[match(cluster_name, cluster_levels)]
          message("  ", cluster_name, " -> ", color_val)
        }
      }
      
      km_plot <- ggsurvplot(
        fit,
        data = rt_surv_clean,
        pval = FALSE,
        conf.int = FALSE,
        risk.table = TRUE,
        palette = custom_colors,
        title = paste("Overall Survival Analysis -", method_iter_name),  
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
        
        # 字体大小设置
        font.main = c(14, "bold", "black"),        
        font.x = c(14, "plain", "black"),          
        font.y = c(14, "plain", "black"),          
        font.legend = c(14, "plain", "black"),     
        font.tickslab = c(14, "plain", "black"),   
        
        # 风险表字体设置
        tables.col = "strata",
        risk.table.fontsize = 6,  
        risk.table.title.fontsize = 12,
        
        # 精细控制Number at risk部分的字体
        risk.table.text.size = 12,    # 分组标签（cluster_high等）字体大小
        risk.table.digit.size = 6,    # 数字（116, 2等）字体大小
        
        # 额外的字体微调
        risk.table.col.text.size = 12,   # 列标签字体大小
        risk.table.col.text.face = "bold"  # 列标签字体加粗
        
        
      )
      
      # 修改p值表格位置到右上方
      if(!is.null(ps$p.value) && any(!is.na(ps$p.value))) {
        x_pos <- max(rt_surv_clean$futime) * 0.95  # 移到右侧
        y_pos <- 1.0  # 移到顶部
        df <- tibble(x = x_pos,
                     y = y_pos,
                     tb = list(addTab))
        km_plot$plot <- km_plot$plot +
          ggpp::geom_table(data = df,
                           aes(x = x, y = y, label = tb),
                           table.rownames = TRUE,
                           size = 5,
                           hjust = 1,   # 右对齐
                           vjust = 1)   # 顶部对齐
      }
      # 保存图片 - 调整宽度从10改为8
      base_filename <- file.path(output_dir, paste0("survival_", method_iter_name))
      if(!dir.exists(dirname(base_filename))) dir.create(dirname(base_filename), recursive = TRUE)
      
      # 保存为PDF
      pdf_filename <- paste0(base_filename, ".pdf")
      tryCatch({
        pdf(pdf_filename, width = 7, height = 7)
        print(km_plot)
        dev.off()
        if(verbose) cat("   📄 PDF保存成功: ", pdf_filename, "\n")
      }, error = function(e) {
        if(verbose) warning("PDF保存失败: ", e$message)
      })
      
      # 2. 保存为PNG
      png_filename <- paste0(base_filename, ".png")
      tryCatch({
        png(png_filename, width = 1000, height =800, res = 300, bg = "white")
        print(km_plot)
        dev.off()
        if(verbose) cat("   🖼️ PNG保存成功: ", png_filename, "\n")
      }, error = function(e) {
        if(verbose) warning("PNG保存失败: ", e$message)
      })
      
      # 3. 尝试保存为TIFF（如果支持）
      tiff_filename <- paste0(base_filename, ".tiff")
      tryCatch({
        if(capabilities("tiff")) {
          tiff(tiff_filename, width = 1000, height =800, res = 300, bg = "white", compression = "lzw")
          print(km_plot)
          dev.off()
          if(verbose) cat("   📊 TIFF保存成功: ", tiff_filename, "\n")
        }
      }, error = function(e) {
        if(verbose) warning("TIFF保存失败: ", e$message)
      })
      
      survival_results_list[[method_iter_name]] <- list(
        fit = fit,
        diff = diff_result,
        pValue_numeric = pValue_val,
        pValue_string = pValue_str,
        plot_object = km_plot,
        sample_counts = cluster_counts
      )
      
    }, error = function(e_surv) {
      if(verbose) print(paste("❌ 生存分析出错 for method '", method_iter_name, "': ", e_surv$message, sep=""))
      survival_results_list[[method_iter_name]] <- list(
        error_message = e_surv$message,
        data_head_at_error = if(exists("rt_surv_clean")) head(rt_surv_clean) else NULL
      )
    })
  }
  
  # 保存完整生存分析结果
  survival_output_file <- file.path(output_dir, "survival_output_complete.rds")
  tryCatch({
    saveRDS(survival_results_list, survival_output_file)
    if(verbose) cat("💾 完整生存分析输出已保存:", survival_output_file, "\n")
  }, error = function(e) {
    warning(paste("保存完整生存分析输出失败:", e$message))
  })
  
  return(survival_results_list)
} 

# ====================================================================
# 辅助函数
# ====================================================================

#' 计算cluster均值 - 通用版本
calculate_cluster_means <- function(gene_regulon, cluster_results_df, verbose = FALSE) { 
  cluster_means <- list()
  cluster_results_df$cluster <- as.character(cluster_results_df$cluster)
  unique_clusters <- unique(cluster_results_df$cluster)
  
  if(ncol(gene_regulon) == 0) {
    if(verbose) warning("gene_regulon数据为空，无法计算聚类均值。")
    return(cluster_means)
  }
  
  # 获取所有数值型列
  numeric_cols <- colnames(gene_regulon)[sapply(gene_regulon, is.numeric)]
  if(length(numeric_cols) == 0) {
    if(verbose) warning("gene_regulon中没有数值型列，无法计算聚类均值。")
    return(cluster_means)
  }
  
  for(i in unique_clusters) {
    cluster_samples <- rownames(cluster_results_df)[cluster_results_df$cluster == i]
    if (length(cluster_samples) == 0) {
      cluster_means[[as.character(i)]] <- list("overall_mean" = 0, "sample_count" = 0) 
      next
    }
    
    # 找到在gene_regulon中存在的样本
    samples_in_regulon_data <- cluster_samples[cluster_samples %in% rownames(gene_regulon)]
    if (length(samples_in_regulon_data) == 0) {
      cluster_means[[as.character(i)]] <- list("overall_mean" = 0, "sample_count" = 0)
      next
    }
    
    # 计算该聚类在所有regulon上的总体均值
    cluster_data <- gene_regulon[samples_in_regulon_data, numeric_cols, drop = FALSE]
    if(nrow(cluster_data) > 0 && ncol(cluster_data) > 0) {
      # 计算每个样本的所有regulon平均值
      sample_means <- rowMeans(cluster_data, na.rm = TRUE)
      # 计算该聚类的总体均值
      overall_cluster_mean <- mean(sample_means, na.rm = TRUE)
      
      cluster_means[[as.character(i)]] <- list(
        "overall_mean" = if(is.finite(overall_cluster_mean)) overall_cluster_mean else 0,
        "sample_count" = length(samples_in_regulon_data),
        "regulon_count" = ncol(cluster_data)
      )
    } else {
      cluster_means[[as.character(i)]] <- list("overall_mean" = 0, "sample_count" = 0, "regulon_count" = 0)
    }
    
    if(verbose) {
      cat("聚类", i, ": 样本数=", length(samples_in_regulon_data), 
          ", Regulon总体均值=", round(cluster_means[[as.character(i)]]$overall_mean, 4), "\n")
    }
  }
  return(cluster_means)
}

#' 分配cluster类型 - 通用版本
assign_cluster_types <- function(cluster_results_df, cluster_means, n_target_clusters, verbose = FALSE) { 
  cluster_results_df$cluster <- as.character(cluster_results_df$cluster) 
  
  if (length(cluster_means) == 0 || n_target_clusters == 0) {
    if(verbose) warning("聚类均值为空或目标聚类数为0，无法分配类型。返回原始聚类。")
    cluster_results_df$cluster <- factor(cluster_results_df$cluster)
    return(cluster_results_df)
  }
  
  if (n_target_clusters == 1) {
    if(length(unique(cluster_results_df$cluster)) == 1) {
      cluster_results_df$cluster <- "cluster_1"
      cluster_results_df$cluster <- factor(cluster_results_df$cluster)
      return(cluster_results_df)
    }
  }
  
  # 提取所有聚类的总体均值
  overall_means <- sapply(cluster_means, function(x) {
    if(is.list(x) && "overall_mean" %in% names(x) && is.finite(x$overall_mean)) {
      return(x$overall_mean)
    } else {
      return(NA_real_)
    }
  })
  
  valid_means <- overall_means[!is.na(overall_means)]
  if(length(valid_means) == 0) {
    if(verbose) warning("没有有效的聚类均值，使用原始聚类ID。")
    cluster_results_df$cluster <- factor(paste0("cluster_", cluster_results_df$cluster))
    return(cluster_results_df)
  }
  
  original_cluster_ids <- names(cluster_means)
  n_clusters_from_means <- length(original_cluster_ids)
  
  # 根据总体均值进行分型
  if (n_target_clusters == 2) {
    # 二分类：高、低
    if(length(valid_means) >= 2) {
      high_cluster <- names(valid_means)[which.max(valid_means)]
      low_cluster <- names(valid_means)[which.min(valid_means)]
      
      new_cluster_names_map <- stats::setNames(rep("cluster_medium", n_clusters_from_means), original_cluster_ids)
      new_cluster_names_map[high_cluster] <- "cluster_high"
      new_cluster_names_map[low_cluster] <- "cluster_low"
    } else {
      new_cluster_names_map <- stats::setNames(paste0("cluster_", 1:n_clusters_from_means), original_cluster_ids)
    }
  } else if (n_target_clusters == 3) {
    # 三分类：低、中、高
    if(length(valid_means) >= 3) {
      sorted_means <- sort(valid_means, decreasing = TRUE)
      high_cluster <- names(sorted_means)[1]
      low_cluster <- names(sorted_means)[length(sorted_means)]
      new_cluster_names_map <- stats::setNames(rep("cluster_medium", n_clusters_from_means), original_cluster_ids)
      new_cluster_names_map[high_cluster] <- "cluster_high"
      new_cluster_names_map[low_cluster] <- "cluster_low"
    } else if(length(valid_means) == 2) {
      high_cluster <- names(valid_means)[which.max(valid_means)]
      low_cluster <- names(valid_means)[which.min(valid_means)]
      new_cluster_names_map <- stats::setNames(rep("cluster_medium", n_clusters_from_means), original_cluster_ids)
      new_cluster_names_map[high_cluster] <- "cluster_high"
      new_cluster_names_map[low_cluster] <- "cluster_low"
    } else {
      new_cluster_names_map <- stats::setNames(paste0("cluster_", 1:n_clusters_from_means), original_cluster_ids)
    }
  }
  
  if(verbose) {
    cat("🔢 聚类分型结果：\n")
    for(cluster_id in names(new_cluster_names_map)) {
      mean_val <- if(cluster_id %in% names(overall_means)) round(overall_means[cluster_id], 4) else "NA"
      cat("   原始聚类", cluster_id, "-> ", new_cluster_names_map[cluster_id], " (均值:", mean_val, ")\n")
    }
  }
  
  # 应用新的聚类名称
  original_clusters_in_df <- as.character(cluster_results_df$cluster)
  mapped_names <- new_cluster_names_map[original_clusters_in_df] 
  
  na_indices <- is.na(mapped_names)
  if (any(na_indices)) {
    if(verbose) warning("部分聚类ID在映射中未找到，使用默认命名。")
    mapped_names[na_indices] <- paste0("cluster_", original_clusters_in_df[na_indices]) 
  }
  
  cluster_results_df$cluster <- factor(mapped_names)
  return(cluster_results_df)
}

#' 创建合并数据
create_merged_data <- function(gene_regulon, cluster_results_df_typed, verbose = FALSE) { 
  if(!("cluster" %in% colnames(cluster_results_df_typed))) {
    stop("在 create_merged_data: cluster_results_df_typed 必须包含 'cluster' 列。")
  }
  if(is.null(rownames(gene_regulon)) || is.null(rownames(cluster_results_df_typed))) {
    stop("在 create_merged_data: gene_regulon 和 cluster_results_df_typed 必须有行名 (样本ID)。")
  }
  
  cluster_info_to_merge <- cluster_results_df_typed[, "cluster", drop = FALSE]
  
  common_samples <- intersect(rownames(gene_regulon), rownames(cluster_info_to_merge))
  if(length(common_samples) == 0) {
    if(verbose) warning("create_merged_data: gene_regulon 和 cluster_info 无共同样本。返回空的合并数据。")
    return(data.frame()) 
  }
  
  merged_regulon <- merge(gene_regulon[common_samples, , drop=FALSE], 
                          cluster_info_to_merge[common_samples, "cluster", drop = FALSE], 
                          by = "row.names", all = FALSE) 
  
  if(nrow(merged_regulon) > 0) { 
    rownames(merged_regulon) <- merged_regulon$Row.names 
    merged_regulon <- merged_regulon[, -which(colnames(merged_regulon) == "Row.names")] 
  } else {
    if(verbose) warning("create_merged_data: 合并后无数据。")
    return(data.frame())
  }
  
  if("cluster.y" %in% colnames(merged_regulon) && !"cluster" %in% colnames(merged_regulon)){
    colnames(merged_regulon)[colnames(merged_regulon) == "cluster.y"] <- "cluster"
  } else if ("cluster.x" %in% colnames(merged_regulon) && "cluster.y" %in% colnames(merged_regulon)){
    if(verbose) warning("create_merged_data: 'cluster' 列名冲突。检查输入数据。")
    merged_regulon$cluster <- merged_regulon$cluster.y
    merged_regulon$cluster.x <- NULL
    merged_regulon$cluster.y <- NULL
  }
  
  return(merged_regulon)
}

#' 绘制热图 - LANCET期刊风格版本
plot_heatmap <- function(merged_data, method_name, output_dir, verbose = FALSE) {
  if (is.null(merged_data) || nrow(merged_data) == 0 || ncol(merged_data) <= 1) {
    if (verbose) warning("热图: 合并数据为空或无特征数据，跳过热图。")
    return(NULL)
  }
  if (!"cluster" %in% colnames(merged_data)) {
    if (verbose) warning("热图: 'cluster' 列缺失，跳过热图。")
    return(NULL)
  }
  
  ## ---- 数据整理（与原代码一致） ----
  merged_data$cluster <- factor(merged_data$cluster,
                                levels = c("cluster_low", "cluster_medium", "cluster_high"))
  sort_index           <- order(merged_data$cluster)
  gene_regulon_sorted  <- merged_data[sort_index, , drop = FALSE]
  feature_cols         <- setdiff(colnames(gene_regulon_sorted), "cluster")
  numeric_feature_cols <- feature_cols[sapply(gene_regulon_sorted[, feature_cols, drop = FALSE], is.numeric)]
  expr_data            <- gene_regulon_sorted[, numeric_feature_cols, drop = FALSE]
  Type_regulon         <- gene_regulon_sorted[, "cluster", drop = FALSE]
  data_h               <- as.matrix(t(expr_data))
  
  if (nrow(data_h) == 0 || ncol(data_h) == 0) {
    if (verbose) warning("热图: 转置后数据为空，跳过热图。")
    return(NULL)
  }
  
  ## ---- 行标准化（手动 z-score） ----
  data_h <- t(scale(t(data_h)))          # 行 z-score
  data_h <- pmin(pmax(data_h, -3), 3)   # 截断到 [-3, 3]
  
  ## ---- 颜色与注释（与原代码一致） ----
  custom_colors <- LANCET_COLORS$heatmap_colors
  cluster_levels      <- levels(Type_regulon$cluster)
  num_clusters_in_data <- length(cluster_levels)
  
  final_palette <- character(0)
  if (num_clusters_in_data <= length(LANCET_COLORS$cluster_colors)) {
    final_palette <- LANCET_COLORS$cluster_colors[cluster_levels]
  } else {
    final_palette <- c(LANCET_COLORS$cluster_colors,
                       LANCET_COLORS$extended_colors[seq_len(num_clusters_in_data -
                                                               length(LANCET_COLORS$cluster_colors))])
  }
  cluster_color_map <- stats::setNames(final_palette, cluster_levels)
  ann_colors <- list(cluster = cluster_color_map)
  
  ## ---- gap 计算（与原代码一致） ----
  cluster_counts_table <- table(Type_regulon$cluster)
  gap_positions <- NULL
  if (length(cluster_counts_table) > 1 && ncol(data_h) > length(cluster_counts_table)) {
    gap_positions <- cumsum(cluster_counts_table)[-length(cluster_counts_table)]
  }
  
  ## ---- 绘图 ----
  p1 <- pheatmap::pheatmap(
    data_h,
    annotation_col    = Type_regulon,
    annotation_colors = ann_colors,
    color             = colorRampPalette(c("#2E8B57", "#FFFFFF", "#DC143C"))(100),
    breaks            = seq(-3, 3, length.out = 101),  # 关键：固定颜色条范围
    cluster_cols      = FALSE,
    cluster_rows      = FALSE,
    show_colnames     = FALSE,
    show_rownames     = ifelse(nrow(data_h) < 50, TRUE, FALSE),
    fontsize          = 12,
    fontsize_row      = max(10, min(14, 14 - floor(nrow(data_h)/15))),
    fontsize_col      = 12,
    main              = paste("Regulon Activation Pattern -", method_name),
    fontface_row      = "bold",
    gaps_col          = gap_positions,
    scale             = "none",          # 已手动 scale，不再让 pheatmap 处理
    na_col            = "grey80"
  )
  
  ## ---- 保存 ----
  pdf_file <- file.path(output_dir, paste0("heatmap_", method_name, ".pdf"))
  plot_height <- max(4, min(10, nrow(data_h) * 0.10 ))
  pdf(pdf_file, width = 7, height = plot_height)
  grid::grid.draw(p1$gtable)
  dev.off()
  if (verbose) cat("🎨 热图已保存到:", pdf_file, "\n")
  
  invisible(p1)
}