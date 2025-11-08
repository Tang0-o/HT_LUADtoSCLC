# SCLC Differential Expression and Pathway Analysis Functions
# Version 2.0
#################################################################################

# Load required packages
load_required_packages <- function() {
  required_packages <- c(
    "ggplot2", "dplyr", "tidyr", "openxlsx", "RColorBrewer", "pheatmap",
    "limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", 
    "DOSE", "enrichplot", "ggrepel", "msigdbr", "fgsea",
    "readr", "stringr", "tibble", "igraph"  # Added igraph for network visualization
  )
  
  for(pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", 
                     "DOSE", "enrichplot", "msigdbr", "fgsea")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
    library(pkg, character.only = TRUE)
  }
}

# Get default SCLC signatures
get_default_sclc_signatures <- function() {
  list(
    # 1. SCLC-specific Neuroendocrine markers and pathways
    "SCLC_specific" = c(
      # Core NE markers
      "CHGA", "CHGB", "SYP", "NCAM1", "ENO2",
      # Key transcription factors
      "ASCL1", "NEUROD1", "INSM1",
      # Additional NE markers
      "DLL3", "GRP", "CALCA", "SCG2",
      # Notch pathway suppression
      "HES1",
      # DNA damage response
      "PARP1", "SLFN11",
      # Cell cycle and survival
      "MYC", "MYCL", "MYCN", "BCL2",
      # Immune evasion
      "CD47", "PVR"
    ),
    
    # 2. LUAD-specific markers
    "LUAD_specific" = c(
      # RTK/RAS/RAF pathway
      "EGFR", "KRAS", "BRAF", "NF1", "ERBB2",
      # Differentiation markers
      "NKX2-1", "NAPSA", "SFTPB", "SFTPC", "FOXA1", "FOXA2",
      # Additional markers
      "CEACAM5", "MUC5AC", "KRT7",
      # Growth and survival
      "IGF1R", "MET", "YAP1", "MDM2"
    ),
    
    # 3. Transformation markers
    "transformation_up" = c(
      # Cell cycle
      "E2F", "CCNE1", "WEE1",
      # Neuronal differentiation
      "ASCL1", "INSM1", "SEZ6", "CHGA",
      # DNA replication
      "MCM2", "POLE", "RFC",
      # Neural activity
      "NLRI", "GRIK3", "GRM4",
      # Epigenetic regulators
      "EZH2", "NUSAP1", "TTK", "UBE2C", "DNMT3A", "HDAC7", "MYC"
    ),
    
    # 4. Epigenetic reprogramming
    "Epigenetic_reprogramming" = c(
      # DNA methylation
      "DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET3",
      # Histone modification
      "KAT2A", "HDAC2", "EZH2",
      # Chromatin remodeling
      "SMARCA4"
    )
  )
}

# Read and validate input data
read_input_data <- function(exp_file, group_file, data_type = "FPKM") {
  message("Reading input data...")
  
  # Read expression data
  if (is.character(exp_file)) {
    if (endsWith(exp_file, ".rds")) {
      predictions <- readRDS(exp_file)
      if (is.list(predictions) && all(c("expression_matrix", "class") %in% names(predictions))) {
        exp_matrix <- predictions$expression_matrix
        group_data <- data.frame(
          group = predictions$class,
          row.names = colnames(predictions$expression_matrix)
        )
        message("Successfully read predictions RDS file with expression matrix and class labels")
        return(list(exp_matrix = exp_matrix, group_data = group_data))
      } else {
        exp_matrix <- predictions
      }
    } else if (endsWith(exp_file, ".csv")) {
      exp_matrix <- read.csv(exp_file, row.names = 1)
    } else if (endsWith(exp_file, ".txt") || endsWith(exp_file, ".tsv")) {
      exp_matrix <- read.delim(exp_file, row.names = 1)
    } else {
      stop("Unsupported expression file format. Please use .rds, .csv, or .txt/.tsv")
    }
  } else if (is.matrix(exp_file) || is.data.frame(exp_file)) {
    exp_matrix <- as.matrix(exp_file)
  } else {
    stop("Invalid expression data format")
  }
  
  # Read group data if not already read from predictions RDS
  if (!exists("group_data")) {
    if (is.character(group_file)) {
      if (endsWith(group_file, ".rds")) {
        group_data <- readRDS(group_file)
        if ("cluster_results" %in% names(group_data) && "hierarchical" %in% names(group_data$cluster_results)) {
          group_data <- data.frame(
            group = group_data$cluster_results$hierarchical,
            row.names = names(group_data$cluster_results$hierarchical)
          )
          message("Successfully read clustering output RDS with hierarchical clustering results")
        }
      } else if (endsWith(group_file, ".csv")) {
        group_data <- read.csv(group_file, row.names = 1)
      } else if (endsWith(group_file, ".txt") || endsWith(group_file, ".tsv")) {
        group_data <- read.delim(group_file, row.names = 1)
      } else {
        stop("Unsupported group file format. Please use .rds, .csv, or .txt/.tsv")
      }
    } else if (is.data.frame(group_file)) {
      group_data <- group_file
    } else if (is.factor(group_file) || is.character(group_file)) {
      group_data <- data.frame(
        group = group_file,
        row.names = names(group_file)
      )
    } else {
      stop("Invalid group data format")
    }
  }
  
  # Basic validation
  if (ncol(exp_matrix) == 0 || nrow(exp_matrix) == 0) {
    stop("Expression matrix is empty")
  }
  
  if (nrow(group_data) == 0) {
    stop("Group data is empty")
  }
  
  # Convert group data to proper format if needed
  if (!"group" %in% colnames(group_data)) {
    if ("cluster" %in% colnames(group_data)) {
      group_data$group <- group_data$cluster
    } else {
      group_data$group <- group_data[,1]
    }
  }
  
  # Print data information
  message("\nData reading complete:")
  message("Expression matrix dimensions: ", nrow(exp_matrix), " x ", ncol(exp_matrix))
  message("Number of samples in group data: ", nrow(group_data))
  message("Group levels: ", paste(unique(group_data$group), collapse = ", "))
  message("Samples per group:")
  print(table(group_data$group))
  
  return(list(
    exp_matrix = exp_matrix,
    group_data = group_data
  ))
} 

# Enhanced data preprocessing function
preprocess_expression_data <- function(exp_matrix, data_type = "FPKM", min_samples = 3) {
  message("Preprocessing expression data...")
  
  # 1. 移除包含NaN的行
  exp_matrix <- exp_matrix[complete.cases(exp_matrix), ]
  message("Genes after NaN filtering: ", nrow(exp_matrix))
  
  # 2. 移除零值基因
  non_zero_samples <- rowSums(exp_matrix > 0)
  keep_genes <- non_zero_samples >= min_samples
  exp_matrix <- exp_matrix[keep_genes, , drop = FALSE]
  message("Genes after zero filtering: ", sum(keep_genes))
  
  if (data_type == "FPKM" || data_type == "TPM") {
    # Log2转换并增加伪计数
    exp_matrix <- log2(exp_matrix + 1)
    
    # 低表达基因过滤
    mean_expr <- rowMeans(exp_matrix, na.rm = TRUE)
    keep_genes <- mean_expr > quantile(mean_expr, 0.1, na.rm = TRUE)
    exp_matrix <- exp_matrix[keep_genes, , drop = FALSE]
    message("Genes after expression filtering: ", sum(keep_genes))
    
  } else if (data_type == "ZSCORE") {
    # 低方差基因过滤
    gene_var <- apply(exp_matrix, 1, function(x) var(x, na.rm = TRUE))
    keep_genes <- gene_var > quantile(gene_var, 0.1, na.rm = TRUE)
    exp_matrix <- exp_matrix[keep_genes, , drop = FALSE]
    message("Genes after variance filtering: ", sum(keep_genes))
    
    # Z分数标准化 - 使用手动方法处理NA
    exp_matrix_t <- t(exp_matrix)
    exp_matrix_scaled <- t(apply(exp_matrix_t, 2, function(x) {
      x_mean <- mean(x, na.rm = TRUE)
      x_sd <- sd(x, na.rm = TRUE)
      if (x_sd == 0) x_sd <- 1  # 防止除0
      (x - x_mean) / x_sd
    }))
    
    exp_matrix <- exp_matrix_scaled
  }
  
  # 最终检查并移除NA/Inf值
  na_genes <- rowSums(is.na(exp_matrix)) > 0
  inf_genes <- rowSums(is.infinite(exp_matrix)) > 0
  exp_matrix <- exp_matrix[!na_genes & !inf_genes, , drop = FALSE]
  
  message("Final number of genes: ", nrow(exp_matrix))
  return(exp_matrix)
}

# Enhanced sample name matching function
match_samples <- function(exp_matrix, group_data, verbose = TRUE) {
  if(verbose) message("=== Starting sample name matching ===")
  
  # Get original sample names
  exp_samples <- colnames(exp_matrix)
  group_samples <- rownames(group_data)
  
  if(verbose) {
    message("Expression matrix samples: ", length(exp_samples))
    message("Group data samples: ", length(group_samples))
    message("\nInitial group distribution:")
    print(table(group_data$group))
  }
  
  # Define matching strategies
  strategies <- list(
    "direct" = function(exp, group) list(exp = exp, group = group),
    "dash_to_dot" = function(exp, group) list(exp = gsub("-", ".", exp), group = gsub("-", ".", group)),
    "dot_to_dash" = function(exp, group) list(exp = gsub("\\.", "-", exp), group = gsub("\\.", "-", group)),
    "remove_special" = function(exp, group) list(
      exp = gsub("[^[:alnum:]]", "", exp),
      group = gsub("[^[:alnum:]]", "", group)
    ),
    "lowercase" = function(exp, group) list(exp = tolower(exp), group = tolower(group)),
    "remove_suffix" = function(exp, group) list(
      exp = sub("_[^_]+$", "", exp),
      group = sub("_[^_]+$", "", group)
    )
  )
  
  # Try each strategy
  best_match <- list(
    strategy = "none",
    n_matched = 0,
    exp_names = exp_samples,
    group_names = group_samples,
    matched_indices = integer(0)
  )
  
  for(strategy_name in names(strategies)) {
    strategy <- strategies[[strategy_name]]
    result <- strategy(exp_samples, group_samples)
    common <- intersect(result$exp, result$group)
    
    if(verbose) {
      message(sprintf("Strategy %s: matched %d samples", strategy_name, length(common)))
    }
    
    if(length(common) > best_match$n_matched) {
      # Get the indices of matched samples
      exp_indices <- match(common, result$exp)
      group_indices <- match(common, result$group)
      
      best_match <- list(
        strategy = strategy_name,
        n_matched = length(common),
        exp_names = result$exp,
        group_names = result$group,
        matched_indices = list(
          exp = exp_indices,
          group = group_indices
        )
      )
    }
  }
  
  if(best_match$n_matched == 0) {
    if(verbose) {
      message("\nSample name examples:")
      message("Expression matrix (first 5):", paste(head(exp_samples, 5), collapse = ", "))
      message("Group data (first 5):", paste(head(group_samples, 5), collapse = ", "))
    }
    stop("No matching samples found with any strategy")
  }
  
  if(verbose) {
    message(sprintf("\nBest strategy: %s with %d matches", 
                    best_match$strategy, best_match$n_matched))
  }
  
  # Get matched samples
  matched_exp_matrix <- exp_matrix[, best_match$matched_indices$exp, drop = FALSE]
  matched_group_data <- group_data[best_match$matched_indices$group, , drop = FALSE]
  
  # Verify the match
  if (!identical(colnames(matched_exp_matrix), rownames(matched_group_data))) {
    stop("Sample matching verification failed")
  }
  
  if(verbose) {
    message("\nFinal sample counts by group:")
    print(table(matched_group_data$group))
  }
  
  return(list(
    exp_matrix = matched_exp_matrix,
    group_data = matched_group_data,
    strategy_used = best_match$strategy,
    n_matched = best_match$n_matched
  ))
}

# Enhanced differential gene expression analysis
perform_deg_analysis <- function(exp_matrix, group_labels, data_type = "FPKM", 
                                 output_dir = "./results/DEG/", 
                                 comparison_name = "comparison",
                                 min_samples_per_group = 3,
                                 lfc_threshold = 0.5,
                                 pval_threshold = 0.05,
                                 target_groups = c("cluster_high", "cluster_low")) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Process group labels
  if (is.data.frame(group_labels)) {
    group_col <- if ("cluster" %in% colnames(group_labels)) {
      "cluster"
    } else if ("group" %in% colnames(group_labels)) {
      "group"
    } else {
      colnames(group_labels)[1]
    }
    group_data <- data.frame(
      group = group_labels[[group_col]],
      row.names = rownames(group_labels),
      stringsAsFactors = FALSE
    )
  } else {
    group_data <- data.frame(
      group = group_labels,
      row.names = names(group_labels),
      stringsAsFactors = FALSE
    )
  }
  
  # Print initial group information
  message("\nInitial group information:")
  print(table(group_data$group))
  
  # Filter for target groups only
  if (!all(target_groups %in% unique(group_data$group))) {
    stop(sprintf("Target groups %s not found in data. Available groups: %s", 
                 paste(target_groups, collapse = ", "), 
                 paste(unique(group_data$group), collapse = ", ")))
  }
  
  # Subset data to include only target groups
  keep_samples <- group_data$group %in% target_groups
  exp_matrix <- exp_matrix[, keep_samples, drop = FALSE]
  group_data <- group_data[keep_samples, , drop = FALSE]
  group_data$group <- factor(group_data$group, levels = target_groups)
  
  message("\nGroup information after filtering for target groups:")
  print(table(group_data$group))
  
  # Debug information
  message("\nSample matching debug info:")
  message("Expression matrix dimensions: ", nrow(exp_matrix), " x ", ncol(exp_matrix))
  message("Number of group labels: ", nrow(group_data))
  message("\nFirst few sample names:")
  message("Expression matrix: ", paste(head(colnames(exp_matrix), 5), collapse = ", "))
  message("Group data: ", paste(head(rownames(group_data), 5), collapse = ", "))
  
  # Ensure sample names match
  if (!identical(colnames(exp_matrix), rownames(group_data))) {
    message("\nSample names don't match exactly. Attempting to find common samples...")
    common_samples <- intersect(colnames(exp_matrix), rownames(group_data))
    if (length(common_samples) == 0) {
      message("\nDetailed sample name comparison:")
      message("Expression matrix unique samples: ", length(unique(colnames(exp_matrix))))
      message("Group data unique samples: ", length(unique(rownames(group_data))))
      message("\nSample name examples that don't match:")
      exp_samples <- colnames(exp_matrix)
      group_samples <- rownames(group_data)
      message("In expression matrix but not in groups: ", 
              paste(head(setdiff(exp_samples, group_samples)), collapse = ", "))
      message("In groups but not in expression matrix: ", 
              paste(head(setdiff(group_samples, exp_samples)), collapse = ", "))
      stop("No matching samples found between expression matrix and group labels")
    }
    
    message(sprintf("\nFound %d matching samples", length(common_samples)))
    exp_matrix <- exp_matrix[, common_samples, drop = FALSE]
    group_data <- group_data[common_samples, , drop = FALSE]
  }
  
  message(sprintf("\nPerforming DEG analysis for %s", comparison_name))
  message("Number of samples: ", ncol(exp_matrix))
  message("Number of genes: ", nrow(exp_matrix))
  message("\nFinal group distribution:")
  print(table(group_data$group))
  
  # Check group sizes
  group_sizes <- table(group_data$group)
  if (any(group_sizes < min_samples_per_group)) {
    stop(sprintf("Some groups have fewer than %d samples:\n%s", 
                 min_samples_per_group, 
                 paste(names(group_sizes), group_sizes, sep = ": ", collapse = "\n")))
  }
  
  # Initialize results
  deg_results <- list()
  
  # Perform comparison between target groups
  group1 <- target_groups[1]  # cluster_high
  group2 <- target_groups[2]  # cluster_low
  
  # Get samples for each group
  group1_samples <- rownames(group_data)[group_data$group == group1]
  group2_samples <- rownames(group_data)[group_data$group == group2]
  
  message(sprintf("\nAnalyzing %s vs %s", group1, group2))
  message(sprintf("Number of samples in %s: %d", group1, length(group1_samples)))
  message(sprintf("Number of samples in %s: %d", group2, length(group2_samples)))
  
  # Subset data
  samples_subset <- c(group1_samples, group2_samples)
  exp_subset <- exp_matrix[, samples_subset, drop = FALSE]
  labels_subset <- factor(group_data$group[match(samples_subset, rownames(group_data))], 
                          levels = c(group1, group2))
  
  # Perform DEG analysis
  # 统一使用limma方法
  deg_result <- perform_limma_deg(exp_subset, labels_subset, group1, group2)
  
  comparison_key <- paste(group1, "vs", group2, sep = "_")
  deg_results[[comparison_key]] <- deg_result
  
  # Save results
  result_file <- file.path(output_dir, paste0(comparison_name, "_", comparison_key, "_DEG.xlsx"))
  write.xlsx(deg_result, result_file)
  
  # Create volcano plot
  create_volcano_plot(deg_result, comparison_key, 
                      file.path(output_dir, paste0(comparison_name, "_", comparison_key, "_volcano.pdf")))
  
  # Create MA plot
  create_ma_plot(deg_result, comparison_key,
                 file.path(output_dir, paste0(comparison_name, "_", comparison_key, "_MA.pdf")))
  
  # Print summary
  message(sprintf("\nComparison: %s", comparison_key))
  message("Total genes tested: ", nrow(deg_result))
  message("Up-regulated genes: ", sum(deg_result$significance == "Up"))
  message("Down-regulated genes: ", sum(deg_result$significance == "Down"))
  
  return(deg_results)
}

# Enhanced limma analysis for FPKM/TPM data
perform_limma_deg <- function(exp_matrix, group_labels, group1, group2) {
  # Create design matrix
  design <- model.matrix(~ 0 + group_labels)
  colnames(design) <- levels(group_labels)
  
  # Fit model
  fit <- lmFit(exp_matrix, design)
  
  # Create contrast
  contrast_matrix <- makeContrasts(
    contrasts = paste(group1, "-", group2, sep = ""),
    levels = design
  )
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE)  # Use trend=TRUE for better variance estimation
  
  # Get results
  results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  results$gene <- rownames(results)
  
  # Add mean expression for each group
  results$mean_group1 <- rowMeans(exp_matrix[, group_labels == group1])
  results$mean_group2 <- rowMeans(exp_matrix[, group_labels == group2])
  
  # Add significance categories
  results$significance <- "NS"
  results$significance[results$adj.P.Val < 0.05 & results$logFC > 0.5] <- "Up"
  results$significance[results$adj.P.Val < 0.05 & results$logFC < -0.5] <- "Down"
  
  return(results)
}

# Enhanced t-test analysis for Z-score data
perform_ttest_deg <- function(exp_matrix, group_labels, group1, group2) {
  # Get samples for each group
  group1_samples <- names(group_labels)[group_labels == group1]
  group2_samples <- names(group_labels)[group_labels == group2]
  
  # Initialize results
  results <- data.frame(
    gene = rownames(exp_matrix),
    logFC = numeric(nrow(exp_matrix)),
    mean_group1 = numeric(nrow(exp_matrix)),
    mean_group2 = numeric(nrow(exp_matrix)),
    t = numeric(nrow(exp_matrix)),
    P.Value = numeric(nrow(exp_matrix)),
    adj.P.Val = numeric(nrow(exp_matrix)),
    stringsAsFactors = FALSE
  )
  
  # Perform t-tests with variance estimation
  for (i in 1:nrow(exp_matrix)) {
    gene_exp <- exp_matrix[i, ]
    group1_exp <- gene_exp[group1_samples]
    group2_exp <- gene_exp[group2_samples]
    
    results$mean_group1[i] <- mean(group1_exp)
    results$mean_group2[i] <- mean(group2_exp)
    results$logFC[i] <- results$mean_group1[i] - results$mean_group2[i]
    
    if (length(group1_exp) > 1 && length(group2_exp) > 1) {
      # Use Welch's t-test (doesn't assume equal variances)
      t_test <- t.test(group1_exp, group2_exp, var.equal = FALSE)
      results$t[i] <- t_test$statistic
      results$P.Value[i] <- t_test$p.value
    } else {
      results$t[i] <- NA
      results$P.Value[i] <- NA
    }
  }
  
  # Adjust p-values
  results$adj.P.Val <- p.adjust(results$P.Value, method = "BH")
  
  # Add significance categories
  results$significance <- "NS"
  results$significance[results$adj.P.Val < 0.05 & results$logFC > 0.5] <- "Up"
  results$significance[results$adj.P.Val < 0.05 & results$logFC < -0.5] <- "Down"
  
  results <- results[order(results$P.Value), ]
  
  return(results)
}

# Enhanced volcano plot with signature genes highlighted
create_volcano_plot <- function(deg_result, comparison_name, output_file) {
  # Get SCLC signatures
  sclc_signatures <- get_default_sclc_signatures()
  all_signature_genes <- unique(unlist(sclc_signatures))
  
  # Filter out NA values
  deg_result <- deg_result[!is.na(deg_result$P.Value), ]
  
  # Calculate point size based on absolute log fold change
  deg_result$point_size <- 0.5 + abs(deg_result$logFC) * 0.2
  deg_result$point_size <- pmin(deg_result$point_size, 2)  # Cap maximum size
  
  # Add signature information
  deg_result$is_signature <- deg_result$gene %in% all_signature_genes
  
  # Get significant signature genes
  sig_signature_genes <- deg_result[deg_result$is_signature & 
                                      deg_result$significance != "NS", ]
  
  # Increase size for signature genes
  sig_signature_genes$point_size <- sig_signature_genes$point_size * 2
  
  #分别获取上调和下调的显著性签名基因
  # 按照-log10(P.Value)排序，确保最显著的基因优先标注
  # 分别获取上调和下调的显著性签名基因
  # up_sig_signature_genes <- sig_signature_genes[sig_signature_genes$logFC > 0, ]
  # down_sig_signature_genes <- sig_signature_genes[sig_signature_genes$logFC < 0, ]
  
  # 对上调和下调基因分别进行排序和分组
  up_sig_genes <- sig_signature_genes[sig_signature_genes$logFC > 0, ] %>%
    arrange(desc(-log10(P.Value))) %>%
    mutate(label_y_pos = rank(-log10(P.Value)),  # 根据显著性创建垂直排序
           label_x_pos = logFC + 0.2 * (label_y_pos %% 3))  # 水平轻微偏移
  
  down_sig_genes <- sig_signature_genes[sig_signature_genes$logFC < 0, ] %>%
    arrange(desc(-log10(P.Value))) %>%
    mutate(label_y_pos = rank(-log10(P.Value)),
           label_x_pos = logFC - 0.2 * (label_y_pos %% 3))
  
  # Create plot with enhanced aesthetics
  p <- ggplot(deg_result, aes(x = logFC, y = -log10(P.Value))) +
    # Add background points
    geom_point(data = subset(deg_result, !is_signature), 
               aes(color = significance, size = point_size), 
               alpha = 0.4) +  # Reduced alpha for background points
    # Add signature genes with different shape and larger size
    geom_point(data = sig_signature_genes, 
               aes(color = significance, size = point_size),
               shape = 17, alpha = 0.9) +  # Increased alpha for signature genes
    scale_color_manual(values = c("Up" = "#DC143C", "Down" = "#2E8B57", "NS" = "#ADB6CA")) +
    scale_size_identity() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    labs(title = paste("Volcano Plot:", comparison_name),
         subtitle = "Triangle points indicate signature genes",
         x = "Log2 Fold Change",
         y = "-log10(P-value)",
         color = "Significance") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # Add labels for all significant signature genes
  if (nrow(sig_signature_genes) > 0) {
    p <- p + geom_text_repel(
      data = up_sig_genes,
      aes(x = label_x_pos, y = -log10(P.Value), label = gene),
      size = 3,
      color = "black", 
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50",
      force = 20,  # 增加分散力
      force_pull = 0.5,  # 降低拉力
      xlim = c(0, NA),  # 限制在右侧
      direction = "both",  # 允许垂直和水平方向
      segment.alpha = 0.5,
      min.segment.length = 0.2,
      nudge_x = 0.2,  # 轻微右移
      seed = 42  # 保持随机性一致
    ) +
      geom_text_repel(
        data = down_sig_genes,
        aes(x = label_x_pos, y = -log10(P.Value), label = gene),
        size = 3,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "grey50", 
        force = 20,
        force_pull = 0.5,
        xlim = c(NA, 0),  # 限制在左侧
        direction = "both",
        segment.alpha = 0.5,
        min.segment.length = 0.2,
        nudge_x = -0.3,  # 轻微左移
        seed = 42
      )
  }
  
  ggsave(output_file, 
         plot = p, 
         width = 5, 
         height =4, 
         dpi = 300)
  
  return(p)
}

# Create MA plot
create_ma_plot <- function(deg_result, comparison_name, output_file) {
  # Calculate mean expression
  deg_result$mean_expr <- (deg_result$mean_group1 + deg_result$mean_group2) / 2
  
  # Create MA plot
  p <- ggplot(deg_result, aes(x = mean_expr, y = logFC, color = significance)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#DC143C", "Down" = "#2E8B57", "NS" = "#ADB6CA")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    labs(title = paste("MA Plot:", comparison_name),
         x = "Mean Expression",
         y = "Log2 Fold Change",
         color = "Significance") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # Add labels for top genes
  top_genes <- rbind(
    head(deg_result[deg_result$significance == "Up", ], 5),
    head(deg_result[deg_result$significance == "Down", ], 5)
  )
  
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.1,
      segment.color = "gray50",
      segment.alpha = 0.5
    )
  }
  
  # Save plot
  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# Enhanced pathway analysis
perform_pathway_analysis <- function(deg_results, output_dir = "./results/pathway/",
                                     comparison_name = "comparison",
                                     min_genes = 10,
                                     max_pathways = 30) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message(sprintf("Performing pathway analysis for %s", comparison_name))
  
  # Load MSigDB Hallmark gene sets
  hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  
  pathway_results <- list()
  
  for (comp_name in names(deg_results)) {
    deg_data <- deg_results[[comp_name]]
    
    # Get significant genes
    up_genes <- deg_data$gene[deg_data$significance == "Up"]
    down_genes <- deg_data$gene[deg_data$significance == "Down"]
    
    message(sprintf("\nProcessing comparison: %s", comp_name))
    message(sprintf("Up-regulated genes: %d", length(up_genes)))
    message(sprintf("Down-regulated genes: %d", length(down_genes)))
    
    enrichment_results <- list(up = list(), down = list())
    
    # Process up-regulated genes
    if (length(up_genes) >= min_genes) {
      message("\nAnalyzing up-regulated genes...")
      tryCatch({
        # Convert to ENTREZ IDs for GO and KEGG
        up_entrez <- bitr(up_genes, fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)$ENTREZID
        
        # GO analysis
        enrichment_results$up$GO <- enrichGO(
          gene = up_entrez,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        # KEGG analysis
        enrichment_results$up$KEGG <- enrichKEGG(
          gene = up_entrez,
          organism = "hsa",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        # Hallmark analysis
        enrichment_results$up$Hallmark <- enricher(
          gene = up_genes,
          TERM2GENE = hallmark_gene_sets,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        message("GO terms found: ", nrow(enrichment_results$up$GO@result))
        message("KEGG pathways found: ", nrow(enrichment_results$up$KEGG@result))
        message("Hallmark gene sets found: ", nrow(enrichment_results$up$Hallmark@result))
      }, error = function(e) {
        message("Error in up-regulated pathway analysis: ", e$message)
      })
    }
    
    # Process down-regulated genes
    if (length(down_genes) >= min_genes) {
      message("\nAnalyzing down-regulated genes...")
      tryCatch({
        # Convert to ENTREZ IDs for GO and KEGG
        down_entrez <- bitr(down_genes, fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)$ENTREZID
        
        # GO analysis
        enrichment_results$down$GO <- enrichGO(
          gene = down_entrez,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        # KEGG analysis
        enrichment_results$down$KEGG <- enrichKEGG(
          gene = down_entrez,
          organism = "hsa",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        # Hallmark analysis
        enrichment_results$down$Hallmark <- enricher(
          gene = down_genes,
          TERM2GENE = hallmark_gene_sets,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        message("GO terms found: ", nrow(enrichment_results$down$GO@result))
        message("KEGG pathways found: ", nrow(enrichment_results$down$KEGG@result))
        message("Hallmark gene sets found: ", nrow(enrichment_results$down$Hallmark@result))
      }, error = function(e) {
        message("Error in down-regulated pathway analysis: ", e$message)
      })
    }
    
    # Store results
    pathway_results[[comp_name]] <- enrichment_results
    
    # Create individual database plots and save results
    tryCatch({
      # Create plots for each database
      create_individual_database_plots(
        enrichment_results = enrichment_results,
        comp_name = comp_name,
        output_dir = output_dir,
        comparison_name = comparison_name
      )
      
      # Create combined pathway plot for upregulated pathways
      create_combined_pathway_plot(
        enrichment_results = enrichment_results,
        comp_name = comp_name,
        output_dir = output_dir,
        comparison_name = comparison_name
      )
      
      # Save detailed results
      save_pathway_results(enrichment_results, comp_name, output_dir, comparison_name)
      save_pathway_results_txt(enrichment_results, comp_name, output_dir, comparison_name)
    }, error = function(e) {
      message("Error in creating plots or saving results: ", e$message)
    })
  }
  
  return(pathway_results)
}

# Define key cancer pathways
get_key_cancer_pathways <- function() {
  list(
    signaling_pathways = c(
      "PI3K",           # PI3K-Akt pathway
      'AKT',
      "MAPK",                # MAPK pathway
      "Wnt",                 # Wnt pathway
      "JAK",          # JAK-STAT pathway
      'STAT',
      "Toll-like",         # Toll-like pathway
      "Ras",                # Ras pathway
      "Rap1",               # Rap1 pathway
      "Hedgehog",           # Hedgehog pathway
      "Hippo",              # Hippo pathway
      "Notch"               # Notch pathway
    )
  )
}

# Function to find key pathways in enrichment results
find_key_pathways <- function(enrichment_result) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    return(NULL)
  }
  
  key_patterns <- get_key_cancer_pathways()$signaling_pathways
  
  # Search for key pathways in the results
  key_pathways <- lapply(key_patterns, function(pattern) {
    matches <- grep(pattern, enrichment_result@result$Description, 
                    ignore.case = TRUE, value = FALSE)
    if (length(matches) > 0) {
      return(enrichment_result@result[matches, ])
    }
    return(NULL)
  })
  
  # Combine all found key pathways
  key_pathways <- do.call(rbind, key_pathways[!sapply(key_pathways, is.null)])
  
  if (!is.null(key_pathways)) {
    key_pathways <- unique(key_pathways)  # Remove any duplicates
  }
  
  return(key_pathways)
}

# Modified create_individual_database_plots function
create_individual_database_plots <- function(enrichment_results, comp_name, output_dir, comparison_name) {
  
  # Function to create plot for a specific database and direction
  create_db_plot <- function(result, database, direction, max_pathways = 30) {
    if (!is.null(result) && nrow(result@result) > 0) {
      # Get top pathways and key pathways
      top_results <- head(result@result, max_pathways)
      key_pathways <- find_key_pathways(result)
      
      # Combine top results with key pathways if direction is "up"
      if (direction == "Up" && !is.null(key_pathways)) {
        # Remove key pathways that are already in top results
        key_pathways <- key_pathways[!key_pathways$Description %in% top_results$Description, ]
        if (nrow(key_pathways) > 0) {
          plot_results <- rbind(top_results, key_pathways)
        } else {
          plot_results <- top_results
        }
      } else {
        plot_results <- top_results
      }
      
      # Calculate enrichment score
      gene_ratio <- as.numeric(sapply(strsplit(plot_results$GeneRatio, "/"), 
                                      function(x) as.numeric(x[1])/as.numeric(x[2])))
      score <- -log10(plot_results$p.adjust) * gene_ratio
      if(direction == "Down") score <- -score
      
      # Create data frame for plotting
      plot_data <- data.frame(
        Pathway = plot_results$Description,
        Score = score,
        P.Value = plot_results$p.adjust,
        Count = plot_results$Count,
        GeneRatio = gene_ratio,
        IsKeyPathway = plot_results$Description %in% key_pathways$Description,
        stringsAsFactors = FALSE
      )
      
      # Create plot
      p <- ggplot(plot_data, 
                  aes(x = Score, 
                      y = reorder(Pathway, Score),
                      size = Count)) +
        geom_point(aes(color = -log10(P.Value), 
                       shape = IsKeyPathway), 
                   alpha = 0.8) +
        scale_color_gradient(low = "blue", high = "red") +
        scale_size_continuous(range = c(3, 10)) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        labs(title = paste(database, "Pathways -", direction, "regulated"),
             subtitle = paste("Comparison:", comp_name),
             x = "Enrichment Score",
             y = "Pathway",
             size = "Gene Count",
             color = "-log10(adj.P.Value)",
             shape = "Key Pathway") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          axis.text.y = element_text(
            size = 14,
            face = "bold",
            angle = 0,  # 改为0度水平显示
            hjust = 1,  # 右对齐
            margin = margin(r = 5)  # 减少右边距
          ),
          axis.title = element_text(size = 16),
          legend.position = "right",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 10),  # 减少整体边距
          panel.spacing = unit(0.5, "lines")
        ) +
        guides(
          size = guide_legend(override.aes = list(size = 3)),  # 调整图例点的大小
          color = guide_legend(override.aes = list(size = 3))
        )
      
      return(list(plot = p, data = plot_data))
    }
    return(NULL)
  }
  
  # Process each database and direction
  databases <- c("GO", "KEGG", "Hallmark")
  directions <- c("up", "down")
  
  for (db in databases) {
    for (dir in directions) {
      if (!is.null(enrichment_results[[dir]]) && !is.null(enrichment_results[[dir]][[db]])) {
        tryCatch({
          result <- create_db_plot(enrichment_results[[dir]][[db]], db, 
                                   stringr::str_to_title(dir))
          
          if (!is.null(result)) {
            # Save plot
            plot_filename <- file.path(output_dir, 
                                       paste0(comparison_name, "_", comp_name, "_", db, "_", dir, ".pdf"))
            
            # Adjust plot height based on number of pathways
            plot_height <- max(10, min(nrow(result$data) * 0.8, 30))  # 增加高度系数
            aspect_ratio <- 1.618  # 黄金比例
            
            ggsave(plot_filename, 
                   plot = result$plot,
                   width = 16,  # 增加宽度
                   height = plot_height, 
                   dpi = 300,
                   limitsize = FALSE,
                   device = cairo_pdf)  # 使用Cairo PDF设备
            
            # 添加自动换行函数（在创建plot_data前）：
            wrap_labels <- function(labels, width = 50) {
              sapply(labels, function(x) {
                paste(strwrap(x, width = width), collapse = "\n")
              })
            }
            
            plot_data$Pathway <- wrap_labels(plot_data$Pathway)
            
            # Save pathway details to text file
            text_filename <- file.path(output_dir, 
                                       paste0(comparison_name, "_", comp_name, "_", db, "_", dir, "_pathways.txt"))
            
            con <- file(text_filename, "w")
            writeLines(paste("Enriched pathways in", db, "-", dir, "regulated genes"), con)
            writeLines(paste("Comparison:", comp_name), con)
            writeLines("", con)
            
            # First write key pathways if any
            key_pathways <- result$data[result$data$IsKeyPathway, ]
            if (nrow(key_pathways) > 0) {
              writeLines("\nKey Cancer Signaling Pathways:", con)
              for(i in 1:nrow(key_pathways)) {
                writeLines(paste0(i, ". ", key_pathways$Pathway[i]), con)
                writeLines(paste0("   Adjusted p-value: ", format(key_pathways$P.Value[i], scientific = TRUE)), con)
                writeLines(paste0("   Gene count: ", key_pathways$Count[i]), con)
                writeLines("", con)
              }
            }
            
            # Then write all pathways
            writeLines("\nAll Enriched Pathways:", con)
            for(i in 1:nrow(result$data)) {
              writeLines(paste0(i, ". ", result$data$Pathway[i]), con)
              writeLines(paste0("   Adjusted p-value: ", format(result$data$P.Value[i], scientific = TRUE)), con)
              writeLines(paste0("   Gene count: ", result$data$Count[i]), con)
              writeLines(paste0("   Key pathway: ", result$data$IsKeyPathway[i]), con)
              writeLines("", con)
            }
            
            close(con)
            message(sprintf("Saved %s pathway details for %s direction to: %s", db, dir, text_filename))
            
            # Save plot data to CSV
            csv_filename <- file.path(output_dir, 
                                      paste0(comparison_name, "_", comp_name, "_", db, "_", dir, "_data.csv"))
            write.csv(result$data, csv_filename, row.names = FALSE)
            message(sprintf("Saved %s data for %s direction to: %s", db, dir, csv_filename))
          }
        }, error = function(e) {
          message(sprintf("Error processing %s %s: %s", db, dir, e$message))
        })
      }
    }
  }
}

# Save pathway results to txt files
save_pathway_results_txt <- function(enrichment_results, comp_name, output_dir, comparison_name) {
  
  # Create main summary file
  summary_file <- file.path(output_dir, paste0(comparison_name, "_", comp_name, "_pathway_summary.txt"))
  con <- file(summary_file, "w")
  
  writeLines(paste("Pathway Analysis Summary for", comp_name), con)
  writeLines(paste(rep("=", 60), collapse = ""), con)
  writeLines(paste("Generated on:", Sys.time()), con)
  writeLines("", con)
  
  # Function to process and write results for each database
  write_database_results <- function(results, direction, database) {
    if (!is.null(results) && nrow(results@result) > 0) {
      writeLines(paste("\n", database, "Database -", stringr::str_to_title(direction), "regulated pathways:"), con)
      writeLines(paste(rep("-", 50), collapse = ""), con)
      writeLines(paste("Total pathways found:", nrow(results@result)), con)
      writeLines("", con)
      
      # Write all pathways with details
      writeLines("Enriched pathways:", con)
      writeLines("", con)
      for (i in 1:nrow(results@result)) {
        writeLines(paste0(i, ". ", results@result$Description[i]), con)
        writeLines(paste0("   P.adjust: ", formatC(results@result$p.adjust[i], format = "e", digits = 2)), con)
        writeLines(paste0("   Gene Ratio: ", results@result$GeneRatio[i]), con)
        writeLines(paste0("   Gene Count: ", results@result$Count[i]), con)
        if (!is.null(results@result$geneID)) {
          genes <- strsplit(results@result$geneID[i], "/")[[1]]
          writeLines(paste0("   Genes: ", paste(genes, collapse = ", ")), con)
        }
        writeLines("", con)
      }
      
      return(TRUE)
    }
    return(FALSE)
  }
  
  # Process all databases and directions
  databases <- c("GO", "KEGG", "Hallmark")
  directions <- c("up", "down")
  
  total_pathways <- 0
  
  for (dir in directions) {
    if (!is.null(enrichment_results[[dir]])) {
      for (db in databases) {
        if (!is.null(enrichment_results[[dir]][[db]])) {
          found <- write_database_results(enrichment_results[[dir]][[db]], dir, db)
          if (found) {
            total_pathways <- total_pathways + nrow(enrichment_results[[dir]][[db]]@result)
          }
        }
      }
    }
  }
  
  # Overall summary
  writeLines("\n\nOVERALL SUMMARY:", con)
  writeLines(paste(rep("=", 30), collapse = ""), con)
  writeLines(paste("Total enriched pathways found:", total_pathways), con)
  
  close(con)
  
  message(sprintf("Detailed pathway analysis saved to: %s", summary_file))
  
  # Create individual database files
  for (dir in directions) {
    if (!is.null(enrichment_results[[dir]])) {
      for (db in databases) {
        if (!is.null(enrichment_results[[dir]][[db]]) && 
            nrow(enrichment_results[[dir]][[db]]@result) > 0) {
          
          db_file <- file.path(output_dir, 
                               paste0(comparison_name, "_", comp_name, "_", db, "_", dir, "_detailed.txt"))
          db_con <- file(db_file, "w")
          
          writeLines(paste(db, "Database -", stringr::str_to_title(dir), "regulated pathways"), db_con)
          writeLines(paste(rep("=", 60), collapse = ""), db_con)
          writeLines(paste("Comparison:", comp_name), db_con)
          writeLines(paste("Generated on:", Sys.time()), db_con)
          writeLines("", db_con)
          
          results <- enrichment_results[[dir]][[db]]@result
          writeLines(paste("Total pathways:", nrow(results)), db_con)
          writeLines("", db_con)
          
          # Write all pathways with details
          for (i in 1:nrow(results)) {
            writeLines(paste0("Pathway ", i, ": ", results$Description[i]), db_con)
            writeLines(paste0("P-value: ", formatC(results$pvalue[i], format = "e", digits = 3)), db_con)
            writeLines(paste0("Adjusted P-value: ", formatC(results$p.adjust[i], format = "e", digits = 3)), db_con)
            writeLines(paste0("Gene Ratio: ", results$GeneRatio[i]), db_con)
            writeLines(paste0("Background Ratio: ", results$BgRatio[i]), db_con)
            writeLines(paste0("Gene Count: ", results$Count[i]), db_con)
            if (!is.null(results$geneID)) {
              writeLines(paste0("Genes: ", results$geneID[i]), db_con)
            }
            writeLines("", db_con)
          }
          
          close(db_con)
        }
      }
    }
  }
}

# Save pathway results to Excel files
save_pathway_results <- function(enrichment_results, comp_name, output_dir, comparison_name) {
  # Create Excel workbook
  wb <- createWorkbook()
  
  # Function to add worksheet if results exist
  add_worksheet_if_exists <- function(wb, name, result) {
    if (!is.null(result) && nrow(result@result) > 0) {
      addWorksheet(wb, name)
      writeData(wb, name, result@result)
      return(TRUE)
    }
    return(FALSE)
  }
  
  # Process up-regulated pathways
  if (!is.null(enrichment_results$up)) {
    # GO results
    if (!is.null(enrichment_results$up$GO)) {
      add_worksheet_if_exists(wb, "GO_up", enrichment_results$up$GO)
    }
    
    # KEGG results
    if (!is.null(enrichment_results$up$KEGG)) {
      add_worksheet_if_exists(wb, "KEGG_up", enrichment_results$up$KEGG)
    }
    
    # Hallmark results
    if (!is.null(enrichment_results$up$Hallmark)) {
      add_worksheet_if_exists(wb, "Hallmark_up", enrichment_results$up$Hallmark)
    }
  }
  
  # Process down-regulated pathways
  if (!is.null(enrichment_results$down)) {
    # GO results
    if (!is.null(enrichment_results$down$GO)) {
      add_worksheet_if_exists(wb, "GO_down", enrichment_results$down$GO)
    }
    
    # KEGG results
    if (!is.null(enrichment_results$down$KEGG)) {
      add_worksheet_if_exists(wb, "KEGG_down", enrichment_results$down$KEGG)
    }
    
    # Hallmark results
    if (!is.null(enrichment_results$down$Hallmark)) {
      add_worksheet_if_exists(wb, "Hallmark_down", enrichment_results$down$Hallmark)
    }
  }
  
  # Save workbook if it has any worksheets
  if (length(names(wb)) > 0) {
    saveWorkbook(wb, 
                 file.path(output_dir, paste0(comparison_name, "_", comp_name, "_pathway.xlsx")),
                 overwrite = TRUE)
    message(sprintf("Pathway results saved to: %s", 
                    file.path(output_dir, paste0(comparison_name, "_", comp_name, "_pathway.xlsx"))))
  } else {
    message("No pathway results to save")
  }
}

# Create combined pathway plot from all databases
create_combined_pathway_plot <- function(enrichment_results, comp_name, output_dir, comparison_name) {
  if (is.null(enrichment_results$up)) {
    message("No up-regulated pathway results found")
    return(NULL)
  }
  
  # Combine results from all databases
  combined_data <- list()
  
  # Process GO results
  if (!is.null(enrichment_results$up$GO) && nrow(enrichment_results$up$GO@result) > 0) {
    go_data <- head(enrichment_results$up$GO@result, 5)  # Take top 10
    go_data$Database <- "GO"
    combined_data$GO <- go_data
  }
  
  # Process KEGG results
  if (!is.null(enrichment_results$up$KEGG) && nrow(enrichment_results$up$KEGG@result) > 0) {
    kegg_data <- head(enrichment_results$up$KEGG@result, 5)  # Take top 10
    kegg_data$Database <- "KEGG"
    combined_data$KEGG <- kegg_data
  }
  
  # Process Hallmark results
  if (!is.null(enrichment_results$up$Hallmark) && nrow(enrichment_results$up$Hallmark@result) > 0) {
    hallmark_data <- head(enrichment_results$up$Hallmark@result, 5)  # Take top 10
    hallmark_data$Database <- "Hallmark"
    combined_data$Hallmark <- hallmark_data
  }
  
  if (length(combined_data) == 0) {
    message("No pathway results found in any database")
    return(NULL)
  }
  
  # Combine all results
  all_results <- do.call(rbind, combined_data)
  
  # Calculate enrichment score
  gene_ratio <- as.numeric(sapply(strsplit(all_results$GeneRatio, "/"), 
                                  function(x) as.numeric(x[1])/as.numeric(x[2])))
  score <- -log10(all_results$p.adjust) * gene_ratio
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Pathway = all_results$Description,
    Score = score,
    P.Value = all_results$p.adjust,
    Count = all_results$Count,
    GeneRatio = gene_ratio,
    Database = all_results$Database,
    stringsAsFactors = FALSE
  )
  
  # Find key pathways
  key_patterns <- get_key_cancer_pathways()$signaling_pathways
  plot_data$IsKeyPathway <- FALSE
  for (pattern in key_patterns) {
    matches <- grep(pattern, plot_data$Pathway, ignore.case = TRUE)
    plot_data$IsKeyPathway[matches] <- TRUE
  }
  
  # Get top 10 pathways from each database
  plot_data_final <- plot_data %>%
    group_by(Database) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Add key pathways that are not in top 10s
  key_pathways <- plot_data[plot_data$IsKeyPathway & !(plot_data$Pathway %in% plot_data_final$Pathway), ]
  if(nrow(key_pathways) > 0) {
    plot_data_final <- rbind(plot_data_final, key_pathways)
  }
  
  # Create the combined plot
  p <- ggplot(plot_data_final, 
              aes(x = Score, 
                  y = reorder(Pathway, Score),
                  size = Count)) +
    geom_point(aes(color = Database,
                   shape = IsKeyPathway),
               alpha = 0.8) +
    scale_color_manual(values = c("GO" = "#1f77b4", 
                                  "KEGG" = "#2ca02c", 
                                  "Hallmark" = "#ff7f0e")) +
    scale_size_continuous(range = c(3, 10)) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18),
                       labels = c("Regular Pathway", "Key Cancer Pathway")) +
    labs(title = "Combined Upregulated Pathways",
         subtitle = paste("Comparison:", comp_name, "\nShowing top 10 pathways from each database and enriched key cancer pathways"),
         x = "Enrichment Score",
         y = "Pathway",
         size = "Gene Count",
         shape = "Pathway Type") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(
        size = 14,
        face = "bold",
        lineheight = 0.3,  # 从0.7→0.5进一步缩小行间距
        margin = margin(r = 2)  # 减少右边距
      ),
      axis.title = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.margin = margin(3, 3, 3, 3),  # 边距从5→3
      panel.spacing = unit(0.1, "lines")  # 面板间距从0.2→0.1
    )
  
  # Save plot
  plot_filename <- file.path(output_dir, 
                             paste0(comparison_name, "_", comp_name, "_combined_up_pathways.pdf"))
  
  # Adjust plot height based on number of pathways
  plot_height <- max(6, min(nrow(plot_data_final) * 0.25, 15))  # 高度系数从0.4→0.25
  ggsave(plot_filename, 
         plot = p, 
         width = 8,  # 宽度从10→8
         height = plot_height, 
         dpi = 300,
         limitsize = FALSE)
  
  # 在自动换行函数中调整参数：
  wrap_labels <- function(labels, width = 35) {  # 从40→35
    sapply(labels, function(x) {
      paste(strwrap(x, width = width), collapse = "\n")
    })
  }
  
  plot_data_final$Pathway <- wrap_labels(plot_data_final$Pathway)
  
  message(sprintf("Saved combined pathway plot to: %s", plot_filename))
  
  # Save detailed results to text file
  text_filename <- file.path(output_dir, 
                             paste0(comparison_name, "_", comp_name, "_combined_up_pathways.txt"))
  
  con <- file(text_filename, "w")
  writeLines("Combined Upregulated Pathways Analysis", con)
  writeLines(paste("Comparison:", comp_name), con)
  writeLines("", con)
  
  # Write key pathways
  key_pathways <- plot_data_final[plot_data_final$IsKeyPathway, ]
  if (nrow(key_pathways) > 0) {
    writeLines("\nEnriched Key Cancer Signaling Pathways:", con)
    for(i in 1:nrow(key_pathways)) {
      writeLines(paste0(i, ". ", key_pathways$Pathway[i]), con)
      writeLines(paste0("   Database: ", key_pathways$Database[i]), con)
      writeLines(paste0("   Adjusted p-value: ", format(key_pathways$P.Value[i], scientific = TRUE)), con)
      writeLines(paste0("   Gene count: ", key_pathways$Count[i]), con)
      writeLines(paste0("   Enrichment score: ", format(key_pathways$Score[i], digits = 3)), con)
      writeLines("", con)
    }
  }
  
  # Write all pathways by database
  writeLines("\nTop 10 Pathways by Database:", con)
  for(db in unique(plot_data_final$Database)) {
    writeLines(paste0("\n", db, " Pathways:"), con)
    db_pathways <- plot_data_final[plot_data_final$Database == db & !plot_data_final$IsKeyPathway, ]
    for(i in 1:nrow(db_pathways)) {
      writeLines(paste0(i, ". ", db_pathways$Pathway[i]), con)
      writeLines(paste0("   Adjusted p-value: ", format(db_pathways$P.Value[i], scientific = TRUE)), con)
      writeLines(paste0("   Gene count: ", db_pathways$Count[i]), con)
      writeLines(paste0("   Enrichment score: ", format(db_pathways$Score[i], digits = 3)), con)
      writeLines("", con)
    }
  }
  
  close(con)
  message(sprintf("Saved combined pathway details to: %s", text_filename))
  
  # Save data to CSV
  csv_filename <- file.path(output_dir, 
                            paste0(comparison_name, "_", comp_name, "_combined_up_pathways.csv"))
  write.csv(plot_data_final, csv_filename, row.names = FALSE)
  message(sprintf("Saved combined pathway data to: %s", csv_filename))
  
  return(p)
}

# Main analysis wrapper
run_deg_pathway_analysis <- function(exp_file = NULL, 
                                     group_file = NULL,
                                     exp_matrix = NULL,
                                     group_labels = NULL,
                                     data_type = "FPKM", 
                                     dataset_name = "analysis",
                                     output_base_dir = "./results/",
                                     min_samples_per_group = 3,
                                     lfc_threshold = 0.5,
                                     pval_threshold = 0.05,
                                     target_groups = c("cluster_high", "cluster_low")) {
  
  # Load required packages
  load_required_packages()
  
  # Create output directories
  output_dir <- file.path(output_base_dir, dataset_name)
  deg_output_dir <- file.path(output_dir, "DEG")
  pathway_output_dir <- file.path(output_dir, "pathway")
  
  message(sprintf("\n=== Starting analysis for dataset: %s ===\n", dataset_name))
  
  # Handle both old and new parameter formats
  if (!is.null(exp_matrix) && !is.null(group_labels)) {
    # Using old parameter format
    message("Using directly provided expression matrix and group labels")
    
    # Convert group_labels to proper format if needed
    if (is.data.frame(group_labels)) {
      if ("cluster" %in% colnames(group_labels)) {
        group_data <- data.frame(
          group = group_labels$cluster,
          row.names = rownames(group_labels),
          stringsAsFactors = FALSE
        )
      } else if ("group" %in% colnames(group_labels)) {
        group_data <- data.frame(
          group = group_labels$group,
          row.names = rownames(group_labels),
          stringsAsFactors = FALSE
        )
      } else {
        group_data <- data.frame(
          group = group_labels[,1],
          row.names = rownames(group_labels),
          stringsAsFactors = FALSE
        )
      }
    } else {
      group_data <- data.frame(
        group = group_labels,
        row.names = names(group_labels),
        stringsAsFactors = FALSE
      )
    }
    
    input_data <- list(
      exp_matrix = as.matrix(exp_matrix),
      group_data = group_data
    )
  } else if (!is.null(exp_file) || !is.null(group_file)) {
    # Using new parameter format
    message("Reading data from files")
    input_data <- read_input_data(exp_file, group_file)
  } else {
    stop("Either provide exp_matrix and group_labels directly, or provide exp_file and group_file paths")
  }
  
  # Extract data
  exp_matrix <- input_data$exp_matrix
  group_data <- input_data$group_data
  
  # Debug information
  message("\nInitial data check:")
  message("Expression matrix dimensions: ", nrow(exp_matrix), " x ", ncol(exp_matrix))
  message("Group data dimensions: ", nrow(group_data), " x ", ncol(group_data))
  message("\nFirst few sample names:")
  message("Expression matrix: ", paste(head(colnames(exp_matrix), 5), collapse = ", "))
  message("Group data: ", paste(head(rownames(group_data), 5), collapse = ", "))
  message("\nGroup distribution:")
  print(table(group_data$group))
  
  # Match samples
  matched_data <- match_samples(exp_matrix, group_data)
  exp_matrix <- matched_data$exp_matrix
  group_data <- matched_data$group_data
  
  # Preprocess expression data
  processed_exp <- preprocess_expression_data(exp_matrix, data_type)
  
  # Perform DEG analysis
  deg_results <- perform_deg_analysis(
    exp_matrix = processed_exp,
    group_labels = group_data,
    data_type = data_type,
    output_dir = deg_output_dir,
    comparison_name = dataset_name,
    min_samples_per_group = min_samples_per_group,
    lfc_threshold = lfc_threshold,
    pval_threshold = pval_threshold,
    target_groups = target_groups
  )
  
  # Perform pathway analysis
  pathway_results <- perform_pathway_analysis(
    deg_results = deg_results,
    output_dir = pathway_output_dir,
    comparison_name = dataset_name
  )
  
  # Create summary report
  create_analysis_summary(
    dataset_name = dataset_name,
    deg_results = deg_results,
    pathway_results = pathway_results,
    output_dir = output_dir
  )
  
  message(sprintf("\n=== Analysis complete for dataset: %s ===\n", dataset_name))
  
  return(list(
    DEG = deg_results,
    pathway = pathway_results,
    data = list(
      exp_matrix = processed_exp,
      group_data = group_data
    )
  ))
}

# Create analysis summary
create_analysis_summary <- function(dataset_name, deg_results, pathway_results, output_dir) {
  summary_file <- file.path(output_dir, paste0(dataset_name, "_analysis_summary.txt"))
  
  con <- file(summary_file, "w")
  
  writeLines(sprintf("Analysis Summary for %s\n", dataset_name), con)
  writeLines(paste(rep("=", 50), collapse = ""), con)
  
  # DEG summary
  writeLines("\nDifferential Expression Analysis:", con)
  writeLines(paste(rep("-", 30), collapse = ""), con)
  
  for (comp_name in names(deg_results)) {
    deg_data <- deg_results[[comp_name]]
    up_genes <- sum(deg_data$significance == "Up")
    down_genes <- sum(deg_data$significance == "Down")
    
    writeLines(sprintf("\nComparison: %s", comp_name), con)
    writeLines(sprintf("Total genes tested: %d", nrow(deg_data)), con)
    writeLines(sprintf("Up-regulated genes: %d", up_genes), con)
    writeLines(sprintf("Down-regulated genes: %d", down_genes), con)
  }
  
  # Pathway analysis summary
  writeLines("\n\nPathway Enrichment Analysis:", con)
  writeLines(paste(rep("-", 30), collapse = ""), con)
  
  for (comp_name in names(pathway_results)) {
    pathway_data <- pathway_results[[comp_name]]
    
    writeLines(sprintf("\nComparison: %s", comp_name), con)
    
    # GO BP (up-regulated)
    if (!is.null(pathway_data$up$GO) && nrow(pathway_data$up$GO@result) > 0) {
      writeLines(sprintf("GO BP (up-regulated) terms: %d", 
                         nrow(pathway_data$up$GO@result)), con)
    }
    
    # GO BP (down-regulated)
    if (!is.null(pathway_data$down$GO) && nrow(pathway_data$down$GO@result) > 0) {
      writeLines(sprintf("GO BP (down-regulated) terms: %d", 
                         nrow(pathway_data$down$GO@result)), con)
    }
    
    # KEGG
    if (!is.null(pathway_data$up$KEGG) && nrow(pathway_data$up$KEGG@result) > 0) {
      writeLines(sprintf("KEGG pathways: %d", 
                         nrow(pathway_data$up$KEGG@result)), con)
    }
    
    # Hallmark
    if (!is.null(pathway_data$up$Hallmark) && nrow(pathway_data$up$Hallmark@result) > 0) {
      writeLines(sprintf("Hallmark gene sets: %d", 
                         nrow(pathway_data$up$Hallmark@result)), con)
    }
  }
  
  close(con)
  message(sprintf("Analysis summary saved to: %s", summary_file))
}

# Combine pathway results from all databases into a single file
combine_pathway_results <- function(pathway_results, output_dir, comparison_name) {
  message("Combining pathway results from all databases...")
  
  for (comp_name in names(pathway_results)) {
    # Initialize combined results for up and down regulated pathways
    combined_up <- data.frame()
    combined_down <- data.frame()
    
    # Process up-regulated pathways
    if (!is.null(pathway_results[[comp_name]]$up)) {
      # Add GO results
      if (!is.null(pathway_results[[comp_name]]$up$GO) && 
          nrow(pathway_results[[comp_name]]$up$GO@result) > 0) {
        go_data <- pathway_results[[comp_name]]$up$GO@result
        go_data$Database <- "GO"
        combined_up <- rbind(combined_up, go_data)
      }
      
      # Add KEGG results
      if (!is.null(pathway_results[[comp_name]]$up$KEGG) && 
          nrow(pathway_results[[comp_name]]$up$KEGG@result) > 0) {
        kegg_data <- pathway_results[[comp_name]]$up$KEGG@result
        kegg_data$Database <- "KEGG"
        combined_up <- rbind(combined_up, kegg_data)
      }
      
      # Add Hallmark results
      if (!is.null(pathway_results[[comp_name]]$up$Hallmark) && 
          nrow(pathway_results[[comp_name]]$up$Hallmark@result) > 0) {
        hallmark_data <- pathway_results[[comp_name]]$up$Hallmark@result
        hallmark_data$Database <- "Hallmark"
        combined_up <- rbind(combined_up, hallmark_data)
      }
    }
    
    # Process down-regulated pathways
    if (!is.null(pathway_results[[comp_name]]$down)) {
      # Add GO results
      if (!is.null(pathway_results[[comp_name]]$down$GO) && 
          nrow(pathway_results[[comp_name]]$down$GO@result) > 0) {
        go_data <- pathway_results[[comp_name]]$down$GO@result
        go_data$Database <- "GO"
        combined_down <- rbind(combined_down, go_data)
      }
      
      # Add KEGG results
      if (!is.null(pathway_results[[comp_name]]$down$KEGG) && 
          nrow(pathway_results[[comp_name]]$down$KEGG@result) > 0) {
        kegg_data <- pathway_results[[comp_name]]$down$KEGG@result
        kegg_data$Database <- "KEGG"
        combined_down <- rbind(combined_down, kegg_data)
      }
      
      # Add Hallmark results
      if (!is.null(pathway_results[[comp_name]]$down$Hallmark) && 
          nrow(pathway_results[[comp_name]]$down$Hallmark@result) > 0) {
        hallmark_data <- pathway_results[[comp_name]]$down$Hallmark@result
        hallmark_data$Database <- "Hallmark"
        combined_down <- rbind(combined_down, hallmark_data)
      }
    }
    
    # Save combined results if they exist
    if (nrow(combined_up) > 0) {
      # Add direction column
      combined_up$Direction <- "Up"
      
      # Calculate enrichment score
      gene_ratio <- as.numeric(sapply(strsplit(combined_up$GeneRatio, "/"), 
                                      function(x) as.numeric(x[1])/as.numeric(x[2])))
      combined_up$EnrichmentScore <- -log10(combined_up$p.adjust) * gene_ratio
      
      # Sort by enrichment score
      combined_up <- combined_up[order(-combined_up$EnrichmentScore), ]
      
      # Save to Excel
      xlsx_file <- file.path(output_dir, 
                             paste0(comparison_name, "_", comp_name, "_all_up_pathways.xlsx"))
      write.xlsx(combined_up, xlsx_file)
      message(sprintf("Combined up-regulated pathways saved to: %s", xlsx_file))
      
      # Save to CSV
      csv_file <- file.path(output_dir, 
                            paste0(comparison_name, "_", comp_name, "_all_up_pathways.csv"))
      write.csv(combined_up, csv_file, row.names = FALSE)
    }
    
    if (nrow(combined_down) > 0) {
      # Add direction column
      combined_down$Direction <- "Down"
      
      # Calculate enrichment score
      gene_ratio <- as.numeric(sapply(strsplit(combined_down$GeneRatio, "/"), 
                                      function(x) as.numeric(x[1])/as.numeric(x[2])))
      combined_down$EnrichmentScore <- -log10(combined_down$p.adjust) * gene_ratio
      
      # Sort by enrichment score
      combined_down <- combined_down[order(-combined_down$EnrichmentScore), ]
      
      # Save to Excel
      xlsx_file <- file.path(output_dir, 
                             paste0(comparison_name, "_", comp_name, "_all_down_pathways.xlsx"))
      write.xlsx(combined_down, xlsx_file)
      message(sprintf("Combined down-regulated pathways saved to: %s", xlsx_file))
      
      # Save to CSV
      csv_file <- file.path(output_dir, 
                            paste0(comparison_name, "_", comp_name, "_all_down_pathways.csv"))
      write.csv(combined_down, csv_file, row.names = FALSE)
    }
    
    # Combine both up and down regulated pathways
    if (nrow(combined_up) > 0 || nrow(combined_down) > 0) {
      all_pathways <- rbind(combined_up, combined_down)
      
      # Save combined results
      xlsx_file <- file.path(output_dir, 
                             paste0(comparison_name, "_", comp_name, "_all_pathways.xlsx"))
      write.xlsx(all_pathways, xlsx_file)
      message(sprintf("All pathways combined and saved to: %s", xlsx_file))
      
      # Create comprehensive visualization for all pathways
      create_comprehensive_pathway_plot(all_pathways, comp_name, output_dir, comparison_name)
    }
  }
}

# Create comprehensive visualization for all pathways
create_comprehensive_pathway_plot <- function(all_pathways, comp_name, output_dir, comparison_name) {
  # Select top pathways by database and direction
  top_pathways <- all_pathways %>%
    group_by(Database, Direction) %>%
    slice_max(order_by = EnrichmentScore, n = 10) %>%
    ungroup()
  
  # Modify EnrichmentScore for down-regulated pathways to be negative
  top_pathways$PlotScore <- top_pathways$EnrichmentScore
  top_pathways$PlotScore[top_pathways$Direction == "Down"] <- 
    -top_pathways$EnrichmentScore[top_pathways$Direction == "Down"]
  
  # Create plot
  p <- ggplot(top_pathways, 
              aes(x = PlotScore, 
                  y = reorder(Description, PlotScore),
                  size = Count,
                  color = Database,
                  shape = Direction)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("GO" = "#1f77b4", 
                                  "KEGG" = "#2ca02c", 
                                  "Hallmark" = "#ff7f0e")) +
    scale_shape_manual(values = c("Up" = 16, "Down" = 17)) +
    scale_size_continuous(range = c(3, 10)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Top Enriched Pathways",
         subtitle = paste("Comparison:", comp_name),
         x = "Enrichment Score (negative = downregulated)",
         y = "Pathway",
         size = "Gene Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # Save plot
  plot_filename <- file.path(output_dir, 
                             paste0(comparison_name, "_", comp_name, "_comprehensive_pathways.pdf"))
  
  # Adjust plot height based on number of pathways
  plot_height <- max(8, min(nrow(top_pathways) * 0.3, 20))
  ggsave(plot_filename, 
         plot = p, 
         width = 12, 
         height = plot_height, 
         dpi = 300,
         limitsize = FALSE)
  
  message(sprintf("Comprehensive pathway plot saved to: %s", plot_filename))
}


# Example usage:
# results <- run_deg_pathway_analysis(
#   exp_matrix = exp_tcga,
#   group_labels = group_labels_tcga,
#   data_type = "FPKM",
#   dataset_name = "TCGA",
#   output_base_dir = "./results_deg_pathway/",
#   target_groups = c("cluster_high", "cluster_low")
# ) 