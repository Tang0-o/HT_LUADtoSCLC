# ================================================================================
# 多队列差异基因分析与风险预测模型构建（优化版 - 基于TF Regulon）
# Multi-cohort Differential Gene Expression Analysis and Risk Prediction Model (Enhanced - TF Regulon Based)
# 
# 主要改进：
# 1. 使用TF regulon基因集替代DEG交集分析
# 2. 增加表达过滤：去除在少于30个样本中有表达的基因
# 3. 加强Lasso筛选：使用lambda.1se规则和系数阈值筛选，获得更少更特异的特征
# ================================================================================

# 加载必需的R包
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(glmnet)
library(pROC)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(rms)
library(nomogramEx)
library(survival)
library(survminer)
library(corrplot)
library(gridExtra)
library(cowplot)

# 定义风险评分标准化函数
normalize_risk_scores <- function(scores) {
  # 首先处理无限值和NA值
  scores_clean <- scores
  scores_clean[is.infinite(scores_clean) | is.na(scores_clean)] <- NA
  
  # 检查是否还有有效值
  if(all(is.na(scores_clean))) {
    stop("所有风险评分都是无效值")
  }
  
  # 计算有效值的范围
  valid_scores <- scores_clean[!is.na(scores_clean)]
  min_score <- min(valid_scores)
  
  # 标准化到正值
  if(min_score < 0) {
    normalized <- valid_scores - min_score + 1
  } else {
    normalized <- valid_scores + 1
  }
  
  # 将标准化后的值放回原始位置
  scores_clean[!is.na(scores_clean)] <- normalized
  
  # 用最小值替换NA
  scores_clean[is.na(scores_clean)] <- min(normalized)
  
  return(scores_clean)
}

# 加载GMT读取函数
source('/home/data/tmh_project/SCLC/source_global/write_read_gmt.R')

# 设置工作目录和路径
base_dir <- "/home/data/tmh_project/SCLC"
tcga_cluster_file <- file.path(base_dir, "Fig3_multicohort_v1/1_TCGA/ssgsea/TCGA_hierarchical_clusters.csv")
tcga_expr_file <- file.path(base_dir, "Fig3_multicohort_v1/0_Fig3_cohort_data/1_TCGA_LUAD/LUAD_fpkm_new.csv")
oak_poplar_expr_file <- file.path(base_dir, "Fig3_multicohort_v1/0_Fig3_cohort_data/5_OAK_POPLAR/OAK_POPLAR_NOS_tpm.csv")
oak_poplar_cluster_file <- file.path(base_dir, "Fig3_multicohort_v1/3_predictions/OAK_POPLAR/ssgsea/hierarchical/knn/OAK_POPLAR_hierarchical_clusters.csv")
nt_expr_file <- file.path(base_dir, "Fig2_NTbulk_human/1_NTbulkdata_E-MTAB-10399/GeneFPKM_new.csv")
nt_metadata_file <- file.path(base_dir, "Fig2_NTbulk_human/1_NTbulkdata_E-MTAB-10399/metadata_all_48samples.csv")
regulon_file <- file.path(base_dir, "Fig2_NTbulk_human/regulon13_human_v1.gmt")

# 创建输出目录
output_dir <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/1_risk_prediction"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 定义配色方案
color_palette <- c(
  "NT_bulk" = "#FF8C00",  # DarkOrange
  "TCGA"    = "#3F51B5",  # Indigo
  "OAK"     = "#ED0000",  # Red
  "POPLAR"  = "#42B540",  # Green
  "ONCOSG"  = "#0099B4"   # Cyan
)

cluster_colors <- c(
  "cluster_low" = "#00468B",
  "cluster_medium" = "#ADB6CA",
  "cluster_high" = "#FDAF91"
)

# ================================================================================
# 步骤1: 读取TCGA样本分组信息
# ================================================================================
cat("步骤1: 读取TCGA样本分组信息...\n")

tcga_clusters <- read.csv(tcga_cluster_file, row.names = 1, stringsAsFactors = FALSE)
cat("TCGA样本分组统计:\n")
print(table(tcga_clusters$cluster))

# ================================================================================
# 步骤2: 读取TF regulon基因集并整理
# ================================================================================
cat("\n步骤2: 读取TF regulon基因集并整理...\n")

# 读取regulon基因集
regulon.list <- read_gmt_file(regulon_file)
cat("读取到", length(regulon.list), "个TF regulon\n")

# 显示每个regulon的基因数量
for (tf_name in names(regulon.list)) {
  cat("  ", tf_name, ": ", length(regulon.list[[tf_name]]), "个基因\n")
}

# 整理所有TF及其下游targets为一个基因集
# 包含TF本身和所有target基因
all_regulon_genes_raw <- unique(c(names(regulon.list), unlist(regulon.list)))

# 清理基因名称，去除(+)符号
all_regulon_genes <- gsub("\\(\\+\\)", "", all_regulon_genes_raw)
all_regulon_genes <- unique(all_regulon_genes)  # 去重，因为可能有重复

cat("TF regulon总基因数量（清理前）:", length(all_regulon_genes_raw), "\n")
cat("TF regulon总基因数量（清理后）:", length(all_regulon_genes), "\n")

# 创建基因集信息表
regulon_gene_info <- data.frame(
  Gene = all_regulon_genes,
  Type = ifelse(all_regulon_genes %in% names(regulon.list), "TF", "Target"),
  stringsAsFactors = FALSE
)

# 为每个基因添加所属的TF信息
regulon_gene_info$Associated_TFs <- sapply(regulon_gene_info$Gene, function(gene) {
  # 如果是TF本身
  if (gene %in% names(regulon.list)) {
    return(gene)
  }
  # 如果是target基因，找到包含它的TF
  tfs <- names(regulon.list)[sapply(regulon.list, function(targets) gene %in% targets)]
  return(paste(tfs, collapse = ";"))
})

cat("TF数量:", sum(regulon_gene_info$Type == "TF"), "\n")
cat("Target基因数量:", sum(regulon_gene_info$Type == "Target"), "\n")

# 保存regulon基因集信息
write.csv(regulon_gene_info, file.path(output_dir, "regulon_gene_info.csv"), row.names = FALSE)

# ================================================================================
# 步骤3: 基于regulon基因集生成TCGA表达矩阵并过滤低表达基因
# ================================================================================
cat("\n步骤3: 基于regulon基因集生成TCGA表达矩阵并过滤低表达基因...\n")

# 读取TCGA表达矩阵
cat("读取TCGA表达矩阵...\n")
tcga_expr <- read.csv(tcga_expr_file, row.names = 1, check.names = FALSE)
cat("TCGA表达矩阵维度:", dim(tcga_expr), "\n")

# 统一格式：把"-"换成"."并去掉首尾空白
fmt <- function(x) trimws(gsub("-", "\\.", x))
colnames(tcga_expr) <- fmt(colnames(tcga_expr))
rownames(tcga_clusters) <- fmt(rownames(tcga_clusters))

cat("格式化后TCGA表达矩阵维度:", dim(tcga_expr), "\n")
cat("格式化后TCGA分组信息维度:", dim(tcga_clusters), "\n")

# 生成基于regulon基因集的表达矩阵
regulon_genes_in_expr <- intersect(all_regulon_genes, rownames(tcga_expr))
tcga_expr_regulon <- tcga_expr[regulon_genes_in_expr, ]
cat("regulon基因在表达矩阵中的数量:", length(regulon_genes_in_expr), "/", length(all_regulon_genes), "\n")

# 过滤在少于30个样本中有表达的基因
cat("过滤在少于30个样本中有表达的基因...\n")
min_samples <- 30
expressed_samples_count <- rowSums(tcga_expr_regulon > 0)
genes_to_keep <- names(expressed_samples_count)[expressed_samples_count >= min_samples]

tcga_expr_regulon_filtered <- tcga_expr_regulon[genes_to_keep, ]
cat("过滤前基因数量:", nrow(tcga_expr_regulon), "\n")
cat("过滤后基因数量:", nrow(tcga_expr_regulon_filtered), "\n")
cat("被过滤的基因数量:", nrow(tcga_expr_regulon) - nrow(tcga_expr_regulon_filtered), "\n")

# 保存过滤后的表达矩阵
write.csv(tcga_expr_regulon_filtered, file.path(output_dir, "TCGA_expr_regulon_filtered.csv"))

# 更新基因列表为过滤后的基因列表
regulon_genes_filtered <- genes_to_keep

# 更新基因信息表，只保留过滤后的基因
regulon_gene_info_filtered <- regulon_gene_info[regulon_gene_info$Gene %in% regulon_genes_filtered, ]
write.csv(regulon_gene_info_filtered, file.path(output_dir, "regulon_gene_info_filtered.csv"), row.names = FALSE)

cat("过滤后TF数量:", sum(regulon_gene_info_filtered$Type == "TF"), "\n")
cat("过滤后Target基因数量:", sum(regulon_gene_info_filtered$Type == "Target"), "\n")

# ================================================================================
# 步骤4: 构建Lasso + Logistic回归预测模型
# ================================================================================
cat("\n步骤4: 构建Lasso + Logistic回归预测模型...\n")

# 准备建模数据
# 匹配样本
common_samples <- intersect(colnames(tcga_expr_regulon_filtered), rownames(tcga_clusters))
cat("共同样本数量:", length(common_samples), "\n")

# 准备特征矩阵和标签
X <- t(tcga_expr_regulon_filtered[, common_samples])
y <- ifelse(tcga_clusters[common_samples, "cluster"] == "cluster_high", 1, 0)

cat("特征矩阵维度:", dim(X), "\n")
cat("标签分布:\n")
print(table(y))

# 设置随机种子
set.seed(123)

# Lasso回归进行特征选择（加强筛选条件）
cat("进行Lasso回归特征选择（加强筛选条件）...\n")

# 使用更严格的交叉验证
cv_lasso <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10, 
                      type.measure = "auc")  # 使用AUC作为评价指标

# 使用lambda.1se而不是lambda.min来获得更稀疏的模型
optimal_lambda <- cv_lasso$lambda.1se  # 使用1se规则，获得更稀疏的模型
cat("最优lambda (1se规则):", optimal_lambda, "\n")
cat("lambda.min:", cv_lasso$lambda.min, "\n")

# 获取选中的特征
lasso_coef <- coef(cv_lasso, s = optimal_lambda)
selected_features_initial <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]  # 移除截距
cat("Lasso初步选择的特征数量:", length(selected_features_initial), "\n")

# 进一步筛选：只保留系数绝对值较大的特征
# 正确提取系数值
coef_matrix <- as.matrix(lasso_coef)
coef_values <- coef_matrix[selected_features_initial, 1]  # 直接使用特征名索引
names(coef_values) <- selected_features_initial

# 计算系数绝对值的阈值（保留前50%或至少5个特征）
coef_abs <- abs(coef_values)
threshold_percentile <- 0.5  # 保留系数绝对值前50%的特征
min_features <- 5  # 至少保留5个特征
max_features <- 15  # 最多保留15个特征

if (length(coef_abs) > max_features) {
  # 如果特征太多，取系数绝对值最大的前15个
  selected_features <- names(sort(coef_abs, decreasing = TRUE)[1:max_features])
} else if (length(coef_abs) > min_features) {
  # 如果特征数量适中，使用百分位数阈值
  threshold <- quantile(coef_abs, 1 - threshold_percentile)
  selected_features <- names(coef_abs[coef_abs >= threshold])
  
  # 确保至少有min_features个特征
  if (length(selected_features) < min_features) {
    selected_features <- names(sort(coef_abs, decreasing = TRUE)[1:min_features])
  }
} else {
  # 如果特征数量已经很少，保留所有特征
  selected_features <- selected_features_initial
}

cat("最终选择的特征数量:", length(selected_features), "\n")
cat("选择的特征:", paste(selected_features, collapse = ", "), "\n")
cat("特征系数绝对值范围:", round(min(abs(coef_values[selected_features])), 4), 
    "到", round(max(abs(coef_values[selected_features])), 4), "\n")

# 显示选中特征的TF/Target信息
selected_gene_info <- regulon_gene_info_filtered[regulon_gene_info_filtered$Gene %in% selected_features, ]
cat("\n选中特征的详细信息:\n")
print(selected_gene_info)

# 使用选中的特征构建逻辑回归模型
X_selected <- X[, selected_features, drop = FALSE]
logistic_model <- glm(y ~ ., data = data.frame(y = y, X_selected), family = "binomial")

# 输出模型系数
cat("\n逻辑回归模型系数:\n")
model_coef <- summary(logistic_model)$coefficients
print(model_coef)

# 构建风险评分公式
intercept <- model_coef["(Intercept)", "Estimate"]
gene_coefs <- model_coef[-1, "Estimate"]
names(gene_coefs) <- selected_features

# 生成公式字符串
formula_parts <- paste(round(gene_coefs, 4), "*", names(gene_coefs), collapse = " + ")
risk_formula <- paste("Risk Score =", round(intercept, 4), "+", formula_parts)
cat("\n风险评分公式:\n")
cat(risk_formula, "\n")

# 保存公式到文件
writeLines(risk_formula, file.path(output_dir, "risk_prediction_formula.txt"))

# 保存选中特征的详细信息
write.csv(selected_gene_info, file.path(output_dir, "selected_features_info.csv"), row.names = FALSE)

# ================================================================================
# 步骤4.1: 模型评价
# ================================================================================
cat("\n步骤4.1: 模型评价...\n")

# 计算预测概率
pred_prob <- predict(logistic_model, type = "response")
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

# 混淆矩阵
confusion_matrix <- table(Predicted = pred_class, Actual = y)
cat("混淆矩阵:\n")
print(confusion_matrix)

# 计算性能指标
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
sensitivity <- confusion_matrix[2,2] / sum(confusion_matrix[,2])
specificity <- confusion_matrix[1,1] / sum(confusion_matrix[,1])
precision <- confusion_matrix[2,2] / sum(confusion_matrix[2,])

cat("准确率 (Accuracy):", round(accuracy, 3), "\n")
cat("敏感性 (Sensitivity/Recall):", round(sensitivity, 3), "\n")
cat("特异性 (Specificity):", round(specificity, 3), "\n")
cat("精确度 (Precision):", round(precision, 3), "\n")

# ROC曲线分析
roc_obj <- roc(y, pred_prob)
auc_value <- auc(roc_obj)
cat("AUC值:", round(auc_value, 3), "\n")

# 绘制ROC曲线
pdf(file.path(output_dir, "ROC_curve.pdf"), width = 4, height = 3)
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"),
     col = "#E31A1C", lwd = 3, cex.main = 1.2, cex.lab = 1)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

# 绘制Lasso路径图
pdf(file.path(output_dir, "Lasso_path.pdf"), width = 5, height = 3)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(cv_lasso$glmnet.fit, "lambda", label = TRUE, main = "Lasso Regression Path")
plot(cv_lasso, main = "Cross-Validation for Lambda Selection")
dev.off()

# ================================================================================
# 步骤5: 对TCGA样本进行风险评分并绘制热图
# ================================================================================
cat("\n步骤5: 对TCGA样本进行风险评分并绘制热图...\n")

# 计算所有TCGA样本的风险评分
tcga_risk_scores <- intercept + as.matrix(t(tcga_expr_regulon_filtered[selected_features, common_samples])) %*% gene_coefs
tcga_risk_scores <- as.vector(tcga_risk_scores)
names(tcga_risk_scores) <- common_samples

# 标准化TCGA风险评分
tcga_risk_scores <- normalize_risk_scores(tcga_risk_scores)

# 按风险评分排序
tcga_risk_scores_sorted <- sort(tcga_risk_scores)
sorted_samples <- names(tcga_risk_scores_sorted)

# 准备热图数据
heatmap_data <- tcga_expr_regulon_filtered[selected_features, sorted_samples]

# 检查并清理热图数据
cat("TCGA热图数据检查:\n")
cat("  - 数据维度:", dim(heatmap_data), "\n")
cat("  - 缺失值数量:", sum(is.na(heatmap_data)), "\n")
cat("  - 无限值数量:", sum(is.infinite(as.matrix(heatmap_data))), "\n")

# 处理缺失值和无限值
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data[is.infinite(as.matrix(heatmap_data))] <- 0
heatmap_data <- as.matrix(heatmap_data)

if (!all(heatmap_data == 0)) {
  cat("  - 数据范围:", round(range(heatmap_data), 4), "\n")
}

# 准备注释信息
annotation_col <- data.frame(
  Cluster = tcga_clusters[sorted_samples, "cluster"],
  Risk_Score = tcga_risk_scores_sorted,
  row.names = sorted_samples
)

# 设置注释颜色
ann_colors <- list(
  Cluster = cluster_colors,
  Risk_Score = colorRamp2(c(min(tcga_risk_scores_sorted), 0, max(tcga_risk_scores_sorted)), 
                          c("blue", "white", "red"))
)

# 绘制热图
pdf(file.path(output_dir, "TCGA_risk_heatmap.pdf"), width = 8, height = 5)
pheatmap(heatmap_data,
         annotation_col = annotation_col,
         annotation_colors = list(Cluster = cluster_colors),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "TCGA Samples Risk Heatmap (Sorted by Risk Score)",
         fontsize = 10,           # 主要字体大小
         fontsize_row = 10,       # 行名字体大小
         fontsize_col = 10,       # 列名字体大小
         annotation_names_row = TRUE,
         annotation_names_col = TRUE)
dev.off()

# 保存风险评分结果
tcga_risk_results <- data.frame(
  Sample = names(tcga_risk_scores),
  Original_Risk_Score = tcga_risk_scores_original <- intercept + as.matrix(t(tcga_expr_regulon_filtered[selected_features, common_samples])) %*% gene_coefs,
  Normalized_Risk_Score = tcga_risk_scores,
  Cluster = tcga_clusters[names(tcga_risk_scores), "cluster"],
  stringsAsFactors = FALSE
)
write.csv(tcga_risk_results, file.path(output_dir, "TCGA_risk_scores.csv"), row.names = FALSE)

# ================================================================================
# 步骤6: 对OAK_POPLAR样本进行风险评分并绘制热图
# ================================================================================
cat("\n步骤6: 对OAK_POPLAR样本进行风险评分并绘制热图...\n")

# 读取OAK_POPLAR表达数据和分组信息
oak_poplar_expr <- read.csv(oak_poplar_expr_file, row.names = 1, check.names = FALSE)
oak_poplar_clusters <- read.csv(oak_poplar_cluster_file, row.names = 1, stringsAsFactors = FALSE)

cat("OAK_POPLAR表达矩阵维度:", dim(oak_poplar_expr), "\n")
cat("OAK_POPLAR分组统计:\n")
print(table(oak_poplar_clusters$cluster))

# 匹配样本和基因
common_samples_oak <- intersect(colnames(oak_poplar_expr), rownames(oak_poplar_clusters))
common_genes_oak <- intersect(selected_features, rownames(oak_poplar_expr))

cat("OAK_POPLAR共同样本数量:", length(common_samples_oak), "\n")
cat("OAK_POPLAR共同基因数量:", length(common_genes_oak), "/", length(selected_features), "\n")

# 计算风险评分（只使用可用的基因）
available_coefs <- gene_coefs[common_genes_oak]
oak_risk_scores <- intercept + as.matrix(t(oak_poplar_expr[common_genes_oak, common_samples_oak])) %*% available_coefs
oak_risk_scores <- as.vector(oak_risk_scores)
names(oak_risk_scores) <- common_samples_oak

# 标准化OAK_POPLAR风险评分
oak_risk_scores <- normalize_risk_scores(oak_risk_scores)

# 按风险评分排序
oak_risk_scores_sorted <- sort(oak_risk_scores)
oak_sorted_samples <- names(oak_risk_scores_sorted)

# 准备热图数据
oak_heatmap_data <- oak_poplar_expr[common_genes_oak, oak_sorted_samples]

# 检查并清理热图数据
cat("OAK_POPLAR热图数据检查:\n")
cat("  - 数据维度:", dim(oak_heatmap_data), "\n")
cat("  - 缺失值数量:", sum(is.na(oak_heatmap_data)), "\n")
cat("  - 无限值数量:", sum(is.infinite(as.matrix(oak_heatmap_data))), "\n")

# 处理缺失值和无限值
oak_heatmap_data[is.na(oak_heatmap_data)] <- 0
oak_heatmap_data[is.infinite(as.matrix(oak_heatmap_data))] <- 0
oak_heatmap_data <- as.matrix(oak_heatmap_data)

if (!all(oak_heatmap_data == 0)) {
  cat("  - 数据范围:", round(range(oak_heatmap_data), 4), "\n")
}

# 准备注释信息
oak_annotation_col <- data.frame(
  Cluster = oak_poplar_clusters[oak_sorted_samples, "cluster"],
  Risk_Score = oak_risk_scores_sorted,
  row.names = oak_sorted_samples
)

# 绘制热图
pdf(file.path(output_dir, "OAK_POPLAR_risk_heatmap.pdf"), width = 8, height = 5)
pheatmap(oak_heatmap_data,
         annotation_col = oak_annotation_col,
         annotation_colors = list(Cluster = cluster_colors),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "OAK_POPLAR Samples Risk Heatmap (Sorted by Risk Score)",
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE)
dev.off()

# 保存风险评分结果
oak_risk_results <- data.frame(
  Sample = names(oak_risk_scores),
  Original_Risk_Score = oak_risk_scores_original <- intercept + as.matrix(t(oak_poplar_expr[common_genes_oak, common_samples_oak])) %*% available_coefs,
  Normalized_Risk_Score = oak_risk_scores,
  Cluster = oak_poplar_clusters[names(oak_risk_scores), "cluster"],
  stringsAsFactors = FALSE
)
write.csv(oak_risk_results, file.path(output_dir, "OAK_POPLAR_risk_scores.csv"), row.names = FALSE)

# ================================================================================
# 步骤7: 对NT样本进行风险评分并绘制热图
# ================================================================================
cat("\n步骤7: 对NT样本进行风险评分并绘制热图...\n")

# 读取NT表达数据和分组信息
nt_expr <- read.csv(nt_expr_file, row.names = 1, check.names = FALSE)
nt_metadata <- read.csv(nt_metadata_file, stringsAsFactors = FALSE)

cat("NT表达矩阵维度:", dim(nt_expr), "\n")
cat("NT分组统计:\n")
print(table(nt_metadata$group))

# 匹配样本和基因
common_samples_nt <- intersect(colnames(nt_expr), nt_metadata$X)
common_genes_nt <- intersect(selected_features, rownames(nt_expr))

cat("NT共同样本数量:", length(common_samples_nt), "\n")
cat("NT共同基因数量:", length(common_genes_nt), "/", length(selected_features), "\n")

# 计算风险评分
available_coefs_nt <- gene_coefs[common_genes_nt]
nt_risk_scores <- intercept + as.matrix(t(nt_expr[common_genes_nt, common_samples_nt])) %*% available_coefs_nt
nt_risk_scores <- as.vector(nt_risk_scores)
names(nt_risk_scores) <- common_samples_nt

# 标准化NT风险评分
nt_risk_scores <- normalize_risk_scores(nt_risk_scores)

# 按风险评分排序
nt_risk_scores_sorted <- sort(nt_risk_scores)
nt_sorted_samples <- names(nt_risk_scores_sorted)

# 准备热图数据
nt_heatmap_data <- nt_expr[common_genes_nt, nt_sorted_samples]

# 检查并清理热图数据
cat("NT热图数据检查:\n")
cat("  - 数据维度:", dim(nt_heatmap_data), "\n")
cat("  - 缺失值数量:", sum(is.na(nt_heatmap_data)), "\n")
cat("  - 无限值数量:", sum(is.infinite(as.matrix(nt_heatmap_data))), "\n")

# 处理缺失值和无限值
nt_heatmap_data[is.na(nt_heatmap_data)] <- 0
nt_heatmap_data[is.infinite(as.matrix(nt_heatmap_data))] <- 0

# 确保数据为数值矩阵
nt_heatmap_data <- as.matrix(nt_heatmap_data)

# 检查是否还有问题
if (all(nt_heatmap_data == 0)) {
  cat("警告: 热图数据全为0，可能存在数据问题\n")
} else {
  cat("  - 数据范围:", round(range(nt_heatmap_data), 4), "\n")
}

# 准备注释信息
nt_annotation_col <- data.frame(
  Group = nt_metadata$group[match(nt_sorted_samples, nt_metadata$X)],
  Risk_Score = nt_risk_scores_sorted,
  row.names = nt_sorted_samples,
  stringsAsFactors = FALSE
)

# 检查风险评分的分布
cat("\nNT风险评分检查:\n")
cat("唯一值数量:", length(unique(nt_risk_scores_sorted)), "\n")
cat("是否有NA:", any(is.na(nt_risk_scores_sorted)), "\n")
cat("是否有Inf:", any(is.infinite(nt_risk_scores_sorted)), "\n")
cat("分位数:\n")
print(quantile(nt_risk_scores_sorted, probs = seq(0, 1, 0.25)))

# 设置风险评分的断点
# 使用四分位数而不是中位数，以确保有足够的不同值
risk_score_breaks <- unique(c(
  min(nt_risk_scores_sorted),
  quantile(nt_risk_scores_sorted, probs = c(0.25, 0.5, 0.75)),
  max(nt_risk_scores_sorted)
))

# 检查断点值
cat("\n风险评分断点值:\n")
print(risk_score_breaks)
cat("断点数量:", length(risk_score_breaks), "\n")
cat("断点是否唯一:", length(unique(risk_score_breaks)) == length(risk_score_breaks), "\n")

# 如果断点值太少，使用等间距的断点
if(length(unique(risk_score_breaks)) < 3) {
  cat("警告：断点值太少，使用等间距断点\n")
  score_range <- range(nt_risk_scores_sorted)
  risk_score_breaks <- seq(score_range[1], score_range[2], length.out = 5)
}

# 生成对应的颜色
risk_score_colors <- colorRampPalette(c("blue", "white", "red"))(length(risk_score_breaks))

# 检查并清理热图数据
cat("\nNT热图数据检查:\n")
cat("  - 数据维度:", dim(nt_heatmap_data), "\n")
cat("  - 缺失值数量:", sum(is.na(nt_heatmap_data)), "\n")
cat("  - 无限值数量:", sum(is.infinite(as.matrix(nt_heatmap_data))), "\n")

# 处理缺失值和无限值
nt_heatmap_data[is.na(nt_heatmap_data)] <- 0
nt_heatmap_data[is.infinite(as.matrix(nt_heatmap_data))] <- 0
nt_heatmap_data <- as.matrix(nt_heatmap_data)

# 设置NT分组颜色
nt_group_colors <- c(
  "C_LUAD" = "#ADB6CA",
  "T_LUAD" = "#FDAF91",
  "C_SCLC" = "#00468B",
  "T_SCLC" = "#ED0000"
)

# 打印分组统计信息
cat("\nNT分组统计（热图样本）:\n")
print(table(nt_annotation_col$Group, useNA = "ifany"))

# 绘制热图
pdf(file.path(output_dir, "NT_risk_heatmap.pdf"), width = 8, height = 5)
tryCatch({
  # 标准化前检查
  cat("\n标准化前的风险评分检查:\n")
  cat("总样本数:", length(nt_risk_scores), "\n")
  cat("NA值数量:", sum(is.na(nt_risk_scores)), "\n")
  cat("无限值数量:", sum(is.infinite(nt_risk_scores)), "\n")
  if(sum(is.finite(nt_risk_scores)) > 0) {
    cat("有效值范围:", range(nt_risk_scores[is.finite(nt_risk_scores)]), "\n")
  } else {
    cat("警告: 没有有效的风险评分值\n")
  }
  
  # 保存原始风险评分
  write.csv(data.frame(
    Sample = names(nt_risk_scores),
    Original_Risk_Score = nt_risk_scores,
    Group = nt_metadata$group[match(names(nt_risk_scores), nt_metadata$X)]
  ), file.path(output_dir, "NT_risk_scores_original.csv"), row.names = FALSE)
  
  # 标准化风险评分
  nt_risk_scores_normalized <- normalize_risk_scores(nt_risk_scores)
  
  # 标准化后检查
  cat("\n标准化后的风险评分检查:\n")
  cat("总样本数:", length(nt_risk_scores_normalized), "\n")
  cat("NA值数量:", sum(is.na(nt_risk_scores_normalized)), "\n")
  cat("无限值数量:", sum(is.infinite(nt_risk_scores_normalized)), "\n")
  cat("值的范围:", range(nt_risk_scores_normalized), "\n")
  
  # 保存标准化后的风险评分
  write.csv(data.frame(
    Sample = names(nt_risk_scores_normalized),
    Normalized_Risk_Score = nt_risk_scores_normalized,
    Group = nt_metadata$group[match(names(nt_risk_scores_normalized), nt_metadata$X)]
  ), file.path(output_dir, "NT_risk_scores_normalized.csv"), row.names = FALSE)
  
  # 按风险评分排序
  nt_risk_scores_sorted <- sort(nt_risk_scores_normalized)
  nt_sorted_samples <- names(nt_risk_scores_sorted)
  
  # 准备热图数据
  nt_heatmap_data <- nt_expr[common_genes_nt, nt_sorted_samples]
  
  # 准备注释信息
  nt_annotation_col <- data.frame(
    Group = nt_metadata$group[match(nt_sorted_samples, nt_metadata$X)],
    Risk_Score = nt_risk_scores_sorted,
    row.names = nt_sorted_samples,
    stringsAsFactors = FALSE
  )
  
  # 设置风险评分的断点（使用5个等间距的点）
  score_min <- min(nt_risk_scores_sorted)
  score_max <- max(nt_risk_scores_sorted)
  risk_score_breaks <- seq(score_min, score_max, length.out = 5)
  
  # 生成对应的颜色
  risk_score_colors <- colorRampPalette(c("blue", "white", "red"))(5)
  
  # 检查并清理热图数据
  cat("\nNT热图数据检查:\n")
  cat("  - 数据维度:", dim(nt_heatmap_data), "\n")
  cat("  - 缺失值数量:", sum(is.na(nt_heatmap_data)), "\n")
  cat("  - 无限值数量:", sum(is.infinite(as.matrix(nt_heatmap_data))), "\n")
  cat("  - 数据范围:", range(nt_heatmap_data), "\n")
  
  # 处理缺失值和无限值
  nt_heatmap_data[is.na(nt_heatmap_data)] <- 0
  nt_heatmap_data[is.infinite(as.matrix(nt_heatmap_data))] <- 0
  nt_heatmap_data <- as.matrix(nt_heatmap_data)
  
  # 设置NT分组颜色
  nt_group_colors <- c(
    "C_LUAD" = "#ADB6CA",
    "T_LUAD" = "#FDAF91",
    "C_SCLC" = "#00468B",
    "T_SCLC" = "#ED0000"
  )
  
  # 打印分组统计信息
  cat("\nNT分组统计（热图样本）:\n")
  print(table(nt_annotation_col$Group, useNA = "ifany"))
  
  # 绘制热图
  print(pheatmap(nt_heatmap_data,
                 annotation_col = nt_annotation_col,
                 annotation_colors = list(
                   Group = nt_group_colors,
                   Risk_Score = circlize::colorRamp2(
                     breaks = risk_score_breaks,
                     colors = risk_score_colors
                   )
                 ),
                 cluster_cols = FALSE,
                 cluster_rows = TRUE,
                 show_colnames = FALSE,
                 show_rownames = TRUE,
                 scale = "row",
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 main = "NT Samples Risk Heatmap (Sorted by Risk Score)",
                 fontsize = 10,
                 fontsize_row = 10,
                 fontsize_col = 10,
                 annotation_names_row = TRUE,
                 annotation_names_col = TRUE))
  dev.off()
  
  # 保存完整的风险评分结果
  write.csv(data.frame(
    Sample = names(nt_risk_scores),
    Original_Risk_Score = nt_risk_scores,
    Normalized_Risk_Score = nt_risk_scores_normalized,
    Group = nt_metadata$group[match(names(nt_risk_scores), nt_metadata$X)]
  ), file.path(output_dir, "NT_risk_scores_complete.csv"), row.names = FALSE)
  
}, error = function(e) {
  cat("错误:", conditionMessage(e), "\n")
})

# 风险评分数据已在tryCatch块中保存

# 打印维度检查信息
cat("\nNT数据维度检查:\n")
cat("热图数据维度:", dim(nt_heatmap_data), "\n")
cat("注释数据维度:", dim(nt_annotation_col), "\n")
cat("样本数量匹配:", ncol(nt_heatmap_data) == nrow(nt_annotation_col), "\n")
cat("样本顺序匹配:", all(colnames(nt_heatmap_data) == rownames(nt_annotation_col)), "\n")

# 保存风险评分结果
nt_risk_results <- data.frame(
  Sample = names(nt_risk_scores),
  Original_Risk_Score = nt_risk_scores_original <- intercept + as.matrix(t(nt_expr[common_genes_nt, common_samples_nt])) %*% available_coefs_nt,
  Normalized_Risk_Score = nt_risk_scores,
  Group = nt_metadata$group[match(names(nt_risk_scores), nt_metadata$X)],
  stringsAsFactors = FALSE
)
write.csv(nt_risk_results, file.path(output_dir, "NT_risk_scores.csv"), row.names = FALSE)

# ================================================================================
# 步骤7.1: 对NT样本中的C_LUAD和T_LUAD分组单独绘制热图
# ================================================================================
cat("\n步骤7.1: 对NT样本中的C_LUAD和T_LUAD分组单独绘制热图...\n")

# 筛选C_LUAD和T_LUAD样本
luad_groups <- c("C_LUAD", "T_LUAD")

# 获取所有NT样本的分组信息
all_nt_sample_groups <- nt_metadata$group[match(nt_sorted_samples, nt_metadata$X)]

# 筛选LUAD样本
luad_mask <- all_nt_sample_groups %in% luad_groups
luad_samples <- nt_sorted_samples[luad_mask]
luad_sample_groups <- all_nt_sample_groups[luad_mask]

cat("所有NT样本分组统计:\n")
print(table(all_nt_sample_groups, useNA = "ifany"))
cat("LUAD样本数量:", length(luad_samples), "\n")
cat("LUAD分组统计:\n")
print(table(luad_sample_groups, useNA = "ifany"))

# 检查是否有LUAD样本
if (length(luad_samples) > 0) {
  # 准备LUAD热图数据
  cat("准备LUAD热图数据...\n")
  cat("  - 可用基因数量:", length(common_genes_nt), "\n")
  cat("  - LUAD样本列表:", paste(luad_samples, collapse = ", "), "\n")
  
  # 检查样本是否在表达矩阵中
  samples_in_expr <- luad_samples %in% colnames(nt_expr)
  cat("  - 样本在表达矩阵中:", sum(samples_in_expr), "/", length(luad_samples), "\n")
  
  if (sum(samples_in_expr) == 0) {
    cat("错误: 没有LUAD样本在表达矩阵中找到\n")
    cat("表达矩阵列名示例:", paste(head(colnames(nt_expr)), collapse = ", "), "\n")
  } else {
    # 获取有效样本
    luad_samples_valid <- luad_samples[samples_in_expr]
    luad_sample_groups_valid <- luad_sample_groups[samples_in_expr]
    
    # 提取热图数据
    luad_heatmap_data <- nt_expr[common_genes_nt, luad_samples_valid]
    luad_risk_scores_subset <- nt_risk_scores_sorted[luad_samples_valid]
    
    # 检查并清理热图数据
    cat("LUAD热图数据检查:\n")
    cat("  - 数据维度:", dim(luad_heatmap_data), "\n")
    cat("  - 缺失值数量:", sum(is.na(luad_heatmap_data)), "\n")
    cat("  - 无限值数量:", sum(is.infinite(as.matrix(luad_heatmap_data))), "\n")
    
    # 处理缺失值和无限值
    luad_heatmap_data[is.na(luad_heatmap_data)] <- 0
    luad_heatmap_data[is.infinite(as.matrix(luad_heatmap_data))] <- 0
    luad_heatmap_data <- as.matrix(luad_heatmap_data)
    
    # 检查数据范围
    if (!all(luad_heatmap_data == 0)) {
      cat("  - 数据范围:", round(range(luad_heatmap_data), 4), "\n")
    }
    
    # 准备注释信息
    luad_annotation_col <- data.frame(
      Group = as.character(luad_sample_groups_valid),
      row.names = luad_samples_valid,
      stringsAsFactors = FALSE
    )
    
    # 设置分组颜色
    luad_group_colors <- c(
      "C_LUAD" = "#ADB6CA",
      "T_LUAD" = "#FDAF91"
    )
    
    # 绘制热图
    pdf(file.path(output_dir, "NT_LUAD_risk_heatmap.pdf"), width = 12, height = 10)
    cat("绘制LUAD热图...\n")
    print(pheatmap(
      luad_heatmap_data,
      annotation_col = luad_annotation_col,
      annotation_colors = list(Group = luad_group_colors),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      scale = "row",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      main = "NT LUAD Samples Risk Heatmap (C_LUAD vs T_LUAD)"
    ))
    dev.off()
    
    # 保存风险评分结果
    luad_risk_results <- data.frame(
      Sample = luad_samples_valid,
      Original_Risk_Score = luad_risk_scores_original <- intercept + as.matrix(t(nt_expr[common_genes_nt, luad_samples_valid])) %*% available_coefs_nt,
      Normalized_Risk_Score = luad_risk_scores_subset,
      Group = luad_sample_groups_valid,
      stringsAsFactors = FALSE
    )
    write.csv(luad_risk_results, 
              file.path(output_dir, "NT_LUAD_risk_scores.csv"), 
              row.names = FALSE)
    
    cat("LUAD热图和风险评分结果已保存\n")
  }
} else {
  cat("警告: 没有找到C_LUAD或T_LUAD样本\n")
}

# ================================================================================
# 步骤8: 绘制CNS级别的列线图
# ================================================================================
cat("\n步骤8: 绘制CNS级别的列线图...\n")

# 准备列线图数据
nomogram_data <- data.frame(
  y = y,
  X_selected
)

# 使用rms包重新拟合模型
dd <- datadist(nomogram_data)
options(datadist = "dd")

# 拟合logistic回归模型
lrm_model <- lrm(y ~ ., data = nomogram_data, x = TRUE, y = TRUE)

# 创建列线图
pdf(file.path(output_dir, "Nomogram_CNS_level.pdf"), width = 16, height = 12)

# 设置绘图参数，增加底部边距
par(mar = c(6, 4, 4, 2))

# 创建列线图对象
tryCatch({
  nom <- nomogram(lrm_model, 
                  fun = plogis,
                  fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),  # 减少刻度点，避免拥挤
                  funlabel = "Risk of cluster_high",
                  lp = FALSE)
  
  # 绘制列线图，调整参数以改善显示
  plot(nom, 
       cex.axis = 0.9,     # 增大轴标签字体
       cex.var = 1.1,      # 增大变量名字体
       lmgp = 0.4,         # 调整标签边距
       xfrac = 0.35)       # 减小刻度线长度，给标签更多空间
  
  # 添加标题
  title(main = "Nomogram for Predicting cluster_high Risk", 
        cex.main = 1.4, 
        line = 1.5)
  
  cat("列线图绘制成功\n")
  
}, error = function(e) {
  cat("列线图绘制失败:", e$message, "\n")
  
  # 如果nomogram失败，绘制一个简单的系数图作为替代
  plot(1:length(gene_coefs), gene_coefs,
       type = "h",
       xlab = "Gene Index",
       ylab = "Coefficient",
       main = "Gene Coefficients in Risk Model",
       xaxt = "n")
  axis(1, at = 1:length(gene_coefs), labels = names(gene_coefs), las = 2)
  abline(h = 0, lty = 2, col = "gray")
})

dev.off()

# ================================================================================
# 步骤8.1: 添加基因评分的其他可视化
# ================================================================================
cat("\n步骤8.1: 添加基因评分的其他可视化...\n")

# 提取每个基因的系数
gene_coef_data <- data.frame(
  Gene = names(gene_coefs),
  Coefficient = as.numeric(gene_coefs),
  AbsCoef = abs(gene_coefs),
  Direction = ifelse(gene_coefs > 0, "Positive", "Negative")
)

# 按系数绝对值排序
gene_coef_data <- gene_coef_data[order(-gene_coef_data$AbsCoef), ]

# 1. 瀑布图：显示基因的贡献度及方向
pdf(file.path(output_dir, "Gene_contribution_waterfall.pdf"), width = 10, height = 8)
par(mar = c(8, 4, 4, 2))
barplot(gene_coef_data$Coefficient,
        names.arg = gene_coef_data$Gene,
        col = ifelse(gene_coef_data$Coefficient > 0, "#E31A1C", "#3F51B5"),
        las = 2,
        main = "Gene Contribution to Risk Score",
        ylab = "Coefficient",
        border = NA)
abline(h = 0, lty = 2, col = "gray50")
legend("topright",
       legend = c("Positive", "Negative"),
       fill = c("#E31A1C", "#3F51B5"),
       border = NA)
dev.off()

# 2. 点图：显示基因系数的分布
pdf(file.path(output_dir, "Gene_coefficient_dotplot.pdf"), width = 10, height = 8)
par(mar = c(8, 4, 4, 2))
plot(1:nrow(gene_coef_data), gene_coef_data$Coefficient,
     pch = 19,
     col = ifelse(gene_coef_data$Coefficient > 0, "#E31A1C", "#3F51B5"),
     xlab = "",
     ylab = "Coefficient",
     main = "Gene Coefficient Distribution",
     xaxt = "n")
axis(1, at = 1:nrow(gene_coef_data), labels = gene_coef_data$Gene, las = 2)
abline(h = 0, lty = 2, col = "gray50")
points(1:nrow(gene_coef_data), gene_coef_data$Coefficient,
       pch = 21,
       bg = ifelse(gene_coef_data$Coefficient > 0, "#E31A1C", "#3F51B5"))
dev.off()

# 3. 热力图：显示基因间的相关性
gene_expr_selected <- t(tcga_expr_regulon_filtered[gene_coef_data$Gene, ])
gene_cor <- cor(gene_expr_selected, method = "spearman")

pdf(file.path(output_dir, "Gene_correlation_heatmap.pdf"), width = 10, height = 10)
pheatmap(gene_cor,
         color = colorRampPalette(c("#3F51B5", "white", "#E31A1C"))(100),
         main = "Gene Expression Correlation",
         fontsize = 10,
         fontsize_row = 8,
         fontsize_col = 8,
         clustering_method = "ward.D2")
dev.off()

# 4. 箱线图：比较不同cluster中基因的表达水平
pdf(file.path(output_dir, "Gene_expression_by_cluster.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))

# 获取共同样本
common_samples_for_plot <- intersect(colnames(tcga_expr_regulon_filtered), rownames(tcga_clusters))

for(gene in gene_coef_data$Gene[1:min(4, nrow(gene_coef_data))]) {
  # 提取基因表达数据
  gene_expr <- as.numeric(tcga_expr_regulon_filtered[gene, common_samples_for_plot])
  # 提取对应的cluster信息
  gene_clusters <- tcga_clusters[common_samples_for_plot, "cluster"]
  
  # 创建数据框
  plot_data <- data.frame(
    Expression = gene_expr,
    Cluster = gene_clusters
  )
  
  # 绘制箱线图
  boxplot(Expression ~ Cluster, 
          data = plot_data,
          main = paste("Expression of", gene),
          xlab = "Cluster",
          ylab = "Expression Level",
          col = cluster_colors[unique(gene_clusters)],
          las = 2,
          cex.main = 1.2)
  
  # 添加样本数量信息
  cluster_counts <- table(gene_clusters)
  mtext(paste("n =", paste(cluster_counts, collapse = ", ")), 
        side = 1, line = 4, cex = 0.8)
}
dev.off()

cat("额外的可视化图表已保存\n")

# 绘制校准曲线
pdf(file.path(output_dir, "Calibration_plot.pdf"), width = 10, height = 8)
cal <- calibrate(lrm_model, method = "boot", B = 200)
plot(cal, 
     main = "Calibration Plot for Risk Prediction Model",
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     cex.main = 1.5,
     cex.lab = 1.3,
     cex.axis = 1.2,
     col = "#E31A1C",
     lwd = 3)
abline(0, 1, lty = 2, col = "gray", lwd = 2)
dev.off()

# 绘制决策曲线分析（DCA）
# 创建决策曲线数据
threshold_range <- seq(0.01, 0.99, by = 0.01)  # 避免0和1的极值
net_benefit <- numeric(length(threshold_range))
net_benefit_all <- numeric(length(threshold_range))
net_benefit_none <- numeric(length(threshold_range))

for (i in seq_along(threshold_range)) {
  threshold <- threshold_range[i]
  
  # 避免除零错误
  if (threshold >= 1) threshold <- 0.99
  if (threshold <= 0) threshold <- 0.01
  
  # 模型的净效益
  tp <- sum(pred_prob >= threshold & y == 1)
  fp <- sum(pred_prob >= threshold & y == 0)
  fn <- sum(pred_prob < threshold & y == 1)
  tn <- sum(pred_prob < threshold & y == 0)
  
  # 计算净效益，避免无限值
  net_benefit[i] <- (tp - fp * threshold / (1 - threshold)) / length(y)
  
  # 全部治疗的净效益
  net_benefit_all[i] <- (sum(y) - (length(y) - sum(y)) * threshold / (1 - threshold)) / length(y)
  
  # 不治疗的净效益
  net_benefit_none[i] <- 0
}

# 处理无限值和NaN值
net_benefit[!is.finite(net_benefit)] <- 0
net_benefit_all[!is.finite(net_benefit_all)] <- 0
net_benefit_none[!is.finite(net_benefit_none)] <- 0

# 计算合理的y轴范围
all_benefits <- c(net_benefit, net_benefit_all, net_benefit_none)
y_min <- max(-0.1, min(all_benefits, na.rm = TRUE))  # 设置下限
y_max <- min(0.5, max(all_benefits, na.rm = TRUE))   # 设置上限

# 绘制决策曲线
pdf(file.path(output_dir, "Decision_curve_analysis.pdf"), width = 10, height = 8)
plot(threshold_range, net_benefit, 
     type = "l", lwd = 3, col = "#E31A1C",
     xlab = "Threshold Probability", 
     ylab = "Net Benefit",
     main = "Decision Curve Analysis",
     cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2,
     ylim = c(y_min, y_max),
     xlim = c(0, 1))

lines(threshold_range, net_benefit_all, lwd = 2, col = "blue", lty = 2)
lines(threshold_range, net_benefit_none, lwd = 2, col = "gray", lty = 3)

legend("topright", 
       legend = c("Prediction Model", "Treat All", "Treat None"),
       col = c("#E31A1C", "blue", "gray"),
       lty = c(1, 2, 3),
       lwd = c(3, 2, 2),
       cex = 1.2)
dev.off()

# 绘制特征重要性图
feature_importance <- abs(gene_coefs)
names(feature_importance) <- selected_features

pdf(file.path(output_dir, "Feature_importance.pdf"), width = 12, height = 8)
par(mar = c(8, 5, 4, 2))
barplot(sort(feature_importance, decreasing = TRUE),
        las = 2,
        col = "#3F51B5",
        main = "Feature Importance in Risk Prediction Model",
        ylab = "Absolute Coefficient Value",
        cex.main = 1.5,
        cex.lab = 1.3,
        cex.axis = 1.2,
        cex.names = 1.1)
dev.off()

# 绘制风险评分分布图
pdf(file.path(output_dir, "Risk_score_distribution.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2))

# TCGA
hist(tcga_risk_scores, 
     breaks = 30, 
     col = color_palette["TCGA"], 
     main = "TCGA Risk Score Distribution",
     xlab = "Risk Score",
     cex.main = 1.3)

# OAK_POPLAR
hist(oak_risk_scores, 
     breaks = 30, 
     col = color_palette["OAK"], 
     main = "OAK_POPLAR Risk Score Distribution",
     xlab = "Risk Score",
     cex.main = 1.3)

# NT
hist(nt_risk_scores, 
     breaks = 30, 
     col = color_palette["NT_bulk"], 
     main = "NT Risk Score Distribution",
     xlab = "Risk Score",
     cex.main = 1.3)

# 合并分布
all_scores <- c(tcga_risk_scores, oak_risk_scores, nt_risk_scores)
hist(all_scores, 
     breaks = 30, 
     col = "lightgray", 
     main = "Combined Risk Score Distribution",
     xlab = "Risk Score",
     cex.main = 1.3)
dev.off()

# ================================================================================
# 步骤8.2: 不同组间风险评分的统计分析和可视化
# ================================================================================
cat("\n步骤8.2: 进行组间风险评分的统计分析和可视化...\n")
stats_summary <- list()
# 初始化报告内容容器
report_content <- character()

# 加载所需的包
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggpubr)) install.packages("ggpubr")
library(ggplot2)
library(ggpubr)

# 创建可视化输出目录
viz_dir <- file.path(output_dir, "group_comparisons")
if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

# 函数：执行统计检验并返回结果
perform_statistical_tests <- function(data, group_var, value_var) {
  # 对所有组进行Kruskal-Wallis检验
  kw_test <- kruskal.test(formula(paste(value_var, "~", group_var)), data = data)
  
  # 两两比较的Wilcoxon检验
  groups <- unique(data[[group_var]])
  pairwise_tests <- list()
  fold_changes <- list()
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      
      # Wilcoxon检验
      test <- wilcox.test(data[[value_var]][data[[group_var]] == g1],
                          data[[value_var]][data[[group_var]] == g2])
      
      # 计算倍数变化（使用中位数）
      fc <- median(data[[value_var]][data[[group_var]] == g2]) / 
        median(data[[value_var]][data[[group_var]] == g1])
      
      pairwise_tests[[paste(g1, "vs", g2)]] <- test$p.value
      fold_changes[[paste(g1, "vs", g2)]] <- fc
    }
  }
  
  return(list(
    kruskal = kw_test,
    pairwise = pairwise_tests,
    fold_changes = fold_changes
  ))
}

# 添加检查数据分布的函数
check_normality <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  results <- list()
  
  for (group in groups) {
    group_data <- data[[value_var]][data[[group_var]] == group]
    # Shapiro-Wilk正态性检验
    sw_test <- shapiro.test(group_data)
    results[[group]] <- list(
      n = length(group_data),
      shapiro_p = sw_test$p.value
    )
  }
  return(results)
}

# 修改create_violin_boxplot函数，添加统计注释
create_violin_boxplot <- function(data, group_var, value_var, title, group_colors) {
  # 检查数据分布
  normality_results <- check_normality(data, group_var, value_var)
  cat("\n数据分布检验结果 (", title, "):\n")
  for (group in names(normality_results)) {
    cat(sprintf("%s组 (n=%d): Shapiro-Wilk p值 = %.4f %s\n",
                group,
                normality_results[[group]]$n,
                normality_results[[group]]$shapiro_p,
                ifelse(normality_results[[group]]$shapiro_p < 0.05,
                       "（不符合正态分布）",
                       "（符合正态分布）")))
  }
  
  # 根据数据类型设置正确的排序
  if (all(c("cluster_low", "cluster_medium", "cluster_high") %in% unique(data[[group_var]]))) {
    data[[group_var]] <- factor(data[[group_var]], 
                                levels = c("cluster_high", "cluster_medium", "cluster_low"))
  } else if (all(c("C_LUAD", "T_LUAD", "T_SCLC") %in% unique(data[[group_var]]))) {
    data[[group_var]] <- factor(data[[group_var]], 
                                levels = c("T_SCLC", "T_LUAD", "C_LUAD"))
  }
  
  # 获取排序后的组别
  groups <- levels(factor(data[[group_var]]))
  
  # 计算每个组的最大值和整体统计
  group_max <- tapply(data[[value_var]], data[[group_var]], max)
  overall_max <- max(group_max)
  y_range <- range(data[[value_var]])
  y_span <- diff(y_range)
  y_max <- overall_max + y_span * 0.8
  
  # 计算中位数和倍数变化
  medians <- tapply(data[[value_var]], data[[group_var]], median)
  means <- tapply(data[[value_var]], data[[group_var]], mean)
  sds <- tapply(data[[value_var]], data[[group_var]], sd)
  fold_changes <- list()
  stat_results <- list()
  
  # 创建详细的统计结果
  detailed_stats <- data.frame(
    Group = names(medians),
    N = tapply(data[[value_var]], data[[group_var]], length),
    Mean = means,
    SD = sds,
    Median = medians,
    Q1 = tapply(data[[value_var]], data[[group_var]], function(x) quantile(x, 0.25)),
    Q3 = tapply(data[[value_var]], data[[group_var]], function(x) quantile(x, 0.75)),
    Min = tapply(data[[value_var]], data[[group_var]], min),
    Max = tapply(data[[value_var]], data[[group_var]], max)
  )
  
  # 生成所有可能的两两比较组合
  comparisons <- list()
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      comparisons[[length(comparisons) + 1]] <- c(groups[i], groups[j])
    }
  }
  
  # 计算所有组间的统计检验结果和倍数变化
  comparison_results <- list()
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      
      # 计算倍数变化
      med1 <- medians[g1]
      med2 <- medians[g2]
      if (med1 >= med2) {
        fc <- med1 / med2
        direction <- "上调"
      } else {
        fc <- med2 / med1
        direction <- "下调"
      }
      
      # Wilcoxon检验
      wilcox_test <- wilcox.test(
        data[[value_var]][data[[group_var]] == g1],
        data[[value_var]][data[[group_var]] == g2]
      )
      
      # t检验
      t_test <- t.test(
        data[[value_var]][data[[group_var]] == g1],
        data[[value_var]][data[[group_var]] == g2]
      )
      
      comparison_name <- paste(g1, "vs", g2)
      fold_changes[[comparison_name]] <- sprintf("%.2f (%s)", fc, direction)
      
      comparison_results[[comparison_name]] <- list(
        group1 = g1,
        group2 = g2,
        median1 = med1,
        median2 = med2,
        fold_change = fc,
        direction = direction,
        wilcox_p = wilcox_test$p.value,
        t_test_p = t_test$p.value
      )
    }
  }
  
  # 创建比较结果数据框
  comparison_df <- do.call(rbind, lapply(names(comparison_results), function(comp) {
    res <- comparison_results[[comp]]
    data.frame(
      Comparison = comp,
      Group1_Median = res$median1,
      Group2_Median = res$median2,
      Fold_Change = res$fold_change,
      Direction = res$direction,
      Wilcoxon_P = res$wilcox_p,
      T_test_P = res$t_test_p,
      stringsAsFactors = FALSE
    )
  }))
  
  # 保存详细结果
  results_dir <- file.path(output_dir, "statistical_results")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # 生成文件名前缀（去除特殊字符）
  file_prefix <- gsub("[^[:alnum:]]", "_", title)
  
  # 保存基本统计量
  write.csv(detailed_stats, 
            file.path(results_dir, paste0(file_prefix, "_basic_stats.csv")), 
            row.names = TRUE)
  
  # 保存比较结果
  write.csv(comparison_df, 
            file.path(results_dir, paste0(file_prefix, "_comparison_results.csv")), 
            row.names = FALSE)
  
  # 打印统计结果
  cat("\n基本统计量:\n")
  print(detailed_stats)
  
  cat("\n组间比较结果:\n")
  print(comparison_df)
  
  # 创建基础图形
  p <- ggplot(data, aes_string(x = group_var, y = value_var, fill = group_var)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.margin = unit(c(2, 1, 0.5, 1), "cm")
    ) +
    ggtitle(title) +
    xlab("Group") +
    ylab("Risk Score") +
    coord_cartesian(ylim = c(y_range[1], y_max))
  
  # 添加所有两两比较的显著性标记
  p <- p + stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.signif",
    step.increase = 0.12,
    vjust = 2,
    bracket.size = 0.7,
    tip.length = 0.01,
    bracket.shorten = 0.05,
    hide.ns = FALSE,  # 显示非显著性结果
    symnum.args = list(  # 自定义显著性标记的阈值
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  )
  
  # 添加Kruskal-Wallis检验结果
  kw_test <- kruskal.test(formula(paste(value_var, "~", group_var)), data = data)
  kw_text <- sprintf("Kruskal-Wallis, p = %.2e", kw_test$p.value)
  p <- p + annotate("text",
                    x = length(groups)/2,
                    y = y_max * 0.95,
                    label = kw_text,
                    size = 4)
  
  # 格式化倍数变化文本
  fc_text <- paste(names(fold_changes), 
                   sprintf(": %s", unlist(fold_changes)), 
                   collapse = "\n")
  
  # 添加倍数变化注释
  p <- p + annotate("text",
                    x = length(groups),
                    y = y_max * 0.85,
                    label = fc_text,
                    size = 4,
                    hjust = 1)
  
  return(p)
}

# 1. TCGA数据分析
tcga_data <- data.frame(
  Group = tcga_clusters[names(tcga_risk_scores), "cluster"],
  Risk_Score = tcga_risk_scores,
  stringsAsFactors = FALSE
)

tcga_stats <- perform_statistical_tests(tcga_data, "Group", "Risk_Score")
cat("\nTCGA组间统计分析结果:\n")
cat("Kruskal-Wallis检验 p值:", tcga_stats$kruskal$p.value, "\n")
cat("组间倍数变化:\n")
print(unlist(tcga_stats$fold_changes))
cat("Wilcoxon检验 p值:\n")
print(unlist(tcga_stats$pairwise))

# 绘制TCGA小提琴箱线图
p_tcga <- create_violin_boxplot(tcga_data, "Group", "Risk_Score", 
                                "TCGA Risk Score Distribution by Cluster", 
                                cluster_colors)
ggsave(file.path(viz_dir, "TCGA_risk_violin_box.pdf"), p_tcga, 
       width = 5, height = 8)  # 增加高度

# 2. OAK_POPLAR数据分析
oak_data <- data.frame(
  Group = oak_poplar_clusters[names(oak_risk_scores), "cluster"],
  Risk_Score = oak_risk_scores,
  stringsAsFactors = FALSE
)

oak_stats <- perform_statistical_tests(oak_data, "Group", "Risk_Score")
cat("\nOAK_POPLAR组间统计分析结果:\n")
cat("Kruskal-Wallis检验 p值:", oak_stats$kruskal$p.value, "\n")
cat("组间倍数变化:\n")
print(unlist(oak_stats$fold_changes))
cat("Wilcoxon检验 p值:\n")
print(unlist(oak_stats$pairwise))

# 绘制OAK_POPLAR小提琴箱线图
p_oak <- create_violin_boxplot(oak_data, "Group", "Risk_Score",
                               "OAK_POPLAR Risk Score Distribution by Cluster",
                               cluster_colors)
ggsave(file.path(viz_dir, "OAK_POPLAR_risk_violin_box.pdf"), p_oak,
       width = 5, height = 8)

# 3. NT数据分析（包括C_LUAD、T_LUAD、T_SCLC三组）
nt_data <- data.frame(
  Group = nt_metadata$group[match(names(nt_risk_scores), nt_metadata$X)],
  Risk_Score = nt_risk_scores,
  stringsAsFactors = FALSE
)

# 只保留需要的组别
nt_data <- nt_data[nt_data$Group %in% c("C_LUAD", "T_LUAD", "T_SCLC"), ]

# 设置NT特定的颜色
nt_colors <- c(
  "C_LUAD" = "#ADB6CA",
  "T_LUAD" = "#FDAF91",
  "T_SCLC" = "#ED0000"
)

# 设置分组顺序
nt_data$Group <- factor(nt_data$Group, levels = c("C_LUAD", "T_LUAD", "T_SCLC"))

nt_stats <- perform_statistical_tests(nt_data, "Group", "Risk_Score")
cat("\nNT组间统计分析结果:\n")
cat("Kruskal-Wallis检验 p值:", nt_stats$kruskal$p.value, "\n")
cat("组间倍数变化:\n")
print(unlist(nt_stats$fold_changes))
cat("Wilcoxon检验 p值:\n")
print(unlist(nt_stats$pairwise))

# 绘制NT小提琴箱线图
p_nt <- create_violin_boxplot(nt_data, "Group", "Risk_Score",
                              "NT Risk Score Distribution by Group",
                              nt_colors)
ggsave(file.path(viz_dir, "NT_risk_violin_box.pdf"), p_nt,
       width = 5, height = 8)

# 更新统计分析结果
stats_summary$NT <- nt_stats

# 创建组合图（更新布局）
combined_plot <- ggpubr::ggarrange(p_tcga, p_oak, p_nt,
                                   ncol = 3, nrow = 1,
                                   labels = c("A", "B", "C"))

ggsave(file.path(viz_dir, "combined_risk_violin_box.pdf"), combined_plot,
       width = 15, height = 8)  # 调整组合图高度

# 将统计结果添加到报告中
add_stats_to_report <- function(stats, group_name) {
  result <- sprintf("\n%s组间比较:\n", group_name)
  result <- paste0(result, "- Kruskal-Wallis检验 p值: ", 
                   format(stats$kruskal$p.value, scientific = TRUE, digits = 3), "\n")
  result <- paste0(result, "- 组间倍数变化:\n")
  for (name in names(stats$fold_changes)) {
    result <- paste0(result, sprintf("  * %s: %.2f倍\n", 
                                     name, stats$fold_changes[[name]]))
  }
  result <- paste0(result, "- 组间Wilcoxon检验 p值:\n")
  for (name in names(stats$pairwise)) {
    result <- paste0(result, sprintf("  * %s: %s\n", 
                                     name, format(stats$pairwise[[name]], 
                                                  scientific = TRUE, digits = 3)))
  }
  return(result)
}

# 更新报告内容
report_content <- c(
  report_content,
  "\n## 7. 组间统计分析详细结果",
  add_stats_to_report(tcga_stats, "TCGA"),
  add_stats_to_report(oak_stats, "OAK_POPLAR"),
  add_stats_to_report(nt_stats, "NT")
)

# 重新保存报告
# 创建报告子目录
report_dir <- file.path(output_dir, "model_report")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
writeLines(report_content, file.path(report_dir, "model_report.md"))

if (require(rmarkdown)) {
  rmarkdown::render(file.path(report_dir, "model_report.md"),
                    output_format = "html_document",
                    output_file = file.path(report_dir, "model_report.html"))
}

# ================================================================================
#  生成完整的分析报告
# ================================================================================
cat("\n步骤9: 生成分析报告...\n")

# 创建报告目录
report_dir <- file.path(output_dir, "model_report")
if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE)
}

# 生成HTML格式的报告
report_content <- c(
  "# 风险预测模型构建与验证报告",
  "\n## 1. 数据预处理",
  sprintf("- 初始TF regulon数量: %d", length(regulon.list)),
  sprintf("- 总基因数量: %d", length(all_regulon_genes)),
  sprintf("- 表达过滤后基因数量: %d", length(regulon_genes_filtered)),
  sprintf("- 过滤后TF数量: %d", sum(regulon_gene_info_filtered$Type == "TF")),
  sprintf("- 过滤后Target数量: %d", sum(regulon_gene_info_filtered$Type == "Target")),
  
  "\n## 2. 特征选择",
  sprintf("- Lasso初步选择特征数: %d", length(selected_features_initial)),
  sprintf("- 最终选择特征数: %d", length(selected_features)),
  "- 选择的特征:",
  paste("  *", selected_features),
  sprintf("- 最优lambda (1se规则): %.4f", optimal_lambda),
  
  "\n## 3. 模型性能",
  sprintf("- AUC: %.3f", auc_value),
  sprintf("- 准确率: %.3f", accuracy),
  sprintf("- 敏感性: %.3f", sensitivity),
  sprintf("- 特异性: %.3f", specificity),
  sprintf("- 精确度: %.3f", precision),
  
  "\n## 4. 风险评分公式",
  risk_formula,
  
  "\n## 5. 验证集结果",
  "\n### TCGA队列",
  sprintf("- 样本数量: %d", length(tcga_risk_scores)),
  sprintf("- Kruskal-Wallis检验 p值: %.3e", tcga_stats$kruskal$p.value),
  "- 组间倍数变化:",
  paste("  *", names(tcga_stats$fold_changes), ": ", sprintf("%.2f", unlist(tcga_stats$fold_changes))),
  
  "\n### OAK_POPLAR队列",
  sprintf("- 样本数量: %d", length(oak_risk_scores)),
  sprintf("- Kruskal-Wallis检验 p值: %.3e", oak_stats$kruskal$p.value),
  "- 组间倍数变化:",
  paste("  *", names(oak_stats$fold_changes), ": ", sprintf("%.2f", unlist(oak_stats$fold_changes))),
  
  "\n### NT队列",
  sprintf("- 样本数量: %d", nrow(nt_data)),
  sprintf("- Kruskal-Wallis检验 p值: %.3e", nt_stats$kruskal$p.value),
  "- 组间倍数变化:",
  paste("  *", names(nt_stats$fold_changes), ": ", sprintf("%.2f", unlist(nt_stats$fold_changes))),
  
  "\n## 6. 输出文件列表",
  paste("- ", list.files(output_dir, recursive = TRUE))
)

# 在计算风险评分后，添加标准化函数
normalize_risk_scores <- function(scores) {
  # 方法1：Min-Max缩放到[0,1]
  # min_score <- min(scores)
  # max_score <- max(scores)
  # normalized <- (scores - min_score) / (max_score - min_score)
  
  # 方法2：将最小值平移到0，然后所有值加1（确保没有0值）
  min_score <- min(scores)
  if(min_score < 0) {
    normalized <- scores - min_score + 1
  } else {
    normalized <- scores + 1
  }
  
  return(normalized)
}

# 修改风险评分计算部分
# TCGA风险评分计算
tcga_risk_scores <- intercept + as.matrix(t(tcga_expr_regulon_filtered[selected_features, common_samples])) %*% gene_coefs
tcga_risk_scores <- as.vector(tcga_risk_scores)
names(tcga_risk_scores) <- common_samples

# 标准化TCGA风险评分
tcga_risk_scores <- normalize_risk_scores(tcga_risk_scores)

# OAK_POPLAR风险评分计算
available_coefs <- gene_coefs[common_genes_oak]
oak_risk_scores <- intercept + as.matrix(t(oak_poplar_expr[common_genes_oak, common_samples_oak])) %*% available_coefs
oak_risk_scores <- as.vector(oak_risk_scores)
names(oak_risk_scores) <- common_samples_oak

# 标准化OAK_POPLAR风险评分
oak_risk_scores <- normalize_risk_scores(oak_risk_scores)

# NT风险评分计算
available_coefs_nt <- gene_coefs[common_genes_nt]
nt_risk_scores <- intercept + as.matrix(t(nt_expr[common_genes_nt, common_samples_nt])) %*% available_coefs_nt
nt_risk_scores <- as.vector(nt_risk_scores)
names(nt_risk_scores) <- common_samples_nt

# 标准化NT风险评分
nt_risk_scores <- normalize_risk_scores(nt_risk_scores)

# 在保存风险评分结果时添加标准化说明
write.csv(data.frame(
  Sample = names(tcga_risk_scores),
  Original_Risk_Score = tcga_risk_scores_original <- intercept + as.matrix(t(tcga_expr_regulon_filtered[selected_features, common_samples])) %*% gene_coefs,
  Normalized_Risk_Score = tcga_risk_scores,
  Cluster = tcga_clusters[names(tcga_risk_scores), "cluster"]
), file.path(output_dir, "TCGA_risk_scores.csv"), row.names = FALSE)

write.csv(data.frame(
  Sample = names(oak_risk_scores),
  Original_Risk_Score = oak_risk_scores_original <- intercept + as.matrix(t(oak_poplar_expr[common_genes_oak, common_samples_oak])) %*% available_coefs,
  Normalized_Risk_Score = oak_risk_scores,
  Cluster = oak_poplar_clusters[names(oak_risk_scores), "cluster"]
), file.path(output_dir, "OAK_POPLAR_risk_scores.csv"), row.names = FALSE)

write.csv(data.frame(
  Sample = names(nt_risk_scores),
  Original_Risk_Score = nt_risk_scores_original <- intercept + as.matrix(t(nt_expr[common_genes_nt, common_samples_nt])) %*% available_coefs_nt,
  Normalized_Risk_Score = nt_risk_scores,
  Group = nt_metadata$group[match(names(nt_risk_scores), nt_metadata$X)]
), file.path(output_dir, "NT_risk_scores.csv"), row.names = FALSE)

# 在报告中添加标准化说明
report_content <- c(
  report_content,
  "\n## 9. 风险评分标准化说明",
  "为了便于比较和解释，对原始风险评分进行了标准化处理：",
  "1. 将所有评分值平移到正数区间",
  "2. 保持了原始评分的相对关系",
  "3. 标准化后的评分均大于1",
  "\n注意：标准化不影响组间的统计比较结果，因为：",
  "- Wilcoxon检验基于秩和，与线性变换无关",
  "- 组间倍数变化的相对关系保持不变"
)

# 保存报告
# 创建报告子目录
report_dir <- file.path(output_dir, "model_report")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

writeLines(report_content, file.path(report_dir, "model_report.md"))

# 如果安装了rmarkdown，也可以生成HTML版本
if (require(rmarkdown)) {
  rmarkdown::render(file.path(report_dir, "model_report.md"),
                    output_format = "html_document",
                    output_file = file.path(report_dir, "model_report.html"))
}

cat("\n分析报告已生成，保存在:", report_dir, "\n")

# ================================================================================
# 总结和输出
# ================================================================================
cat("\n================================================================================\n")
cat("分析完成！结果总结:\n")
cat("================================================================================\n")

cat("1. TF Regulon基因集结果:\n")
cat("   - TF regulon数量:", length(regulon.list), "个\n")
cat("   - 总基因数量:", length(all_regulon_genes), "个\n")
cat("   - 表达过滤前基因:", length(regulon_genes_in_expr), "个\n")
cat("   - 表达过滤后基因:", length(regulon_genes_filtered), "个\n")
cat("   - 过滤后TF数量:", sum(regulon_gene_info_filtered$Type == "TF"), "个\n")
cat("   - 过滤后Target数量:", sum(regulon_gene_info_filtered$Type == "Target"), "个\n")

cat("\n2. 预测模型结果:\n")
cat("   - Lasso初步选择特征数:", length(selected_features_initial), "个\n")
cat("   - 最终选择特征数:", length(selected_features), "个\n")
cat("   - 使用lambda.1se (更稀疏模型)\n")
cat("   - 模型AUC:", round(auc_value, 3), "\n")
cat("   - 模型准确率:", round(accuracy, 3), "\n")
cat("   - 风险评分公式已保存\n")

cat("\n3. 风险评分结果:\n")
cat("   - TCGA样本:", length(tcga_risk_scores), "个\n")
cat("   - OAK_POPLAR样本:", length(oak_risk_scores), "个\n")
cat("   - NT样本:", length(nt_risk_scores), "个\n")
if (exists("luad_samples") && length(luad_samples) > 0) {
  cat("   - NT LUAD样本:", length(luad_samples), "个\n")
}

cat("\n4. 输出文件:\n")
output_files <- list.files(output_dir, full.names = FALSE)
for (file in output_files) {
  cat("   -", file, "\n")
}

cat("\n所有结果已保存到:", output_dir, "\n")
cat("分析脚本执行完成！\n")

# 清理环境
options(datadist = NULL) 

# 添加异常值处理函数
remove_outliers <- function(data, group_var, value_var, threshold = 1.5) {
  result <- data
  for (group in unique(data[[group_var]])) {
    group_data <- data[[value_var]][data[[group_var]] == group]
    Q1 <- quantile(group_data, 0.25)
    Q3 <- quantile(group_data, 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - threshold * IQR
    upper_bound <- Q3 + threshold * IQR
    
    # 标记异常值
    outliers <- data[[value_var]] > upper_bound & data[[group_var]] == group |
      data[[value_var]] < lower_bound & data[[group_var]] == group
    
    if (any(outliers)) {
      cat(sprintf("在%s组中移除了%d个异常值\n", group, sum(outliers)))
      result <- result[!outliers, ]
    }
  }
  return(result)
}

# 处理NT数据中的异常值
nt_data_clean <- remove_outliers(nt_data, "Group", "Risk_Score", threshold = 2)

# 更新NT数据的统计分析
nt_stats <- perform_statistical_tests(nt_data_clean, "Group", "Risk_Score")
cat("\nNT组间统计分析结果 (去除异常值后):\n")
cat("Kruskal-Wallis检验 p值:", nt_stats$kruskal$p.value, "\n")
cat("组间倍数变化:\n")
print(unlist(nt_stats$fold_changes))
cat("Wilcoxon检验 p值:\n")
print(unlist(nt_stats$pairwise))

# 更新NT小提琴箱线图
p_nt <- create_violin_boxplot(nt_data_clean, "Group", "Risk_Score",
                              "NT Risk Score Distribution by Group\n(Outliers Removed)",
                              nt_colors)
ggsave(file.path(viz_dir, "NT_risk_violin_box.pdf"), p_nt,
       width = 5, height = 8)

# 更新组合图
combined_plot <- ggpubr::ggarrange(p_tcga, p_oak, p_nt,
                                   ncol = 3, nrow = 1,
                                   labels = c("A", "B", "C"))

ggsave(file.path(viz_dir, "combined_risk_violin_box.pdf"), combined_plot,
       width = 15, height = 8)  # 调整组合图高度

# 保存清理后的数据
write.csv(nt_data_clean, file.path(viz_dir, "NT_risk_scores_clean.csv"), row.names = FALSE)

# 更新报告中的NT结果部分
report_content <- gsub(
  "\n### NT队列",
  sprintf("\n### NT队列 (异常值处理后)\n- 原始样本数量: %d\n- 清理后样本数量: %d", 
          nrow(nt_data), nrow(nt_data_clean)),
  report_content
) 


###############################################################################
# 追加：绘制 3 条 ROC 曲线（同一坐标系，ggplot 版）
###############################################################################

library(pROC)
library(ggplot2)
library(dplyr)

# 颜色向量（顺序随意，只要名字匹配）
roc_cols <- c(
  TCGA = "#3F51B5",
  OAK  = "#ED0000",
  NT   = "#FF8C00"
)

# ------------------------------------------------------------------
# 1. 构造标签并计算 ROC
# ------------------------------------------------------------------
make_roc_df <- function(score, label, set_name) {
  roc_obj <- roc(label, score, quiet = TRUE)
  coords  <- coords(roc_obj, "all", ret = c("specificity", "sensitivity"))
  data.frame(
    FPR = 1 - coords$specificity,
    TPR = coords$sensitivity,
    Dataset = set_name
  )
}

tcga_label <- factor(ifelse(tcga_data$Group == "cluster_high", "high_risk", "low_risk"),
                     levels = c("low_risk", "high_risk"))
oak_label  <- factor(ifelse(oak_data$Group == "cluster_high", "high_risk", "low_risk"),
                     levels = c("low_risk", "high_risk"))
nt_label   <- factor(ifelse(nt_data_clean$Group == "C_LUAD", "low_risk", "high_risk"),
                     levels = c("low_risk", "high_risk"))

roc_df <- bind_rows(
  make_roc_df(tcga_data$Risk_Score, tcga_label, "TCGA"),
  make_roc_df(oak_data$Risk_Score,  oak_label,  "OAK"),
  make_roc_df(nt_data_clean$Risk_Score, nt_label, "NT")
)

# ------------------------------------------------------------------
# 2. 绘图
# ------------------------------------------------------------------
p <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Dataset)) +
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = roc_cols) +
  labs(title = "Combined ROC Curves",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = .5, face = "bold"),
        legend.title = element_blank())

# ------------------------------------------------------------------
# 3. 在图上标注 AUC
# ------------------------------------------------------------------
auc_vals <- c(
  TCGA = round(auc(roc(tcga_label, tcga_data$Risk_Score, quiet = TRUE)), 3),
  OAK  = round(auc(roc(oak_label,  oak_data$Risk_Score,  quiet = TRUE)), 3),
  NT   = round(auc(roc(nt_label,   nt_data_clean$Risk_Score, quiet = TRUE)), 3)
)

p <- p + annotate(
  "text",
  x = 0.65, y = 0.35,
  label = paste(sprintf("%s AUC = %.3f", names(auc_vals), auc_vals), collapse = "\n"),
  size = 4.5, hjust = 0
)

# ------------------------------------------------------------------
# 4. 保存
# ------------------------------------------------------------------
roc_dir <- file.path(output_dir, "roc_curves")
dir.create(roc_dir, showWarnings = FALSE)
ggsave(file.path(roc_dir, "ROC_combined_single.pdf"), p, width = 7, height = 6)

cat("\nCombined ROC 曲线已保存至:", file.path(roc_dir, "ROC_combined_single.pdf"), "\n")

