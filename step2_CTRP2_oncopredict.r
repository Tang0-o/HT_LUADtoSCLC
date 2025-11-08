# ============================================================================
# 0. 加载库 & 创建输出文件夹
# ============================================================================
library(oncoPredict)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(stringr)          # 新增：用于 CVCL 提取

#这里需要设置为觉得路径，不然后续oncopredict生成的文件会默认在根目录下，新建文件夹。
output_dir <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/2_DrugPredictions_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ============================================================================
# 1. 数据路径定义
# ============================================================================
expression_data_path <- "/home/data/tmh_project/SCLC/Fig3_multicohort_v1/0_Fig3_cohort_data/1_TCGA_LUAD/LUAD_fpkm_new.csv"
risk_score_data_path  <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/1_risk_prediction/TCGA_risk_scores.csv"
ctrp_expr_path        <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/osf.io_trainingdata/CTRP2_Expr (RPKM, not log transformed).rds"
ctrp_res_path         <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/osf.io_trainingdata/CTRP2_Res.rds"

# ============================================================================
# 2. 数据加载
# ============================================================================
print("Loading LUAD expression data...")
luad_expr <- read.csv(expression_data_path, row.names = 1)
print(paste("LUAD expression data dimensions:", nrow(luad_expr), "x", ncol(luad_expr)))

print("Loading risk score data...")
risk_scores <- read.csv(risk_score_data_path, row.names = 1)
print(paste("Risk score data dimensions:", nrow(risk_scores), "x", ncol(risk_scores)))

print("Loading CTRP training data...")
CTRP2_Expr <- readRDS(ctrp_expr_path)
CTRP2_Res  <- readRDS(ctrp_res_path)
print(paste("CTRP2 expression dimensions:", nrow(CTRP2_Expr), "x", ncol(CTRP2_Expr)))
print(paste("CTRP2 response dimensions:", nrow(CTRP2_Res), "x", ncol(CTRP2_Res)))

# ============================================================================
# 3. CTRP 样本名称精确匹配（CVCL-ID 对齐）
# ============================================================================
cat("\n=== 创建细胞系映射（CVCL-ID 对齐） ===\n")

# 1) 从表达矩阵列名提取 CVCL_XXXX
tmp     <- str_extract(colnames(CTRP2_Expr), "(?<=CVCL_)\\d+")
expr_id <- paste0("CVCL_", tmp)

# 2) 响应矩阵行名已经是 CVCL_XXXX
res_id <- rownames(CTRP2_Res)

# 3) 交集
common_id <- intersect(expr_id, res_id)
cat("成功匹配细胞系数:", length(common_id), "\n")
if(length(common_id) < 100){
  stop("匹配样本过少，请检查 CTRP 数据版本或命名格式！")
}

# 4) 对齐两个矩阵
expr_idx <- match(common_id, expr_id)
res_idx  <- match(common_id, res_id)

ctrp_matched <- list(
  expr = CTRP2_Expr[, expr_idx] |> `colnames<-`(common_id),
  res  = CTRP2_Res[res_idx, ]   |> `rownames<-`(common_id),
  method = "CVCL_ID_matched"
)

# ============================================================================
# 4. LUAD 数据预处理
# ============================================================================
preprocess_luad_data <- function(expr_data) {
  tumor_samples <- grep("01A$", colnames(expr_data), value = TRUE)
  if(length(tumor_samples) == 0) tumor_samples <- grep("01$", colnames(expr_data), value = TRUE)
  if(length(tumor_samples) == 0) {
    cat("Warning: No tumor samples found, using all samples\n")
    tumor_samples <- colnames(expr_data)
  }
  expr_tumor <- expr_data[, tumor_samples, drop = FALSE]
  expr_tumor <- expr_tumor[!duplicated(rownames(expr_tumor)), ]
  gene_means <- rowMeans(expr_tumor, na.rm = TRUE)
  expr_filtered <- expr_tumor[gene_means > 1, ]
  print(paste("After filtering:", nrow(expr_filtered), "genes,", ncol(expr_filtered), "samples"))
  return(expr_filtered)
}

luad_expr_processed <- preprocess_luad_data(luad_expr)

# ============================================================================
# 5. 药物敏感性预测（指定输出目录）
# ============================================================================
outDir <- file.path(output_dir, "CTRP2_calcPhenotype_Output")
dir.create(outDir)
owd <- setwd(outDir)

print("\n=== 运行oncoPredict进行药物敏感性预测（CTRP 训练集） ===")
tryCatch({
  drug_predictions <- calcPhenotype(
    trainingExprData = ctrp_matched$expr,
    trainingPtype    = ctrp_matched$res,
    testExprData     = as.matrix(luad_expr_processed),
    batchCorrect     = 'eb',
    minNumSamples    = 10,
    printOutput      = TRUE,
    removeLowVaryingGenes = 0.1,
    removeLowVaringGenesFrom = "homogenizeData"
  )
  cat("药物敏感性预测成功！\n")
  print(paste("预测结果维度:", nrow(drug_predictions), "x", ncol(drug_predictions)))
}, error = function(e) {
  cat("批次校正失败，尝试不使用批次校正...\n")
  drug_predictions <<- calcPhenotype(
    trainingExprData = ctrp_matched$expr,
    trainingPtype    = ctrp_matched$res,
    testExprData     = as.matrix(luad_expr_processed),
    batchCorrect     = 'none',
    minNumSamples    = 5,
    printOutput      = TRUE,
    removeLowVaryingGenes = 0.05,
    removeLowVaringGenesFrom = "rawData"
  )
  cat("药物敏感性预测完成（无批次校正）！\n")
  print(paste("预测结果维度:", nrow(drug_predictions), "x", ncol(drug_predictions)))
})

setwd(owd)