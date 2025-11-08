# ============================================================================
# 基于CTRP2预测结果的药物敏感性可视化分析代码（增强版）
# 适配TCGA_LUAD + CTRP2预测数据 + 总览报告生成
# ============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(readr)

# 设置文件路径
drug_predictions_path <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/2_DrugPredictions_Output/CTRP2_calcPhenotype_Output/CTRP2_DrugPredictions.csv"
risk_score_path <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/1_risk_prediction/TCGA_risk_scores.csv"

# 设置输出目录
output_dir <- "/home/data/tmh_project/SCLC/Fig5_Risk_prediction_model/2_DrugPredictions_Output/CTRP2_drug_analysis_v3"
correlation_dir <- file.path(output_dir, "risk_ic50_correlations")
data_dir <- file.path(output_dir, "correlation_data")

# 创建所需目录
dirs_to_create <- c(
  output_dir, 
  correlation_dir, 
  data_dir,
  file.path(correlation_dir, "resistant_in_high_risk"),
  file.path(correlation_dir, "sensitive_in_high_risk"),
  file.path(correlation_dir, "weak_correlation")
)

for(dir in dirs_to_create) {
  if(!dir.exists(dir)) {
    success <- dir.create(dir, recursive = TRUE)
    cat("创建目录", dir, ":", ifelse(success, "成功", "失败"), "\n")
  } else {
    cat("目录已存在:", dir, "\n")
  }
}

# ============================================================================
# 1. 数据读取和预处理函数
# ============================================================================

load_and_prepare_ctrp2_data <- function() {
  cat("=== 读取CTRP2预测数据 ===\n")
  
  # 读取药物预测结果
  if(!file.exists(drug_predictions_path)) {
    cat("错误: 药物预测文件不存在:", drug_predictions_path, "\n")
    return(NULL)
  }
  
  cat("读取药物预测数据...\n")
  drug_predictions <- read.csv(drug_predictions_path, row.names = 1)
  cat("药物预测数据维度:", dim(drug_predictions), "\n")
  cat("样本数:", nrow(drug_predictions), "\n")
  cat("药物数:", ncol(drug_predictions), "\n")
  
  # 读取风险评分数据
  if(!file.exists(risk_score_path)) {
    cat("错误: 风险评分文件不存在:", risk_score_path, "\n")
    return(NULL)
  }
  
  cat("读取风险评分数据...\n")
  risk_scores <- read.csv(risk_score_path, row.names = 1)
  cat("风险评分数据维度:", dim(risk_scores), "\n")
  
  # 检查数据结构
  cat("风险评分列名:", colnames(risk_scores), "\n")
  cat("前几个样本ID:", head(rownames(risk_scores)), "\n")
  
  # 找到共同样本
  common_samples <- intersect(rownames(drug_predictions), rownames(risk_scores))
  cat("共同样本数:", length(common_samples), "\n")
  
  if(length(common_samples) < 50) {
    cat("警告: 共同样本数太少 (", length(common_samples), ")，可能影响分析质量\n")
  }
  
  # 合并数据
  drug_data_subset <- drug_predictions[common_samples, ]
  risk_data_subset <- risk_scores[common_samples, ]
  
  # 创建组合数据
  combined_data <- data.frame(
    Sample_ID = common_samples,
    Risk_Score = risk_data_subset[, 1], # 假设风险评分在第一列
    stringsAsFactors = FALSE
  )
  
  # 添加药物数据
  combined_data <- cbind(combined_data, drug_data_subset)
  
  # 创建风险分组（使用中位数分割）
  risk_median <- median(combined_data$Risk_Score, na.rm = TRUE)
  combined_data$Risk_Group <- ifelse(combined_data$Risk_Score >= risk_median, "High", "Low")
  
  cat("风险分组分布:\n")
  print(table(combined_data$Risk_Group))
  
  # 检查数据质量
  cat("检查数据质量...\n")
  cat("Risk_Score缺失值:", sum(is.na(combined_data$Risk_Score)), "\n")
  cat("药物数据缺失值比例:", 
      round(mean(is.na(combined_data[, 4:ncol(combined_data)])), 3), "\n")
  
  return(combined_data)
}

# ============================================================================
# 2. 药物分类函数
# ============================================================================
create_ctrp2_drug_categories <- function(drug_names) {
  cat("=== 创建CTRP2药物分类 ===\n")
  
  # ===================================================================
  # 先计算 LUAD 和 SCLC 一线用药（放在 list 外面）
  # ===================================================================
    # 【LUAD 一线用药计算】
    cat("识别 LUAD 一线用药...\n")
    
    # 1️⃣ 匹配 EGFR-TKI 单药
    egfr_tki <- grep(
      "^(osimertinib|gefitinib|erlotinib|afatinib)($|[^a-z])",
      drug_names, ignore.case = TRUE, value = TRUE
    )
    # 排除与研究药物的组合
    egfr_tki <- egfr_tki[!grepl(
      "(PLX|UNC|navitoclax|decitabine|selumetinib|trametinib|vemurafenib|BAY|GSK|AZD|SB)", 
      egfr_tki, ignore.case = TRUE
    )]
    
    # 2️⃣ 匹配 Pemetrexed + Platinum 组合（扩展分隔符支持"."、"+"、" "等）
    pem_plat <- grep(
      "(pemetrexed[\\. +]*(carboplatin|cisplatin))|((carboplatin|cisplatin)[\\. +]*pemetrexed)",
      drug_names, ignore.case = TRUE, value = TRUE
    )
    
    # 3️⃣ 匹配 Gemcitabine + Platinum 组合（类似扩展）
    gem_plat <- grep(
      "(gemcitabine[\\. +]*(carboplatin|cisplatin))|((carboplatin|cisplatin)[\\. +]*gemcitabine)",
      drug_names, ignore.case = TRUE, value = TRUE
    )
    
    # 4️⃣ 匹配单药化疗（包括铂类，但排除与 etoposide 的组合）
    single_chemo <- grep(
      "^(pemetrexed|carboplatin|cisplatin|gemcitabine)($|[^a-z])",
      drug_names, ignore.case = TRUE, value = TRUE
    )
    single_chemo <- single_chemo[!grepl("etoposide", single_chemo, ignore.case = TRUE)]
    
    # 新增：如果单药未匹配，从组合中提取铂类单药（假设数据只有组合）
    if (length(grep("^(carboplatin|cisplatin)$", single_chemo, ignore.case = TRUE)) == 0) {
      cat("  未找到纯单药carboplatin/cisplatin，从组合中提取...\n")
      plat_combos <- grep("(carboplatin|cisplatin)", drug_names, ignore.case = TRUE, value = TRUE)
      # 提取唯一铂类单药名称（模拟添加）
      extracted_plat <- unique(gsub(".*(carboplatin|cisplatin).*", "\\1", plat_combos, ignore.case = TRUE))
      single_chemo <- unique(c(single_chemo, extracted_plat[grepl("^(carboplatin|cisplatin)$", extracted_plat, ignore.case = TRUE)]))
    }
    
    # 5️⃣ 排除含研究性药物的组合
    experimental_pattern <- "(PLX|UNC|GSK|AZD|SB|navitoclax|decitabine|selumetinib|trametinib|vorinostat|BAY|CIL|FQI|triazolothiadiazine)"
    pem_plat <- pem_plat[!grepl(experimental_pattern, pem_plat, ignore.case = TRUE)]
    gem_plat <- gem_plat[!grepl(experimental_pattern, gem_plat, ignore.case = TRUE)]
    single_chemo <- single_chemo[!grepl(experimental_pattern, single_chemo, ignore.case = TRUE)]
    
    # 6️⃣ 合并所有 LUAD 药物
    LUAD_first_line <- unique(c(egfr_tki, pem_plat, gem_plat, single_chemo))
    
    cat("  找到", length(LUAD_first_line), "个 LUAD 一线相关药物\n")
    if(length(LUAD_first_line) <= 10 && length(LUAD_first_line) > 0) {
      cat("  药物列表:", paste(LUAD_first_line, collapse = ", "), "\n")
    } else if(length(LUAD_first_line) > 10) {
      cat("  药物列表示例 (前10个):", paste(head(LUAD_first_line, 10), collapse = ", "), "...\n")
      cat("  包含配伍方案示例: 培美曲塞 + 铂类 (如 pemetrexed + cisplatin), 吉西他滨 + 铂类 (如 gemcitabine + carboplatin)\n")
      cat("  铂类单药匹配:", paste(grep("(carboplatin|cisplatin)", LUAD_first_line, value = TRUE), collapse = ", "), "\n")
    }
  
  # 【SCLC 一线用药计算】
  cat("识别 SCLC 一线用药...\n")
  
  # 1. 匹配 carboplatin/cisplatin + etoposide 组合
  plat_etop <- grep(
    "(carboplatin|cisplatin).*etoposide|etoposide.*(carboplatin|cisplatin)",
    drug_names, ignore.case = TRUE, value = TRUE
  )
  
  # 2. 匹配 irinotecan + platinum 组合
  iri_plat <- grep(
    "(irinotecan.*(carboplatin|cisplatin))|((carboplatin|cisplatin).*irinotecan)",
    drug_names, ignore.case = TRUE, value = TRUE
  )
  
  # 3. 匹配单药
  single_sclc <- grep(
    "^(etoposide|paclitaxel|irinotecan)($|[^a-z])",
    drug_names, ignore.case = TRUE, value = TRUE
  )
  
  # 合并所有 SCLC 药物
  SCLC_first_line <- unique(c(plat_etop, iri_plat, single_sclc))
  
  cat("  找到", length(SCLC_first_line), "个 SCLC 一线相关药物\n")
  if(length(SCLC_first_line) <= 10 && length(SCLC_first_line) > 0) {
    cat("  药物列表:", paste(SCLC_first_line, collapse = ", "), "\n")
  }
  
  # ===================================================================
  # 定义其他药物类别
  # ===================================================================
  
  drug_categories <- list(
    
    # 【临床一线用药分类 - 直接使用上面计算的结果】
    LUAD_first_line_drugs = LUAD_first_line,
    SCLC_first_line_drugs = SCLC_first_line,
    
    # 靶向治疗药物
    EGFR_inhibitors = grep("(gefitinib|erlotinib|afatinib|osimertinib)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    ALK_inhibitors = grep("(crizotinib|alectinib|ceritinib|brigatinib)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    PI3K_inhibitors = grep("(idelalisib|copanlisib|duvelisib|alpelisib|umbralisib|pictilisib|buparlisib|inavolisib|taselisib|zstk474|gdc-0941|gdc0941|gedatolisib|paxalisib|bimiralisib|samotolisib|voxtalisib|bkm120|mln1117|apitolisib)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    AKT_inhibitors = grep("(capivasertib|ipatasertib|afuresertib|mk-2206|mk2206|perifosine|tas-117|tas117|gsk2141795|m2698|azd5363|gsk2110183|uprosertib|azd8186|azd8189|arq092|arq-092|miransertib|azd8055)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    mTOR_inhibitors = grep("(sirolimus|rapamycin|temsirolimus|everolimus|ridaforolimus)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    MEK_inhibitors = grep("(mek|trametinib|cobimetinib|selumetinib)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    MYC_inhibitors = grep("(omomyc|omo-103|omo103|wbc100|wbc-100|mrt-2359|mrt2359|idp-121|idp121|kb-0742|kb0742|otx-2002|otx2002|pc-002|pc002|sepantronium)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    CDK_inhibitors = grep("(cdk|palbociclib|ribociclib|abemaciclib)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    PARP_inhibitors = grep("(parp|olaparib|rucaparib|niraparib|talazoparib)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    # 药物目录中没有
    EZH2_inhibitors = grep("(tazemetostat|valemetostat|GSK126|CPI-1205|PF-06821497)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    # 化疗药物
    Platinum_agents = grep("(cisplatin|carboplatin|oxaliplatin)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    Taxanes = grep("(paclitaxel|docetaxel|cabazixtaxel)",
                   drug_names, ignore.case = TRUE, value = TRUE),
    
    Topoisomerase_inhibitors = grep("(topotecan|irinotecan|etoposide|doxorubicin)",
                                    drug_names, ignore.case = TRUE, value = TRUE),
    
    Antimetabolites = grep("(gemcitabine|pemetrexed|5-fu|capecitabine|methotrexate)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    # Checkpoint 抑制剂
    # 药物目录中没有
    PD1_PDL1_inhibitors = grep("(nivolumab|pembrolizumab|cemiplimab|sintilimab|toripalimab|camrelizumab|dostarlimab|retifanlimab|tislelizumab|envafolimab|kn035|atezolizumab|durvalumab|avelumab|sugemalimab|cs1001|duligotuzumab|mehd7945a)",
                               drug_names, ignore.case = TRUE, value = TRUE),
    
    PD1_inhibitors = grep("(nivolumab|pembrolizumab|cemiplimab|sintilimab|toripalimab|camrelizumab|dostarlimab|retifanlimab|tislelizumab|envafolimab|kn035)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    PDL1_inhibitors = grep("(atezolizumab|durvalumab|avelumab|sugemalimab|cs1001|duligotuzumab|mehd7945a)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    checkpoint_inhibitors = grep("(nivolumab|pembrolizumab|cemiplimab|sintilimab|toripalimab|camrelizumab|dostarlimab|retifanlimab|tislelizumab|envafolimab|kn035|atezolizumab|durvalumab|avelumab|sugemalimab|cs1001|duligotuzumab|mehd7945a|ipilimumab|tremelimumab|relatlimab|tiragolumab|mtig7192a|vibostolimab|mk-7684|mk7684|mbg453|sabstitimab|kn046|m7824|bintrafusp|mcla-128|mcla128)",
                                 drug_names, ignore.case = TRUE, value = TRUE),
    
    # 新兴靶向药物
    Aurora_kinase_inhibitors = grep("(aurora|alisertib|barasertib)",
                                    drug_names, ignore.case = TRUE, value = TRUE),
    
    PLK_inhibitors = grep("(plk|volasertib)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    WEE1_inhibitors = grep("(wee1|adavosertib)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    CHK_inhibitors = grep("(chk|prexasertib)",
                          drug_names, ignore.case = TRUE, value = TRUE),
    
    HDAC_inhibitors = grep("(hdac|dasatinib|vorinostat|romidepsin|panobinostat)",
                           drug_names, ignore.case = TRUE, value = TRUE),
    
    DNMT_inhibitors = grep("(azacitidine|5-azacytidine|decitabine|5-aza-2-deoxycytidine|guadecitabine|sgn-110|zebularine|rg108|rg-108|hydralazine|procainamide|procaine|egcg|epigallocatechin|mg98)",
                           drug_names, ignore.case = TRUE, value = TRUE)
  )
  
  # 移除空类别
  drug_categories <- drug_categories[sapply(drug_categories, length) > 0]
  
  # 创建"All_CTRP2_Drugs"类别
  drug_categories$All_CTRP2_Drugs <- drug_names
  
  # 打印分类结果
  cat("\nCTRP2药物分类结果:\n")
  for(category in names(drug_categories)) {
    cat("-", category, ":", length(drug_categories[[category]]), "个药物\n")
    if(length(drug_categories[[category]]) <= 5 && length(drug_categories[[category]]) > 0) {
      cat("  ", paste(drug_categories[[category]], collapse = ", "), "\n")
    } else if(length(drug_categories[[category]]) > 5) {
      cat("  ", paste(head(drug_categories[[category]], 3), collapse = ", "), "...\n")
    }
  }
  
  return(list(major_categories = drug_categories))
}


# ============================================================================
# 3. 计算相关系数函数
# ============================================================================

calculate_drug_correlations <- function(data, drug_list) {
  correlation_results <- data.frame(
    Drug = character(),
    Spearman_r = numeric(),
    Spearman_p = numeric(),
    Pearson_r = numeric(),
    N_samples = integer(),
    High_mean = numeric(),
    Low_mean = numeric(),
    Mean_diff = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(drug in drug_list) {
    if(!drug %in% colnames(data)) next
    
    tryCatch({
      drug_data <- data[, drug]
      risk_score <- data$Risk_Score
      risk_group <- data$Risk_Group
      
      # 移除缺失值
      valid_idx <- !is.na(drug_data) & !is.na(risk_score) & !is.na(risk_group)
      
      if(sum(valid_idx) < 15) next
      
      # 计算相关系数
      spearman_test <- cor.test(risk_score[valid_idx], drug_data[valid_idx], method = "spearman")
      pearson_r <- cor(risk_score[valid_idx], drug_data[valid_idx], use = "complete.obs")
      
      # 计算组间均值
      high_vals <- drug_data[valid_idx & risk_group == "High"]
      low_vals <- drug_data[valid_idx & risk_group == "Low"]
      
      if(length(high_vals) >= 5 && length(low_vals) >= 5) {
        correlation_results <- rbind(correlation_results, data.frame(
          Drug = drug,
          Spearman_r = as.numeric(spearman_test$estimate),
          Spearman_p = spearman_test$p.value,
          Pearson_r = pearson_r,
          N_samples = sum(valid_idx),
          High_mean = mean(high_vals),
          Low_mean = mean(low_vals),
          Mean_diff = mean(high_vals) - mean(low_vals),
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      cat("计算", drug, "相关系数时出错:", e$message, "\n")
    })
  }
  
  return(correlation_results)
}

# ============================================================================
# 4. 条形图绘制函数
# ============================================================================

plot_correlation_barplot <- function(correlation_results, category_name, 
                                     output_dir = "CTRP2_drug_analysis") {
  
  if(nrow(correlation_results) == 0) {
    cat("警告:", category_name, "没有相关性结果\n")
    return(NULL)
  }
  
  tryCatch({
    # 筛选显著结果
    sig_results <- correlation_results[correlation_results$Spearman_p < 0.05, ]
    
    if(nrow(sig_results) == 0) {
      cat("警告:", category_name, "没有显著的相关性结果\n")
      return(NULL)
    }
    
    # 按相关系数排序
    sig_results <- sig_results[order(sig_results$Spearman_r), ]
    
    # 限制显示药物数量（最多20个）
    if(nrow(sig_results) > 20) {
      sig_results <- sig_results[order(abs(sig_results$Spearman_r), decreasing = TRUE)[1:20], ]
      sig_results <- sig_results[order(sig_results$Spearman_r), ]
    }
    
    # 添加分类标签
    sig_results$Direction <- ifelse(sig_results$Spearman_r > 0, "Resistant_in_High_Risk", "Sensitive_in_High_Risk")
    sig_results$Significance <- ifelse(sig_results$Spearman_p < 0.001, "***",
                                       ifelse(sig_results$Spearman_p < 0.01, "**", "*"))
    
    # 创建条形图
    p <- ggplot(sig_results, aes(x = reorder(Drug, Spearman_r), y = Spearman_r)) +
      geom_col(aes(fill = Direction), alpha = 0.8, width = 0.7) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
      geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50", alpha = 0.7) +
      scale_fill_manual(values = c("Resistant_in_High_Risk" = "#DC143C", 
                                   "Sensitive_in_High_Risk" = "#2E8B57")) +
      coord_flip() +
      labs(
        title = paste("CTRP2 Drug Sensitivity: Risk Score Correlation"),
        subtitle = paste(category_name, "- Significant correlations (p < 0.05)"),
        x = "Drug",
        y = "Spearman Correlation Coefficient (r)",
        fill = "Effect Direction"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # 添加显著性标记
    p <- p + geom_text(aes(label = Significance, 
                           x = Drug, 
                           y = ifelse(Spearman_r > 0, Spearman_r + 0.02, Spearman_r - 0.02)),
                       size = 3, hjust = ifelse(sig_results$Spearman_r > 0, 0, 1))
    
    # 保存条形图
    filename <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", category_name), "_correlation_barplot.pdf"))
    ggsave(filename, plot = p, width = 6, height = max(3, nrow(sig_results) * 0.3))
    
    # 保存数据
    data_filename <- file.path(data_dir, paste0(gsub("[^A-Za-z0-9]", "_", category_name), "_correlation_barplot_data.csv"))
    write.csv(sig_results, data_filename, row.names = FALSE)
    
    cat("保存相关系数条形图:", filename, "\n")
    cat("保存条形图数据:", data_filename, "\n")
    
    return(p)
    
  }, error = function(e) {
    cat("条形图生成失败:", category_name, "- 错误:", e$message, "\n")
    return(NULL)
  })
}

# ============================================================================
# 5. 药物敏感性模式判断函数
# ============================================================================

judge_drug_sensitivity_pattern <- function(plot_data, drug_name) {
  
  # 计算Spearman相关系数
  cor_test <- cor.test(plot_data$Risk_Score, plot_data$IC50, method = "spearman")
  cor_coef <- as.numeric(cor_test$estimate)
  cor_pval <- cor_test$p.value
  
  # 计算Pearson相关系数作为补充
  cor_pearson <- cor(plot_data$Risk_Score, plot_data$IC50, use = "complete.obs")
  
  # 组间均值比较
  group_means <- tapply(plot_data$IC50, plot_data$Risk_Group, mean, na.rm = TRUE)
  mean_diff <- group_means["High"] - group_means["Low"]
  
  # 综合判断逻辑
  direction <- "Unknown"
  confidence <- "Low"
  folder_name <- "weak_correlation"
  
  # 主要判断：相关性为主，均值差异为辅
  if(!is.na(cor_pval) && cor_pval < 0.05) {
    if(cor_coef >= 0.2) {
      direction <- "Resistant_in_High_Risk"
      folder_name <- "resistant_in_high_risk"
      confidence <- ifelse(cor_coef >= 0.4, "High", 
                           ifelse(cor_coef >= 0.3, "Medium", "Low"))
    } else if(cor_coef <= -0.2) {
      direction <- "Sensitive_in_High_Risk"
      folder_name <- "sensitive_in_high_risk"
      confidence <- ifelse(cor_coef <= -0.4, "High", 
                           ifelse(cor_coef <= -0.3, "Medium", "Low"))
    } else if(abs(cor_coef) >= 0.1) {
      direction <- if(cor_coef > 0) "Resistant_in_High_Risk" else "Sensitive_in_High_Risk"
      folder_name <- if(cor_coef > 0) "resistant_in_high_risk" else "sensitive_in_high_risk"
      confidence <- "Low"
    }
  } else {
    if(abs(mean_diff) > 0.5 && !is.na(mean_diff)) {
      direction <- if(mean_diff > 0) "Resistant_in_High_Risk" else "Sensitive_in_High_Risk"
      folder_name <- if(mean_diff > 0) "resistant_in_high_risk" else "sensitive_in_high_risk"
      confidence <- "Very_Low"
    }
  }
  
  return(list(
    direction = direction,
    confidence = confidence,
    folder_name = folder_name,
    spearman_r = cor_coef,
    spearman_p = cor_pval,
    pearson_r = cor_pearson,
    mean_difference = mean_diff,
    high_mean = group_means["High"],
    low_mean = group_means["Low"]
  ))
}

# ============================================================================
# 6. 相关性散点图绘制函数
# ============================================================================

plot_risk_ic50_correlation_enhanced <- function(data, drug_list, category_name, 
                                                correlation_results, 
                                                correlation_dir, data_dir) {
  
  cat("进入相关性分析函数 - 类别:", category_name, "\n")
  
  if(length(drug_list) == 0) {
    cat("警告: 药物列表为空\n")
    return(NULL)
  }
  
  if(nrow(correlation_results) == 0) {
    cat("警告: 该类别无相关性结果\n")
    return(NULL)
  }
  
  # 获取显著药物
  significant_drugs <- correlation_results$Drug[correlation_results$Spearman_p < 0.05]
  cat("显著药物数量:", length(significant_drugs), "\n")
  
  if(length(significant_drugs) == 0) {
    cat("警告: 无显著相关性药物\n")
    return(NULL)
  }
  
  # 筛选可用的显著药物
  available_sig_drugs <- significant_drugs[significant_drugs %in% colnames(data)]
  cat("数据中可用的显著药物数量:", length(available_sig_drugs), "\n")
  
  if(length(available_sig_drugs) == 0) {
    cat("警告: 无可用的显著药物进行相关性分析\n")
    return(NULL)
  }
  
  # 存储所有药物的判断结果
  all_judgments <- data.frame(
    Category = character(),
    Drug = character(),
    Direction = character(),
    Confidence = character(),
    Spearman_r = numeric(),
    Spearman_p = numeric(),
    Pearson_r = numeric(),
    Mean_Diff = numeric(),
    High_Mean = numeric(),
    Low_Mean = numeric(),
    N_samples = integer(),
    stringsAsFactors = FALSE
  )
  
  # 为每个显著药物进行分析
  for(drug in available_sig_drugs) {
    cat("  正在分析药物:", drug, "\n")
    
    tryCatch({
      # 提取数据
      drug_data <- data[, drug]
      risk_score <- data$Risk_Score
      risk_group <- data$Risk_Group
      
      # 移除缺失值
      valid_idx <- !is.na(drug_data) & !is.na(risk_score) & !is.na(risk_group)
      
      if(sum(valid_idx) < 15) {
        cat("    药物", drug, "有效数据点不足（", sum(valid_idx), "< 15），跳过\n")
        next
      }
      
      plot_data <- data.frame(
        Risk_Score = risk_score[valid_idx],
        IC50 = drug_data[valid_idx],
        Risk_Group = risk_group[valid_idx]
      )
      
      # 使用判断函数
      judgment <- judge_drug_sensitivity_pattern(plot_data, drug)
      
      # 记录判断结果
      all_judgments <- rbind(all_judgments, data.frame(
        Category = category_name,
        Drug = drug,
        Direction = judgment$direction,
        Confidence = judgment$confidence,
        Spearman_r = judgment$spearman_r,
        Spearman_p = judgment$spearman_p,
        Pearson_r = judgment$pearson_r,
        Mean_Diff = judgment$mean_difference,
        High_Mean = judgment$high_mean,
        Low_Mean = judgment$low_mean,
        N_samples = nrow(plot_data)
      ))
      
      # 创建散点图（修改标签以适应CTRP2数据）
      p <- ggplot(plot_data, aes(x = Risk_Score, y = IC50)) +
        geom_point(aes(color = Risk_Group), alpha = 0.6, size = 1.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.6) +
        scale_color_manual(values = c("Low" = "#2E8B57", "High" = "#DC143C")) +
        labs(
          title = paste("Risk Score vs", drug, "Sensitivity (CTRP2)"),
          subtitle = paste0(
            category_name, " | ", judgment$direction, " (", judgment$confidence, " confidence)",
            "\nSpearman r = ", round(judgment$spearman_r, 3), 
            ", p = ", format(judgment$spearman_p, scientific = TRUE, digits = 3),
            " | Mean diff = ", round(judgment$mean_difference, 3)
          ),
          x = "Histological Transformation Risk Score",
          y = paste(drug, "Predicted Sensitivity (CTRP2)"),
          color = "Risk Group"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8)
        )
      
      # 添加样本量信息
      p <- p + annotate("text", 
                        x = min(plot_data$Risk_Score) + 0.1 * diff(range(plot_data$Risk_Score)),
                        y = max(plot_data$IC50) - 0.1 * diff(range(plot_data$IC50)),
                        label = paste("n =", nrow(plot_data)),
                        hjust = 0, vjust = 1, size = 3, color = "black")
      
      # 根据判断结果确定保存目录
      save_dir <- file.path(correlation_dir, judgment$folder_name)
      
      # 保存散点图
      plot_filename <- file.path(save_dir, paste0(category_name, "_", drug, "_correlation.pdf"))
      ggsave(plot_filename, plot = p, width = 6, height = 4.5, dpi = 300)
      
      # 保存对应的数据文件
      data_filename <- file.path(save_dir, paste0(category_name, "_", drug, "_plot_data.csv"))
      write.csv(plot_data, data_filename, row.names = FALSE)
      
      # 保存详细的统计信息
      stats_data <- data.frame(
        Drug = drug,
        Category = category_name,
        Direction = judgment$direction,
        Confidence = judgment$confidence,
        Spearman_r = judgment$spearman_r,
        Spearman_p_value = judgment$spearman_p,
        Pearson_r = judgment$pearson_r,
        Mean_Difference = judgment$mean_difference,
        High_Risk_Mean_Sensitivity = judgment$high_mean,
        Low_Risk_Mean_Sensitivity = judgment$low_mean,
        N_samples = nrow(plot_data),
        High_Risk_N = sum(plot_data$Risk_Group == "High"),
        Low_Risk_N = sum(plot_data$Risk_Group == "Low")
      )
      
      stats_filename <- file.path(save_dir, paste0(category_name, "_", drug, "_statistics.csv"))
      write.csv(stats_data, stats_filename, row.names = FALSE)
      
      cat("    ", drug, ":", judgment$direction, "(", judgment$confidence, "), r =", 
          round(judgment$spearman_r, 3), ", p =", format(judgment$spearman_p, digits = 3), 
          " [已保存]\n")
      
    }, error = function(e) {
      cat("    分析", drug, "时出错:", e$message, "\n")
    })
  }
  
  # 保存该类别所有药物的判断汇总
  if(nrow(all_judgments) > 0) {
    category_summary_file <- file.path(data_dir, paste0(category_name, "_all_drugs_summary.csv"))
    write.csv(all_judgments, category_summary_file, row.names = FALSE)
    cat("保存", category_name, "类别汇总:", category_summary_file, "\n")
  }
  
  return(all_judgments)
}

# ============================================================================
# 7. 判断标准文件创建函数
# ============================================================================

create_judgment_criteria_file <- function(data_dir) {
  criteria_file <- file.path(data_dir, "CTRP2_correlation_judgment_criteria.txt")
  
  sink(criteria_file)
  cat("CTRP2药物敏感性相关性判断标准\n")
  cat("===============================\n\n")
  
  cat("1. 数据来源：\n")
  cat("   - 药物敏感性数据：CTRP2数据库预测结果\n")
  cat("   - 风险评分：TCGA_LUAD队列风险评分\n")
  cat("   - 分析对象：TCGA样本的预测药物敏感性与风险评分相关性\n\n")
  
  cat("2. 主要图表：\n")
  cat("   - Spearman相关系数条形图：显示风险评分与药物敏感性相关性\n")
  cat("   - 相关性散点图：可视化风险评分-药物敏感性关系\n")
  cat("   - 散点图尺寸：6×4.5英寸\n\n")
  
  cat("3. 判断标准：\n")
  cat("   - 主要依据：Spearman相关系数\n")
  cat("   - Spearman r >= 0.2 且 p < 0.05 → 高风险组抗性增强\n")
  cat("   - Spearman r <= -0.2 且 p < 0.05 → 高风险组敏感性增强\n")
  cat("   - 0.1 <= |r| < 0.2 且 p < 0.05 → 弱相关性\n")
  cat("   - p >= 0.05 → 无显著相关性\n\n")
  
  cat("4. 置信度等级：\n")
  cat("   - High: |r| >= 0.4 且 p < 0.05\n")
  cat("   - Medium: 0.3 <= |r| < 0.4 且 p < 0.05\n")
  cat("   - Low: 0.2 <= |r| < 0.3 且 p < 0.05\n")
  cat("   - Very_Low: 0.1 <= |r| < 0.2 或均值差异显著但相关性不显著\n\n")
  
  cat("5. 文件组织结构：\n")
  cat("   - resistant_in_high_risk/: 高风险组抗性增强的药物\n")
  cat("   - sensitive_in_high_risk/: 高风险组敏感性增强的药物\n")
  cat("   - weak_correlation/: 弱相关性或无显著性的药物\n\n")
  
  cat("6. 临床意义解释：\n")
  cat("   - Resistant_in_High_Risk: 高风险患者对该药物可能产生抗性\n")
  cat("   - Sensitive_in_High_Risk: 高风险患者对该药物可能更敏感\n")
  cat("   - 数值越大表示预测的IC50越高（即敏感性越低）\n\n")
  
  sink()
  
  cat("已创建判断标准文件:", criteria_file, "\n")
}

# ============================================================================
# 8. 综合统计报告函数
# ============================================================================

generate_comprehensive_summary <- function(all_drug_judgments) {
  if(is.null(all_drug_judgments) || nrow(all_drug_judgments) == 0) {
    cat("警告: 无药物判断数据用于生成综合统计\n")
    return(NULL)
  }
  
  cat("=== 生成综合统计报告 ===\n")
  
  # 创建综合统计文件
  summary_file <- file.path(output_dir, "CTRP2_comprehensive_analysis_summary.txt")
  
  sink(summary_file)
  cat("CTRP2药物敏感性分析综合报告\n")
  cat("==================================\n")
  cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # 总体统计
  cat("1. 总体统计\n")
  cat("----------\n")
  cat("总分析药物数:", nrow(all_drug_judgments), "\n")
  cat("总分类数:", length(unique(all_drug_judgments$Category)), "\n")
  
  # 按方向分类统计
  direction_counts <- table(all_drug_judgments$Direction)
  cat("按效应方向分布:\n")
  for(direction in names(direction_counts)) {
    cat("  -", direction, ":", direction_counts[direction], "个药物\n")
  }
  
  # 按置信度分类统计
  confidence_counts <- table(all_drug_judgments$Confidence)
  cat("按置信度分布:\n")
  for(conf in names(confidence_counts)) {
    cat("  -", conf, ":", confidence_counts[conf], "个药物\n")
  }
  
  cat("\n2. 各类别详细统计\n")
  cat("----------------\n")
  
  for(category in unique(all_drug_judgments$Category)) {
    cat_data <- all_drug_judgments[all_drug_judgments$Category == category, ]
    cat("分类:", category, "\n")
    cat("  药物数:", nrow(cat_data), "\n")
    
    # 按方向统计
    cat_directions <- table(cat_data$Direction)
    for(dir in names(cat_directions)) {
      cat("    -", dir, ":", cat_directions[dir], "\n")
    }
    
    # 显示最强相关性药物
    if(nrow(cat_data) > 0) {
      top_drug <- cat_data[which.max(abs(cat_data$Spearman_r)), ]
      cat("    最强相关药物:", top_drug$Drug, 
          "(r =", round(top_drug$Spearman_r, 3), 
          ", p =", format(top_drug$Spearman_p, digits = 3), ")\n")
    }
    cat("\n")
  }
  
  cat("3. 高置信度药物推荐\n")
  cat("------------------\n")
  
  high_conf_drugs <- all_drug_judgments[all_drug_judgments$Confidence %in% c("High", "Medium"), ]
  if(nrow(high_conf_drugs) > 0) {
    high_conf_drugs <- high_conf_drugs[order(abs(high_conf_drugs$Spearman_r), decreasing = TRUE), ]
    
    cat("高/中等置信度药物 (", nrow(high_conf_drugs), "个):\n")
    for(i in 1:min(10, nrow(high_conf_drugs))) {
      drug_info <- high_conf_drugs[i, ]
      cat("  ", i, ".", drug_info$Drug, " (", drug_info$Category, ")\n")
      cat("      方向:", drug_info$Direction, ", 置信度:", drug_info$Confidence, "\n")
      cat("      相关系数:", round(drug_info$Spearman_r, 3), 
          ", p值:", format(drug_info$Spearman_p, digits = 3), "\n")
    }
  } else {
    cat("未发现高置信度药物相关性\n")
  }
  
  sink()
  
  cat("综合统计报告已保存:", summary_file, "\n")
  
  # 生成CSV格式的汇总数据
  summary_csv <- file.path(output_dir, "CTRP2_all_drugs_analysis_results.csv")
  write.csv(all_drug_judgments, summary_csv, row.names = FALSE)
  cat("所有药物分析结果CSV已保存:", summary_csv, "\n")
}

# ============================================================================
# 9. 新增：生成全部分类的总览报告函数
# ============================================================================

generate_full_class_report <- function(all_drug_judgments) {
  if(is.null(all_drug_judgments) || nrow(all_drug_judgments) == 0) {
    cat("⚠️ 没有可用的 all_drug_judgments 数据，无法生成总报告\n")
    return(NULL)
  }
  
  cat("=== 生成全部分类总览报告 ===\n")
  
  # 将结果按 Category 排序，便于查看
  all_sorted <- all_drug_judgments %>%
    arrange(Category, Direction, desc(abs(Spearman_r)))
  
  # 保存 CSV 格式总报告
  full_csv_file <- file.path(output_dir, "CTRP2_Full_Category_Report.csv")
  write.csv(all_sorted, full_csv_file, row.names = FALSE)
  cat("✅ 已生成全部分类药物敏感性结果 CSV:", full_csv_file, "\n")
  
  # 保存文字格式总览报告
  full_txt_file <- file.path(output_dir, "CTRP2_Full_Category_Report.txt")
  sink(full_txt_file)
  cat("=======================================\n")
  cat("      CTRP2 药物敏感性总览报告\n")
  cat("=======================================\n")
  cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("数据来源: CTRP2 预测数据库\n")
  cat("风险评分: TCGA_LUAD 队列\n")
  cat("总药物数:", nrow(all_sorted), "\n")
  cat("总分类数:", length(unique(all_sorted$Category)), "\n\n")
  
  # 按分类展示详细结果
  for(cat in unique(all_sorted$Category)) {
    sub_data <- all_sorted[all_sorted$Category == cat, ]
    cat("==== 分类:", cat, "====\n")
    cat("该分类药物总数:", nrow(sub_data), "\n")
    
    # 统计各方向药物数
    direction_stat <- table(sub_data$Direction)
    cat("效应方向分布:\n")
    for(dir_name in names(direction_stat)) {
      cat("  -", dir_name, ":", direction_stat[dir_name], "个\n")
    }
    cat("\n")
    
    # 展示详细药物信息表格
    cat("详细药物列表:\n")
    display_data <- sub_data[, c("Drug", "Direction", "Confidence", 
                                 "Spearman_r", "Spearman_p", 
                                 "High_Mean", "Low_Mean", "Mean_Diff")]
    
    # 格式化数值显示
    display_data$Spearman_r <- round(display_data$Spearman_r, 3)
    display_data$Spearman_p <- format(display_data$Spearman_p, scientific = TRUE, digits = 3)
    display_data$High_Mean <- round(display_data$High_Mean, 2)
    display_data$Low_Mean <- round(display_data$Low_Mean, 2)
    display_data$Mean_Diff <- round(display_data$Mean_Diff, 2)
    
    print(display_data, row.names = FALSE)
    cat("\n")
  }
  
  # 添加总结性统计
  cat("=======================================\n")
  cat("                 总结\n")
  cat("=======================================\n")
  
  # 高置信度药物汇总
  high_confidence <- all_sorted[all_sorted$Confidence %in% c("High", "Medium"), ]
  if(nrow(high_confidence) > 0) {
    cat("高/中等置信度药物 (推荐关注):\n")
    for(i in 1:nrow(high_confidence)) {
      drug_row <- high_confidence[i, ]
      cat(sprintf("  %d. %s (%s) - %s [%s置信度]\n", 
                  i, drug_row$Drug, drug_row$Category, 
                  drug_row$Direction, drug_row$Confidence))
      cat(sprintf("     相关性: r=%.3f, p=%s\n", 
                  drug_row$Spearman_r, 
                  format(drug_row$Spearman_p, scientific = TRUE, digits = 3)))
    }
  } else {
    cat("未发现高置信度的药物敏感性相关性\n")
  }
  
  cat("\n")
  
  # 效应方向总体统计
  cat("全部药物效应方向汇总:\n")
  overall_directions <- table(all_sorted$Direction)
  for(dir in names(overall_directions)) {
    percentage <- round(overall_directions[dir] / nrow(all_sorted) * 100, 1)
    cat(sprintf("  %s: %d个 (%.1f%%)\n", dir, overall_directions[dir], percentage))
  }
  
  sink()
  cat("✅ 已生成全部分类药物敏感性总览 TXT:", full_txt_file, "\n")
  
  # 生成简化的临床报告版本
  clinical_report_file <- file.path(output_dir, "CTRP2_Clinical_Summary.csv")
  clinical_data <- all_sorted %>%
    filter(Confidence %in% c("High", "Medium", "Low")) %>%
    arrange(desc(abs(Spearman_r))) %>%
    select(Category, Drug, Direction, Confidence, Spearman_r, Spearman_p, Mean_Diff) %>%
    mutate(
      Clinical_Interpretation = case_when(
        Direction == "Resistant_in_High_Risk" ~ "高风险组可能抗药",
        Direction == "Sensitive_in_High_Risk" ~ "高风险组可能敏感",
        TRUE ~ "相关性不明确"
      ),
      Priority = case_when(
        Confidence == "High" ~ "高优先级",
        Confidence == "Medium" ~ "中优先级", 
        Confidence == "Low" ~ "低优先级",
        TRUE ~ "极低优先级"
      )
    )
  
  write.csv(clinical_data, clinical_report_file, row.names = FALSE)
  cat("✅ 已生成临床简化版报告:", clinical_report_file, "\n")
}

# ============================================================================
# 10. 主分析函数
# ============================================================================

main_ctrp2_analysis <- function() {
  cat("开始CTRP2药物敏感性分析流程\n")
  cat("========================================\n")
  
  # 步骤1: 加载和准备数据
  combined_data <- load_and_prepare_ctrp2_data()
  if(is.null(combined_data)) {
    cat("数据加载失败，终止分析\n")
    return(NULL)
  }
  
  # 步骤2: 获取药物名称并创建分类
  drug_names <- colnames(combined_data)[4:(ncol(combined_data)-1)]
  cat("可用药物总数:", length(drug_names), "\n")
  
  drug_classification <- create_ctrp2_drug_categories(drug_names)
  
  # 步骤3: 创建判断标准文件
  create_judgment_criteria_file(data_dir)
  
  # 存储所有药物判断结果
  all_drug_judgments <- data.frame()
  
  # 步骤4: 对每个药物类别进行分析
  categories <- drug_classification$major_categories
  
  for(category_name in names(categories)) {
    cat("\n=== 分析药物类别:", category_name, "===\n")
    
    category_drugs <- categories[[category_name]]
    available_drugs <- intersect(category_drugs, drug_names)
    
    if(length(available_drugs) == 0) {
      cat("警告:", category_name, "类别中没有可用药物\n")
      next
    }
    
    cat("该类别可用药物:", length(available_drugs), "个\n")
    
    # 计算相关系数
    correlation_results <- calculate_drug_correlations(combined_data, available_drugs)
    
    if(nrow(correlation_results) == 0) {
      cat("警告:", category_name, "没有有效的相关性分析结果\n")
      next
    }
    
    # 生成条形图
    plot_correlation_barplot(correlation_results, category_name, output_dir)
    
    # 生成相关性散点图并获取判断结果
    category_judgments <- plot_risk_ic50_correlation_enhanced(
      combined_data, available_drugs, category_name, 
      correlation_results, correlation_dir, data_dir
    )
    
    # 合并判断结果
    if(!is.null(category_judgments) && nrow(category_judgments) > 0) {
      all_drug_judgments <- rbind(all_drug_judgments, category_judgments)
    }
  }
  
  # 步骤5: 生成综合报告
  if(nrow(all_drug_judgments) > 0) {
    generate_comprehensive_summary(all_drug_judgments)
    generate_full_class_report(all_drug_judgments)  # 新增总览报告
  } else {
    cat("\n警告: 未找到任何显著的药物相关性结果\n")
  }
  
  cat("\n========================================\n")
  cat("CTRP2药物敏感性分析完成！\n")
  cat("结果保存在:", output_dir, "\n")
  cat("========================================\n")
  
  return(all_drug_judgments)
}

# ============================================================================
# 11. 生成最终报告函数
# ============================================================================

generate_ctrp2_report <- function() {
  cat("\n=== 生成CTRP2分析最终报告 ===\n")
  
  report_file <- file.path(output_dir, "CTRP2_Analysis_Final_Report.txt")
  
  sink(report_file)
  cat("CTRP2药物敏感性预测分析最终报告\n")
  cat("=====================================\n")
  cat("分析完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("分析脚本版本: v2.0 (增强版)\n\n")
  
  cat("1. 分析概述\n")
  cat("-----------\n")
  cat("本分析基于CTRP2药物敏感性预测数据库，\n")
  cat("评估TCGA_LUAD队列中不同风险评分患者\n")
  cat("对各类抗肿瘤药物的预测敏感性差异。\n\n")
  
  cat("2. 数据来源\n")
  cat("-----------\n")
  cat("药物预测数据: CTRP2_DrugPredictions.csv\n")
  cat("风险评分数据: TCGA_risk_scores.csv\n")
  cat("分析队列: TCGA_LUAD\n\n")
  
  cat("3. 输出文件说明\n")
  cat("---------------\n")
  cat("主要输出目录:", output_dir, "\n\n")
  cat("3.1 条形图文件 (PDF格式):\n")
  cat("   - 各类别药物Spearman相关系数条形图\n")
  cat("   - 显示风险评分与药物敏感性的相关强度和方向\n\n")
  
  cat("3.2 相关性散点图 (按效应分类保存):\n")
  cat("   - resistant_in_high_risk/: 高风险组抗性增强的药物\n")
  cat("   - sensitive_in_high_risk/: 高风险组敏感性增强的药物\n")
  cat("   - weak_correlation/: 相关性较弱的药物\n\n")
  
  cat("3.3 数据文件:\n")
  cat("   - correlation_data/: 所有分析的原始数据和统计结果\n")
  cat("   - CTRP2_all_drugs_analysis_results.csv: 完整分析结果汇总\n")
  cat("   - CTRP2_comprehensive_analysis_summary.txt: 文字版综合报告\n")
  cat("   - CTRP2_Full_Category_Report.csv: 分类总览报告(CSV)\n")
  cat("   - CTRP2_Full_Category_Report.txt: 分类总览报告(TXT)\n")
  cat("   - CTRP2_Clinical_Summary.csv: 临床简化版报告\n\n")
  
  cat("4. 分析方法\n")
  cat("-----------\n")
  cat("相关性计算: Spearman等级相关\n")
  cat("显著性阈值: p < 0.05\n")
  cat("效应大小阈值: |r| >= 0.2 (中等效应)\n")
  cat("分组方法: 基于风险评分中位数分割\n\n")
  
  cat("5. 结果解释\n")
  cat("-----------\n")
  cat("Resistant_in_High_Risk: 高风险组预测IC50更高(抗性增强)\n")
  cat("Sensitive_in_High_Risk: 高风险组预测IC50更低(敏感性增强)\n")
  cat("置信度等级: High > Medium > Low > Very_Low\n\n")
  
  cat("6. 临床意义\n")
  cat("-----------\n")
  cat("本分析结果可为精准治疗提供参考，\n")
  cat("但需要进一步的实验验证和临床试验确认。\n")
  cat("建议重点关注High和Medium置信度的药物相关性。\n\n")
  
  cat("7. 技术细节\n")
  cat("-----------\n")
  cat("样本量要求: 每个分析至少15个有效数据点\n")
  cat("缺失值处理: 逐对删除法\n")
  cat("多重比较: 未进行校正(探索性分析)\n")
  cat("图表尺寸: PDF格式，300 DPI分辨率\n\n")
  
  cat("分析完成。\n")
  
  sink()
  
  cat("最终报告已保存:", report_file, "\n")
}

# ============================================================================
# 12. 执行主分析
# ============================================================================

# 运行主分析流程
final_ctrp2_results <- main_ctrp2_analysis()

# 生成最终报告
generate_ctrp2_report()

cat("\n🎉 CTRP2药物敏感性分析流程全部完成！\n")
cat("请查看输出目录中的结果文件。\n")

