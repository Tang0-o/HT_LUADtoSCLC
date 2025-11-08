# /home/data/tmh_project/SCLC/Fig3_multicohort_v1/Rscripts/12_mutation_analysis_functions.R

# Source the function script
source("/home/data/tmh_project/SCLC/Fig3_multicohort_v1/Rscripts/12_mutation_analysis_functions.R")

# Parameters
output_dir <- "5_TCGA_mutation"
dir.create(output_dir, showWarnings = FALSE)

clustering_file <- "/home/data/tmh_project/SCLC/Fig3_cohort_bulk_human/1_TCGA/ssgsea/clustering_output_complete.rds"
key_genes_NT <- c("TP53", "RB1", "MYC", "SMAD4", "ETV1", "NOTCH2", "PIK3CA")
key_genes_LUAD_SCLC <- c("KMT2C", "NCOR2", "DNMT3A", "KDM5B", "SOX2", "EZH2", "PTEN", "EGFR", "KRAS")
group_names <- c("Low", "Medium", "High")
group_labels <- c("transcluster_cluster_low", "transcluster_cluster_medium", "transcluster_cluster_high")

# 1. Load and Prepare Data
luad <- load_and_prepare_data(study = "LUAD", source = "MC3", clustering_file = clustering_file)

# Debug: Check transCluster values
message("Unique transCluster values:")
print(table(luad@clinical.data$transCluster))
message("Sample of clinical data:")
print(head(luad@clinical.data[, c("Tumor_Sample_Barcode", "Tumor_Sample_Barcode_min", "transCluster")]))

# 2. Subset MAF Objects
transCluster_low <- safe_subset_maf(luad, "transcluster_cluster_low")
transCluster_medium <- safe_subset_maf(luad, "transcluster_cluster_medium")
transCluster_high <- safe_subset_maf(luad, "transcluster_cluster_high")

# Create a named list
maf_objects <- list(
  Low = transCluster_low,
  Medium = transCluster_medium,
  High = transCluster_high
)

# 3. Mutation Statistics
group_stats <- do.call(rbind, lapply(seq_along(maf_objects), function(i) {
  calculate_group_mutation_stats(maf_objects[[i]], group_names[i])
}))
write.csv(group_stats, file.path(output_dir, "group_mutation_stats.csv"), row.names = FALSE)

# 4. Key Genes Analysis (Separate Plots for NT and LUAD_SCLC)
key_genes_NT_data <- do.call(rbind, lapply(seq_along(maf_objects), function(i) {
  calculate_mutation_frequency(maf_objects[[i]], key_genes_NT, group_names[i])
}))
write.csv(key_genes_NT_data, file.path(output_dir, "key_genes_NT_mutations.csv"), row.names = FALSE)
plot_key_genes(key_genes_NT_data, output_dir, filename = "key_genes_NT_mutations_plot.pdf", gene_order = key_genes_NT, width = 5, x_axis_size = 11)

key_genes_LUAD_SCLC_data <- do.call(rbind, lapply(seq_along(maf_objects), function(i) {
  calculate_mutation_frequency(maf_objects[[i]], key_genes_LUAD_SCLC, group_names[i])
}))
write.csv(key_genes_LUAD_SCLC_data, file.path(output_dir, "key_genes_LUAD_SCLC_mutations.csv"), row.names = FALSE)
plot_key_genes(key_genes_LUAD_SCLC_data, output_dir, filename = "key_genes_LUAD_SCLC_mutations_plot.pdf", gene_order = key_genes_LUAD_SCLC, x_axis_size = 11)

# 5. Statistical Testing
key_genes_NT_stats <- key_genes_statistical_test(key_genes_NT, maf_objects[["High"]], maf_objects[["Low"]])
write.csv(key_genes_NT_stats, file.path(output_dir, "key_genes_NT_statistical_test.csv"), row.names = FALSE)

key_genes_LUAD_SCLC_stats <- key_genes_statistical_test(key_genes_LUAD_SCLC, maf_objects[["High"]], maf_objects[["Low"]])
write.csv(key_genes_LUAD_SCLC_stats, file.path(output_dir, "key_genes_LUAD_SCLC_statistical_test.csv"), row.names = FALSE)

# 6. Somatic Interactions
somatic_summaries <- do.call(rbind, lapply(group_names, function(group) {
  create_somatic_interactions(maf_objects[[group]], c(key_genes_NT, key_genes_LUAD_SCLC), tolower(group), output_dir)
}))
write.csv(somatic_summaries, file.path(output_dir, "somatic_interactions_summary.csv"), row.names = FALSE)

# Extract and plot interacting genes
extract_genes_from_pairs <- function(gene_pairs) {
  genes <- character(0)
  for (pair_str in gene_pairs) {
    if (!is.na(pair_str) && pair_str != "None") {
      pairs <- unlist(strsplit(pair_str, "; "))
      for (pair in pairs) {
        genes <- c(genes, unlist(strsplit(pair, "-")))
      }
    }
  }
  unique(genes)
}

me_genes <- extract_genes_from_pairs(somatic_summaries$Gene_Pairs[somatic_summaries$Pattern == "Mutual Exclusive"])
co_genes <- extract_genes_from_pairs(somatic_summaries$Gene_Pairs[somatic_summaries$Pattern == "Co-occurring"])
significant_genes <- unique(c(me_genes, co_genes))

if (length(significant_genes) == 0) {
  significant_genes <- c(key_genes_NT, key_genes_LUAD_SCLC)
  warning("No significant gene interactions found. Using key genes.")
}

interacting_freq <- do.call(rbind, lapply(seq_along(maf_objects), function(i) {
  calculate_mutation_frequency(maf_objects[[i]], significant_genes, group_names[i])
}))
interacting_freq$InteractionType <- sapply(interacting_freq$Gene, function(gene) {
  if (gene %in% me_genes && gene %in% co_genes) "Both"
  else if (gene %in% me_genes) "Mutual Exclusive"
  else if (gene %in% co_genes) "Co-occurring"
  else "None"
})
interacting_freq$Gene <- factor(interacting_freq$Gene, levels = significant_genes)
interacting_freq$Group <- factor(interacting_freq$Group, levels = group_names)
interacting_freq$InteractionType <- factor(interacting_freq$InteractionType, levels = c("Mutual Exclusive", "Co-occurring", "Both", "None"))
write.csv(interacting_freq, file.path(output_dir, "interacting_genes_mutation_frequency.csv"), row.names = FALSE)
plot_interacting_genes(interacting_freq, output_dir)

# Write somatic interactions report
report_file <- file.path(output_dir, "somatic_interactions_report.txt")
writeLines(paste(
  "Somatic Interactions Analysis Report\n",
  "==================================\n\n",
  "Total Mutual Exclusive Genes: ", length(me_genes), "\n",
  "Total Co-occurring Genes: ", length(co_genes), "\n",
  "Genes showing both patterns: ", paste(intersect(me_genes, co_genes), collapse = ", "), "\n\n",
  paste(sapply(group_names, function(g) {
    paste(
      g, "TransCluster Group:\n",
      "Mutual Exclusive Pairs: ", somatic_summaries$Gene_Pairs[somatic_summaries$Group == tolower(g) & somatic_summaries$Pattern == "Mutual Exclusive"], "\n",
      "Co-occurring Pairs: ", somatic_summaries$Gene_Pairs[somatic_summaries$Group == tolower(g) & somatic_summaries$Pattern == "Co-occurring"], "\n",
      "Average Mutation Frequency: ", sprintf("%.2f%%", mean(interacting_freq$MutationFrequency[interacting_freq$Group == g])), "\n\n"
    )
  }), collapse = ""),
  "Note: Only significant interactions (p < 0.05) are reported.\n"
), report_file)

# 7. Top 30 Genes Analysis
top30_data <- get_top30_genes(maf_objects, group_names)
write.csv(top30_data, file.path(output_dir, "top_30_mutation_frequencies.csv"), row.names = FALSE)
plot_top30_genes(top30_data, output_dir)

# 8. Mutation Burden Analysis
mutation_burden <- calculate_group_mutation_burden(maf_objects, group_names)
write.csv(mutation_burden, file.path(output_dir, "group_mutation_burden.csv"), row.names = FALSE)
plot_mutation_burden(mutation_burden, output_dir)

# 9. Waterfall Plot
plot_waterfall(luad, c(key_genes_NT, key_genes_LUAD_SCLC), output_dir)

# 9b. Individual Group Waterfall Plots
plot_waterfall_by_group(luad, c(key_genes_NT, key_genes_LUAD_SCLC), output_dir)

# 10. Summary Report
summary_report <- paste(
  "TCGA Mutation Analysis Summary Report\n",
  "===================================\n\n",
  "1. Group Sample Sizes:\n",
  paste(sprintf("- %s: %d samples", mutation_burden$Group, mutation_burden$TotalSamples), collapse = "\n"), "\n\n",
  "2. Mutation Burden Comparison:\n",
  paste(sprintf("- %s: %d unique mutated genes (%.2f%% of genome)", 
                mutation_burden$Group, mutation_burden$UniqueMutatedGenes, 
                mutation_burden$GenesWithMutations_Percentage), collapse = "\n"), "\n\n",
  "3. Mutation Rate per 1000 Genes:\n",
  paste(sprintf("- %s: %.3f mutations per 1000 genes", 
                mutation_burden$Group, mutation_burden$MutationRate_per1000genes), collapse = "\n"), "\n\n",
  "4. Average Mutations per Sample:\n",
  paste(sprintf("- %s: %.1f ± %.1f mutations per sample", 
                mutation_burden$Group, mutation_burden$MutationsPerSample_Mean, 
                mutation_burden$MutationsPerSample_SD), collapse = "\n"), "\n\n",
  "5. Key Findings:\n",
  "- Unified mutation frequency calculations\n",
  "- Key genes (NT and LUAD_SCLC) and somatic interactions analyzed\n",
  "- Group-level mutation burden differences identified\n"
)
writeLines(summary_report, file.path(output_dir, "comprehensive_analysis_summary.txt"))

message("Analysis completed successfully! Outputs saved to:", output_dir)
