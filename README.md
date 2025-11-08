# HTRS-LUAD-SCLC Transformation Analysis

A comprehensive computational framework for identifying and validating the Histologic Transformation Regulator Signature (HTRS) driving the transition from Lung Adenocarcinoma (LUAD) to Small Cell Lung Cancer (SCLC).

## Overview

This repository contains code and analysis pipelines for:
- **Single-cell RNA-seq analysis** of mouse GEMM models undergoing LUAD-SCLC transformation
- **Regulon inference** using SCENIC to identify key transcriptional drivers
- **Human cohort validation** across bulk and single-cell RNA-seq datasets
- **Risk stratification** and prognostic modeling in LUAD patients
- **Therapeutic sensitivity prediction** for transformation-prone tumors

## Key Features

### 🔬 Molecular Insights
- Identified 13-core regulon HTRS signature driving lineage plasticity
- Characterized transcriptional reprogramming from alveolar to neuroendocrine fate
- Validated across multiple species (mouse → human) and data types (scRNA-seq → bulk RNA-seq)

### 📊 Risk Stratification
- Three-tier risk classification (low/medium/high) for transformation susceptibility
- Prognostic validation in TCGA, OAK+POPLAR, and ONCOSG cohorts
- Integration of genomic, transcriptomic, and clinical features

### 💊 Therapeutic Implications
- Drug sensitivity predictions using CTRP2 database
- Identification of SCLC-directed therapeutic vulnerabilities
- Epigenetic targeting opportunities in high-risk subgroups

## Dataset Requirements

### Mouse Models
- scRNA-seq data from GEMM (GSE223958)
- Epithelial cell isolation and subtype annotation
- SCENIC regulon inference

### Human Cohorts
- **Bulk RNA-seq**: E-MTAB-10399 (paired pre/post-transformation)
- **Validation cohorts**: TCGA-LUAD, OAK+POPLAR, ONCOSG
- **scRNA-seq atlases**: HLCA, LuCA, HTAN, GSE131907

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/HTRS-LUAD-SCLC.git
cd HTRS-LUAD-SCLC

# Install required R packages
Rscript install_dependencies.R

# Or install manually
install.packages(c("Seurat", "SCENIC", "GSVA", "AUCell", "monocle", "pheatmap", "survival", "ggplot2"))
