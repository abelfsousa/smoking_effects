#!/usr/bin/env bash

# TCGA LUAD
Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"htseq.counts" \
"Primary Tumor" \
"tcga_luad"


Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"FPKM" \
"Primary Tumor" \
"tcga_luad"


Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"htseq.counts" \
"Solid Tissue Normal" \
"tcga_luad"


Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"FPKM" \
"Solid Tissue Normal" \
"tcga_luad"





# TCGA LUSC
Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"htseq.counts" \
"Primary Tumor" \
"tcga_lusc"


Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"FPKM" \
"Primary Tumor" \
"tcga_lusc"


Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"htseq.counts" \
"Solid Tissue Normal" \
"tcga_lusc"


Rscript ./src/01_pre_processing/gdc_rnaseq_build_tables.R \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/gene_expression/" \
"./data/raw/TCGA-GDC/rna_seq/TCGA-LUSC/metadata/gdc_sample_sheet.2021-06-22.tsv" \
"FPKM" \
"Solid Tissue Normal" \
"tcga_lusc"
