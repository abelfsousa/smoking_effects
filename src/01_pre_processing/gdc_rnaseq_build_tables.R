# create gene expression tables from TCGA-GDC rna seq data


# load R packages and functions
library(tidyverse)
source("./src/functions/gdc_rnaseq_merge_files.R")


# load command line arguments
args_list <- commandArgs(trailingOnly = TRUE)
#args_list <- c(
#  "./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/gene_expression/",
#  "./data/raw/TCGA-GDC/rna_seq/TCGA-LUAD/metadata/gdc_sample_sheet.2021-06-22.tsv",
#  "htseq.counts",
#  "Primary Tumor",
#  "tcga_luad")


arg1 <- args_list[1]
arg2 <- args_list[2]
arg3 <- args_list[3]
arg4 <- args_list[4]
arg5 <- args_list[5]


# gene list
gene_list <- read_tsv("./data/raw/TCGA-GDC/rna_seq/gene_annotation/gene_ids.txt", col_names = F) %>%
  rename(gene = X1)


# path to files directory
files_dir <- arg1


# information about the files
files_info <- data.table::fread(file = arg2, check.names = T) %>%
  as_tibble() %>%
  mutate(Expr.Unit = str_extract(File.Name, "(FPKM-UQ|FPKM|htseq.counts)"))


# assemble the gene expression table
expression_data <- merge_files(files_info, files_dir, gene_list, arg3, arg4)


# transform the gene expression table into a "tidy" format
expression_data <- expression_data %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")


# average duplicated samples by gene (coming from multiple files for the same sample)
dup_samples <- expression_data %>%
  select(sample) %>%
  distinct() %>%
  mutate(sample = str_replace(sample, "\\.{1}[xy.]+$", "")) %>%
  group_by(sample) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(duplicated = if_else(n > 1, 1, 0)) %>%
  select(-n) %>%
  arrange(desc(duplicated))

sel_dup <- expression_data %>%
  mutate(sample2 = str_replace(sample, "\\.{1}[xy.]+$", "")) %>%
  filter(sample2 %in% dup_samples[dup_samples$duplicated == 1, ]$sample) %>%
  group_by(sample2) %>%
  nest() %>%
  ungroup() %>%
  mutate(selected = map_dbl(data, corr_dup_samples)) %>%
  select(-data, sample = sample2)

expression_data <- expression_data %>%
  mutate(sample = str_replace(sample, "\\.{1}[xy.]+$", ""))

expression_data2 <- expression_data %>%
  filter(sample %in% sel_dup[sel_dup$selected == 1, ]$sample) %>%
  group_by(sample, gene) %>%
  summarise(expression = round(mean(expression))) %>%
  ungroup()

expression_data <- expression_data %>%
  filter(sample %in% dup_samples[dup_samples$duplicated == 0, ]$sample) %>%
  bind_rows(expression_data2)


# average duplicated samples by gene (same sample with multiple vials)
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
dup_samples <- expression_data %>%
  select(sample) %>%
  distinct() %>%
  mutate(sample = str_replace(sample, "[A-Z]{1}$", "")) %>%
  group_by(sample) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(duplicated = if_else(n > 1, 1, 0)) %>%
  select(-n) %>%
  arrange(desc(duplicated))

sel_dup <- expression_data %>%
  mutate(sample2 = str_replace(sample, "[A-Z]{1}$", "")) %>%
  filter(sample2 %in% dup_samples[dup_samples$duplicated == 1, ]$sample) %>%
  group_by(sample2) %>%
  nest() %>%
  ungroup() %>%
  mutate(selected = map_dbl(data, corr_dup_samples)) %>%
  select(-data, sample = sample2)

expression_data <- expression_data %>%
  mutate(sample = str_replace(sample, "[A-Z]{1}$", ""))

expression_data2 <- expression_data %>%
  filter(sample %in% sel_dup[sel_dup$selected == 1, ]$sample) %>%
  group_by(sample, gene) %>%
  summarise(expression = round(mean(expression))) %>%
  ungroup()

expression_data <- expression_data %>%
  filter(sample %in% dup_samples[dup_samples$duplicated == 0, ]$sample) %>%
  bind_rows(expression_data2)


# transform the gene expression data back to a matrix format
expression_data <- expression_data %>%
  pivot_wider(names_from = "sample", values_from = "expression")


# export gene expression table
file_path <- str_c("./data/processed/TCGA-GDC/rna_seq/", arg5, "_", str_replace_all(arg4, " ", ""), "_", arg3, ".txt.gz", sep = "")

write_tsv(x = expression_data, file = file_path)

