# Differential gene expression analysis between smokers and nonsmokers in TCGA lung cancer cohorts


# load R packages
library(tidyverse)
library(limma)
library(edgeR)


# source R script with functions to estimate differentially expressed genes
source("./src/functions/diff_expr_methods.R")


# source R script with functions to filter out lowly expressed genes
source("./src/functions/filter_genes.R")


# load gene expression data in raw counts format
luad_counts <- read_tsv("./data/processed/TCGA-GDC/rna_seq/tcga_luad_PrimaryTumor_htseq.counts.txt.gz")
lusc_counts <- read_tsv("./data/processed/TCGA-GDC/rna_seq/tcga_lusc_PrimaryTumor_htseq.counts.txt.gz")


# load samples metadata
luad_metadata <- read_tsv("./data/processed/metadata/tcga_luad_metadata_pancancer_atlas.txt")
lusc_metadata <- read_tsv("./data/processed/metadata/tcga_lusc_metadata_pancancer_atlas.txt")


# load gene annotation
gene_annot <- read_tsv("./data/raw/TCGA-GDC/rna_seq/gene_annotation/gencode.gene.info.v22.tsv")


# select only protein coding and lincRNA genes
sel_genes <- gene_annot %>%
  filter(gene_type %in% c("protein_coding", "lincRNA"))

luad_counts <- luad_counts %>%
  semi_join(sel_genes, by = c("gene" = "gene_id"))

lusc_counts <- lusc_counts %>%
  semi_join(sel_genes, by = c("gene" = "gene_id"))


# remove outlier samples based on PCA/UMAP and hierarchical clustering analyses
luad_clusters <- read_tsv(file = "./output/files/02_exploratory_analysis/tcga_luad_outlier_samples_hclust_200_height.txt")
lusc_clusters <- read_tsv(file = "./output/files/02_exploratory_analysis/tcga_lusc_outlier_samples_hclust_190_height.txt")

luad_counts <- luad_counts %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ str_sub(.x, 1, 12)) %>%
  select_if(.predicate = colnames(.) %in% c("gene", luad_clusters[luad_clusters$cluster == 1, ]$sample))

lusc_counts <- lusc_counts %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ str_sub(.x, 1, 12)) %>%
  select_if(.predicate = colnames(.) %in% c("gene", lusc_clusters[lusc_clusters$cluster == 1, ]$sample))


# -- differential gene expression analysis between smokers and nonsmokers in LUAD and LUSC cohorts


# select covariates
luad_sel_covars <- luad_metadata %>%
  filter(sample %in% colnames(luad_counts)[-1]) %>%
  select(sample, gender, age_at_diagnosis, race, ethnicity, tissue_source_site, pathologic_stage, histological_type, smoker) %>%
  filter(!is.na(smoker) & !is.na(pathologic_stage)) %>%
  #filter(if_all(.cols = c("smoker", "pathologic_stage"), ~ !is.na(.x))) %>%
  #filter_at(.vars = c("smoker", "pathologic_stage"), all_vars(!is.na(.))) %>%
  mutate(age_at_diagnosis = replace_na(age_at_diagnosis, median(age_at_diagnosis, na.rm = T))) %>%
  mutate(race = replace_na(race, names(sort(table(race), decreasing = T)[1]))) %>%
  mutate(ethnicity = replace_na(ethnicity, names(sort(table(ethnicity), decreasing = T)[1]))) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample")

lusc_sel_covars <- lusc_metadata %>%
  filter(sample %in% colnames(lusc_counts)[-1]) %>%
  select(sample, gender, age_at_diagnosis, race, ethnicity, tissue_source_site, pathologic_stage, histological_type, smoker) %>%
  filter(!is.na(smoker) & !is.na(pathologic_stage)) %>%
  #filter(if_all(.cols = c("smoker", "pathologic_stage"), ~ !is.na(.x))) %>%
  #filter_at(.vars = c("smoker", "pathologic_stage"), all_vars(!is.na(.))) %>%
  mutate(age_at_diagnosis = replace_na(age_at_diagnosis, median(age_at_diagnosis, na.rm = T))) %>%
  mutate(race = replace_na(race, names(sort(table(race), decreasing = T)[1]))) %>%
  mutate(ethnicity = replace_na(ethnicity, names(sort(table(ethnicity), decreasing = T)[1]))) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample")


# set up a count matrix for each lung cancer cohort
luad_mat <- luad_counts %>%
  select_if(.predicate = colnames(.) %in% c("gene", rownames(luad_sel_covars))) %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()

lusc_mat <- lusc_counts %>%
  select_if(.predicate = colnames(.) %in% c("gene", rownames(lusc_sel_covars))) %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()


# select only the genes with median log2CPM >= 1 across samples
luad_mat_genes <- filter_genes(cpm(luad_mat, log = T, prior.count = 1), rownames(luad_mat), strategy = "by_median_value", median_cutoff = 1)
lusc_mat_genes <- filter_genes(cpm(lusc_mat, log = T, prior.count = 1), rownames(lusc_mat), strategy = "by_median_value", median_cutoff = 1)

luad_mat <- luad_mat[luad_mat_genes, rownames(luad_sel_covars)]
lusc_mat <- lusc_mat[lusc_mat_genes, rownames(lusc_sel_covars)]


# check if the colnames and rownames of the count and covariate matrices are in the same order
all(colnames(luad_mat) == rownames(luad_sel_covars))
all(colnames(lusc_mat) == rownames(lusc_sel_covars))


# calculate the differentially expressed genes between smokers and nonsmokers adjusted for additional covariates
luad_design_matrix <- model.matrix(~ gender + age_at_diagnosis + race + ethnicity + tissue_source_site + pathologic_stage + histological_type + smoker, data = luad_sel_covars)
lusc_design_matrix <- model.matrix(~ gender + age_at_diagnosis + race + ethnicity + tissue_source_site + pathologic_stage + histological_type + smoker, data = lusc_sel_covars)

luad_diff_genes <- limma_diff_expression(luad_design_matrix, luad_mat, gene_annot = gene_annot[, c("gene_id", "gene_name", "gene_type")])
luad_diff_genes <- luad_diff_genes$genes

lusc_diff_genes <- limma_diff_expression(lusc_design_matrix, lusc_mat, gene_annot = gene_annot[, c("gene_id", "gene_name", "gene_type")])
lusc_diff_genes <- lusc_diff_genes$genes


# export differentially expressed genes
write_tsv(luad_diff_genes, "./output/files/03_diff_expression/tcga_luad_degs_smokerVSnonsmoker.txt")
write_tsv(lusc_diff_genes, "./output/files/03_diff_expression/tcga_lusc_degs_smokerVSnonsmoker.txt")
