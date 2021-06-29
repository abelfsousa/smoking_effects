
limma_diff_expression <- function(design, expression_data, gene_annot = NULL){
  
  # differential expression analysis using limma and voom
  # this function expects raw counts as expression values
  
  # adjust for covariates when defining design matrix:
  # model.matrix(~ cov1 + cov2 + cov3 + ..., data)
  
  # DGEList object using edgeR
  y <- DGEList(counts=expression_data)
  
  # scale normalization using TMM method
  y <- calcNormFactors(y)
  
  # transform the counts using voom
  y <- voom(y, design)
  
  # fit the linear model for each gene
  fit <- lmFit(y, design)
  
  # empirical bayes statistics for differential expression
  fit <- eBayes(fit)
  
  # extract gene table with statistics
  diff_expr_table <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "p", p.value = 1, lfc = 0, adjust.method="BH")
  
  if(!is.null(gene_annot)){
    diff_expr_table <- diff_expr_table %>%
      mutate(gene_id = rownames(diff_expr_table)) %>%
      as_tibble() %>%
      inner_join(gene_annot, by = "gene_id") %>%
      select_at(vars(c("gene_id", colnames(gene_annot)[-1]), everything()))
  }
  
  res <- list(genes = diff_expr_table, fit = fit)
  
  return(res)
}
