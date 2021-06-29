
merge_files <- function(files_info, files_dir, gene_ids, expr_unit, sample_type){
  
  # this function binds together several gene expression files per sample into a single table
  # use with GDC gene expression files
  
  # select the file types to collect
  sel_files <- files_info %>%
    filter(Expr.Unit == expr_unit & Sample.Type == sample_type)
  
  expression_file <- gene_ids
  
  for(i in 1:nrow(sel_files)){
    
    file_i <- paste0(files_dir, as.character(sel_files[i,"File.Name"]))
    file_i <- data.table::fread( file = file_i, col.names = c("gene", str_replace_all(as.character(sel_files[i,"Sample.ID"]), "-", ".")) ) %>%
      as_tibble()
    
    expression_file <- inner_join(expression_file, file_i, by = "gene")
    
  }
  
  return(expression_file)
  
}


corr_dup_samples <- function(df){
  
  mat <- df %>%
    pivot_wider(names_from = "sample", values_from = "expression") %>%
    select(-gene) %>%
    as.matrix()
  
  corr_mat <- cor(mat)
  up_triangle <- upper.tri(corr_mat)
  
  corrs <- corr_mat[up_triangle]
  res <- all(corrs > 0.4)
  
  return(res)
}

