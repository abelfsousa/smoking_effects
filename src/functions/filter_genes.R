

filter_genes <- function(expression_data, gene_ids, strategy = "by_median_value", median_cutoff = 1, abs_cutoff = NULL, abs_perct = NULL){

  # filter out low-expressed genes
  # returns a vector that contains the ensembl gene IDs of the genes that passed the filter

  if( strategy == "by_absolute_value" & c(is.null(abs_cutoff) | is.null(abs_perct)))
    stop("Cutoff values for absolute filtering strategy are missing")
  else if( strategy == "by_absolute_value" & c(!is.numeric(abs_cutoff) | !is.numeric(abs_perct)))
    stop("Cutoff values for absolute filtering strategy were wrongly given")

  if (strategy == "by_median_value"){
    # select the genes whose median expression across samples is higher than "median_cutoff"
    keep <- apply(expression_data, 1, function(x) median(x, na.rm=T) >= median_cutoff)
    keep_ids <- gene_ids[keep]
  } else if (strategy == "by_absolute_value"){
    # select genes that have at least "abs_cutoff" expression in at least "abs_perct" of samples
    keep <- rowSums( expression_data >= abs_cutoff ) >= ncol(expression_data)*abs_perct
    keep_ids <- gene_ids[keep]
  } else {
    print("Filtering strategy not recognized!")
    keep_ids <- NULL
  }

  return(keep_ids)
}


# filter_genes <- function(expression_data, gene_ids, strategy = "by_median_value", median_cutoff = 1, abs_cutoff = 1, abs_perct = 0.5){
#   
#   # filter out low-expressed genes
#   # returns a vector that contains the ensembl gene IDs of the genes that passed the filter
#   
#   if (strategy == "by_median_value"){
#     # select the genes whose median expression across samples is higher than "median_cutoff"
#     keep <- apply(expression_data, 1, function(x) median(x, na.rm=T) >= median_cutoff)
#     keep_ids <- gene_ids[keep]
#   } else if (strategy == "by_absolute_value"){
#     # select genes that have at least "abs_cutoff" expression in at least "abs_perct" of samples
#     keep <- rowSums( expression_data >= abs_cutoff ) >= ncol(expression_data)*abs_perct
#     keep_ids <- gene_ids[keep]
#   } else {
#     print("Filtering strategy not recognized!")
#     keep_ids <- NULL
#   }
#   
#   return(keep_ids)
# }
# 
