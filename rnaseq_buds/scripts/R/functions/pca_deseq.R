#' Perform Principal Components Analysis on a DESeqTransform object
#' 
#' This function is based on the `DESeq2::plotPCA()` function, but returns the 
#' results of `prcomp` in a tidy list format. This is more flexible for further 
#' custom plotting and exploring factor loadings of the PCA.
#' 
#' @param x an object of class DESeqTransform
#' @ntop number of most-variable genes to select. Igored if "genes" is specified.
#' @genes character vector of specific genes to use
#' 
#' @return a list with four `data.frame` objects: eigen_vectors, eigen_values, 
#' factor_loadings and the original data.
prcomp.DESeqTransform <- function(x, ntop = 500L, genes = NULL, ...){
  
  # Get sample info
  sample_info <- colData(x) %>% as.data.frame()
  
  # Get counts
  x <- SummarizedExperiment::assay(x)
  
  if(!is.null(genes)){
    message("Only using ", genes, " genes as requested.")
    
    if(!all(genes %in% rownames(x))) stop("Not all provided genes are in the gene count matrix.")
    
    selected_genes <- which(genes %in% rownames(x))
    
  } else if(is.integer(ntop) & ntop < nrow(x)){
    message("Only using ", ntop, " most variable genes.")
    
    # calculate the variance for each gene
    rv <- genefilter::rowVars(x)
    
    # select the ntop genes by variance
    selected_genes <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
  } else {
    message("Using all ", nrow(x), " genes.")
    selected_genes <- 1:nrow(x)
  }
  
  # Get the data for those genes
  selected_expr <- x[selected_genes,]
  
  # perform a PCA on the data in assay(x) for the selected genes
  ## Need to transpose the matrix as prcomp clusters by rows
  pca <- prcomp(t(selected_expr), ...)
  
  #### Eigen vectors table ####
  # Get sample information from DESeq x
  # and bind the PCA vectors
  pc_scores <- sample_info %>% 
    bind_cols(as.data.frame(pca$x))
  
  #### Eigen values table ####
  eigen_values <- data.frame(PC = colnames(pca$x), stdev = pca$sdev) %>% 
    mutate(var = stdev^2, var_pct = var/sum(var), cum_var = cumsum(var_pct),
           PC = fct_inorder(PC))
  
  #### Factor loadings table ####
  factor_loadings <- pca$rotation %>% 
    as.data.frame() %>% 
    mutate(gene = row.names(.)) %>% 
    select(gene, everything())
  
  #### Convert the original data to a data.frame ####
  selected_expr <- selected_expr %>% 
    as.data.frame() %>% 
    rename_all(funs(paste0("sample", .))) %>% 
    mutate(gene = rownames(.)) %>% 
    select(gene, everything())
  
  # Return a list with each of these xs
  return(list(pc_scores = pc_scores, 
              eigen_values = eigen_values, 
              loadings = factor_loadings,
              original = selected_expr))
  
}
