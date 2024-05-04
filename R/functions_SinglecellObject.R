#' 
#' Functions to preprocess single cell RNAseq data
#' @name SC-preprocess
#' @import Seurat
#' @importFrom presto collapse_counts
#' @import magrittr
#' @import dplyr
#' @importFrom rlang dots_list
#' @import edgeR


NULL

#' function to get normalized pseudo bulk RNAseq of single cell 
#' @param counts_mat counts matrix where columns represent cells and rows represent features
#' @param meta_data data.frame containing cell metadata with patient column
#' @param varnames subset of 'meta_data' column names, it should be "patient" here
#' @param norm specify normalization strategy
#' @returns the normalized counts followed by logarithm-base 2 transformation 
#' @noRd

SC_pseudobulk = function(counts_mat, meta_data, varnames, norm = "cpm" ){
  
  ## add harmony to correct batch?
  
  
  arguments = rlang::dots_list(
    counts_mat = counts_mat,
    meta_data = meta_data,
    varnames= "patient", ..., .homonyms = "last"
    )
  
  
  call = rlang::call2( presto::collapse_counts, !!!arguments) 
  data_collapsed = eval(call)
  
  counts = data_collapsed$counts_mat
  counts = counts[rowSums(counts) > 0,]
  ## if filter low expressed genes ?
  
  
  ## add normalization step
  
  switch (norm,
    "cpm" = {
      y = edgeR::DGEList(counts=counts)
      counts_normalized = edgeR::cpm(y,normalized.lib.sizes = TRUE, log = TRUE ) ## log2(); prior.count = 2 
    },
    "fpkm" = { 
      y =edgeR::DGEList(counts=counts)
      counts_normalized = edgeR::rpkm(y,normalized.lib.sizes = TRUE, log = TRUE )
      },
    "tpm" ={
      
    },
    "tmm" = {
      y =edgeR::DGEList(counts=counts)
      y= edgeR::calcNormFactors( y, method = "TMM")
      counts_normalized = edgeR::cpm(y,normalized.lib.sizes = TRUE, log = TRUE )
    },
    "quantile" = {
      counts_normalized = normalizeQuantiles(counts, ties = FALSE)
      counts_normalized =log2(counts_normalized + 2 )
    }
    
  )
  
  return(counts_normalized)
  
}




