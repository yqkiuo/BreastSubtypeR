

#' function for mapping ID and supplementing missing data if necessary
#' @param gene_expression_matrix Gene expression matrix
#' @param featuredata Feature data provided by user. The table should contain at least two column, probe and ENTREZID. 
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select "IQR" for Affy, "mean" for Agilent and max for RNAseq
#' @param impute Logic. Please specify if there are NA data adn want keep them
#' @param verbose Logic. 
#' @export
Mapping = function(gene_expression_matrix ,featuredata = NA, method = "max", impute = TRUE, verbose = TRUE ){

  x = gene_expression_matrix
  y = featuredata
  
  ## check feature data
  if(length(y) == 0) {
    stop("Please provide feature annotation to do probe ID mapping ")
  }
  
  ## loading genes.signature
  genes.signature = BreastSubtypeR$genes.signature
  
  ## first step 
  ## for empty cells. imput or not ?
  if(sum(apply(x,2,is.na))>0 & impute){
    
    library(impute)
    if(verbose){
      probeid_NA = rownames(x)[rowSums(is.na(x))]
      sample_NA = colnames(x)[colSums(is.na(x))]
      print(paste0("The imput objects: ", probeid_NA, " in ", sample_NA))
    }
    
    x = impute.knn(x)
    x = x$data
  }
  

  # 
  ## second step 
  ## if mapping (microarray or transcript )
 
  probeid = rownames(x)
  entrezid = y$ENTREZID
  names(entrezid) = y$probe
  
  ## process probeid in input data
  entrezid = entrezid[probeid]
  ##remove NA
  entrezid = entrezid[!(is.na(entrezid))]
  x = x[names(entrezid),]
  entrezid = factor( entrezid, levels =  unique(entrezid) ) 
  ## names are unique probeid and content are redundant entrezid 
  
  
  ## This is for probeID or transcriptID 
  ## split expression matrix
  split_mat <- split( as.data.frame(x), entrezid, drop = F)
  
  # function to calculate the desired statistic
  calculate_stat <- function(mat, method) {
    switch(method,
           "mean" = apply(mat, 2, mean, na.rm = TRUE),
           "median" = apply(mat, 2, median, na.rm = TRUE),
           "max" = mat[which.max( rowSums(mat) ), ],
           "stdev" = {
             temp = mat
             stdevs = apply(mat, 1, sd, na.rm = TRUE)
             vals =  temp[match(max(stdevs),stdevs),]
             return( vals)
           },
           "iqr" = {
             temp = mat
             iqrs = apply(mat, 1, function(x) { quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE) })
             vals =  temp[match(max(iqrs),iqrs),]
             return( vals)
           }
    )
  }
  
  ## keep processed x
  x = mapply( calculate_stat, split_mat, MoreArgs = list(method = method), SIMPLIFY = T, USE.NAMES = T)
  x = apply(x, 1, unlist)
  
  ##print necessary information
  ##Parker
  missing_ID_parker = setdiff( BreastSubtypeR$genes.sig50$EntrezGene.ID, rownames(x) )
  if( length(missing_ID_parker) == 0 & verbose ){ 
    print("Genes used in NC-based methods are covered")
  } else if(verbose) {
    print("These genes are missing for NC-based methods:")
    print(missing_ID_parker)
  }
  
  ##AIMS
  missing_ID_AIMS = setdiff( genes.signature[ genes.signature$SSP_based == "Yes",]$EntrezGene.ID, rownames(x) )
  if( length(missing_ID_AIMS) == 0 & verbose ){ 
    print("Genes used in SSP-based methods are covered")
  } else if(verbose) {
    print("These genes are missing for SSP-based methods:")
    print(missing_ID_AIMS)
  }
  
  ## get matrix for NC (symbol as rows, sample as col)
  x_NC = x[ match( BreastSubtypeR$genes.sig50$EntrezGene.ID, rownames(x) ) ,]
  rownames(x_NC) = BreastSubtypeR$genes.sig50$Symbol[ match( rownames(x_NC), BreastSubtypeR$genes.sig50$EntrezGene.ID )  ]
  
  ## get matrix for AIMS (entrezID as colnames)
  x_SSP = x[ match(as.character( genes.signature[ genes.signature$SSP_based == "Yes",]$EntrezGene.ID) , rownames(x)  ),]
  
  result = list(x_NC = data.frame( x_NC), x_NC.log = log2(data.frame( x_NC) +1) , x_SSP = data.frame( x_SSP) )
  
  return(result)
}


#' Function for consensus subtype
#' @noRd 
get_consensus_subtype <- function(patient_row) {
  patient_row = unlist(patient_row, use.names = FALSE)
  counts <- table(patient_row)
  max_subtype <- names(counts)[which.max(counts)]
  return(max_subtype)
}


#' Function to get the average correlation and ROR
get_average_subtype = function(res_ihc_iterative, consensus_subtypes) {

  ## correlation and ROR to be averaged
  ## if hasclini ?? need to be added later
  sum_colnames = c("Basal","Her2","LumA", "LumB", "Normal")
  
  all_patients = names(res_ihc_iterative[[1]]$predictions )
  
  sum_cols_list = mapply(function(res_ihc){
    
    #res_ihc = as.data.frame( res_ihc_iterative[[1]])
    
    res_ihc$distances = as.data.frame(res_ihc$distances )
    
    ## if FALSE, make the cell as NULL
    keep = res_ihc$predictions == consensus_subtypes
    res_ihc$distances[ !keep, ] = as.list(rep(NA, 5 ))
    
    res = mutate_at(res_ihc$distances, vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
    
    return(res )
    
  }, res_ihc_iterative, SIMPLIFY = FALSE, USE.NAMES = FALSE )
  
  
  ## count_na
  count_weight_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
    x[!is.na(x)] <- 1
    x[is.na(x)] <- 0
    return(x)
  }))
  
  ## sum all for each cell
  sum_cols_save = Reduce(`+`, lapply(sum_cols_list, function(x) {
    
    ## change NA cell to 0 cell
    x[is.na(x)] = 0
    
    return(x)
  }))
  
  ## get the mean for each cell
  ## only when subtype is supported by consensus_subtypes for each iteration and each patient
  mean_cols_save = sum_cols_save / count_weight_save
  
  ## get the mean testdata
  sum_cols_list.testdata <- mapply(function(res_ihc){
    res_ihc$testData
  }, res_ihc_iterative, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  ## sum all for each cell
  mean_cols_save.testdata = Reduce(`+`, sum_cols_list.testdata) / length(sum_cols_list.testdata)

 
  ## distances.prosigna for ROR 
  sum_cols_list.prosigna = mapply(function(res_ihc){
    
    #res_ihc = res_ihc_iterative[[1]]
    
    res_ihc$distances.prosigna = as.data.frame(res_ihc$distances.prosigna )
    
    ## if FALSE, make the cell as NULL
    keep = res_ihc$predictions == consensus_subtypes
    res_ihc$distances.prosigna[ !keep, ] = as.list(rep(NA, 4 ))
    
    res = mutate_at(res_ihc$distances.prosigna, vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
    
    return(res )
    
  }, res_ihc_iterative , SIMPLIFY = FALSE, USE.NAMES = FALSE )
  
  ## count_na
  count_weight_save.prosigna <- Reduce(`+`, lapply(sum_cols_list.prosigna, function(x) {
    x[!is.na(x)] <- 1
    x[is.na(x)] <- 0
    return(x)
  }))
  
  
  ## sum all for each cell
  sum_cols_save.prosigna = Reduce(`+`, lapply(sum_cols_list.prosigna, function(x) {
    
    ## change NA cell to 0 cell
    x[is.na(x)] = 0
    
    return(x)
  }))
  
  ## get the mean for each cell
  ## only when subtype is supported by consensus_subtypes for each iteration and each patient
  sum_cols_save.prosigna = sum_cols_save.prosigna / count_weight_save.prosigna
  
  
  res = list( mean_distance = mean_cols_save, mean_distance.prosigna = sum_cols_save.prosigna, testdata =  mean_cols_save.testdata) 
  

  return(res)
}

#### Visualization  ####

#' Functions for visualization
#' Function for boxplot of correlation per subtype
#' @param out a data frame includes "patientID" and "Subtype"
#' @param correlations  correlations table from NC-based methods
#' @export
#' 

Vis_boxpot = function(out, correlations ){
  
  df = data.frame( predictions = out$Subtype, cor = apply(correlations, 1, max))
  
  plot =  ggplot( df, aes( x = predictions, y = cor) ) +
    geom_boxplot()
  
 return(plot)
  
}


#' Function for heatmap visualizayion
#' @param x gene expression matrix, log2 transformed
#' @param out a data frame includes "patientID" and "Subtype"
#' @export
#' 

Vis_heatmap = function(x, out){

  # ## test data
  # x = data_input$x_NC
  # out= data.frame(PatientID = res$results$parker.median$BS.all$PatientID,
  #                 Subtype = res$results$parker.median$BS.all$BS )
  # 
  
  scaled_mat = t(scale(t(x)))
  
  col_fun = colorRamp2(c(min(scaled_mat), 0 , max(scaled_mat)), c("green", "black", "red"))

  ## column annotation
  col_anno = data.frame( row.names = out$PatientID, Subtype = out$Subtype )

  anno_col = HeatmapAnnotation(df= col_anno, show_legend = FALSE,col = list(Subtype = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" ) ))
  
  heatmap = Heatmap(scaled_mat, name = "Subtype",
          col = col_fun,
          ## annotation
          top_annotation = anno_col,
          
          ## clustering
          ## as original heatmap plot
          clustering_distance_rows = "pearson",
          clustering_method_rows = "complete",
          
          cluster_column_slices = TRUE,
          column_split = col_anno$Subtype,
          clustering_distance_columns  = "pearson",
          clustering_method_columns = "complete",

          ## general
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          row_names_gp = gpar(fontsize = 8))
  
  return(heatmap)
}
  



#' Function for PCAplot 
#' @param x gene expression matrix, log2 transformed
#' @param out a data table includes "patientID" and "Subtype"
#' @param screeplot Logic. Please specify if show screeplot
#' @export
#' 

## reduce dependency ???
Vis_PCA = function(x, out, Eigen = FALSE){
  
  
  Subtype.color = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" )
  
  #x = data_input$x_parker
  
  x_pca = prcomp(t(x),scale. = TRUE )
  
  screeplot = fviz_eig(x_pca, addlabels = TRUE, ylim = c(0, 50))
  
  pcaplot = fviz_pca_ind(x_pca, label="none",mean.point = FALSE, pointshape = 16 ,
               col.ind = as.factor(out$Subtype )) + 
    scale_color_manual( name = "Subtype", values = Subtype.color )
 
  if(Eigen){
    return(screeplot)
  } else {
    return(pcaplot)
  }

  
 }



