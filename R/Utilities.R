
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



