#' 
#' Functions for BreastSubtypeR package 
#' @name NC-based
#' @import ggplot2
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import ggrepel
#' @import magrittr
#' @import impute
#' @import dplyr 
#' @import circlize
NULL


#' Map Gene IDs and Supplement Missing Data
#'
#' @name MapGeneIDs
#' @description
#' This function maps gene IDs from the provided gene expression matrix to ENTREZ IDs and supplements missing data if necessary. It includes options for deduplicating probes based on the chosen method.
#'
#' @param gene_expression_matrix A matrix of gene expression data with genes in rows and samples in columns. The data should be pre-processed (e.g., log-transformed).
#' @param featuredata A table containing feature data with at least two columns:
#'   - `"probe"`: The probe identifiers (e.g., probe names).
#'   - `"ENTREZID"`: Corresponding ENTREZ gene IDs.
#' @param method A character string specifying the method to deduplicate probes:
#'   - `"IQR"`: For Affymetrix microarray data, uses interquartile range (IQR) to select a representative probe.
#'   - `"mean"`: For Agilent microarrays, selects the mean expression of probes that map to the same gene.
#'   - `"max"`: For RNAseq data, selects the probe with the highest expression.
#' @param impute Logical. If `TRUE`, missing values (NA) in the gene expression matrix are handled (e.g., imputation may occur). If `FALSE`, rows with NA values are excluded.
#' @param verbose Logical. If `TRUE`, detailed information on the process is printed during execution.
#'
#' @return A gene expression matrix with mapped ENTREZ IDs, and optionally imputed missing values. The matrix will have genes in rows and samples in columns, with ENTREZ gene IDs as row names.
#'
#' @examples
#' data("gene_expression_data")
#' data("feature_data")
#' res <- MapGeneIDs(
#'   gene_expression_matrix = gene_expression_data,
#'   featuredata = feature_data,
#'   method = "IQR",
#'   impute = TRUE,
#'   verbose = TRUE
#' )
#'
#' @export
Mapping = function(gene_expression_matrix ,featuredata = NA, method = "max", impute = TRUE, verbose = TRUE ){

  # ## test
  # gene_expression_matrix = OSLO2EMIT0.103.genematrix_noNeg[,clinic.oslo$PatientID]
  # featuredata = anno_feature

  x = gene_expression_matrix
  y = featuredata
  samplenames = colnames(x)

  ## check feature data
  if(length(y) == 0) {
    stop("Please provide feature annotation to do probe ID mapping ")
  }
  
  ## loading genes.signature
  genes.signature = BreastSubtypeR$genes.signature
 
  ## filter by ENTREZID
  y = y[y$ENTREZID %in% genes.signature$EntrezGene.ID,]
  x = x[y$probe,]
  
  
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
    cat("Genes used in NC-based methods are covered. \n")
  } else if(verbose) {
    cat("These genes are missing for NC-based methods: \n")
    cat(missing_ID_parker, sep = "\n")
  }
  
  ##AIMS
  missing_ID_AIMS = setdiff( genes.signature[ genes.signature$SSP_based == "Yes",]$EntrezGene.ID, rownames(x) )
  if( length(missing_ID_AIMS) == 0 & verbose ){ 
    cat("Genes used in SSP-based methods are covered. \n")
  } else if(verbose) {
    cat("These genes are missing for SSP-based methods: \n")
    cat(missing_ID_AIMS, sep="\n")
  }
  
  ## get matrix for NC (symbol as rows, sample as col)
  x_NC = x[ na.omit(match( BreastSubtypeR$genes.sig50$EntrezGene.ID, rownames(x) )) ,]
  rownames(x_NC) = BreastSubtypeR$genes.sig50$Symbol[match(rownames(x_NC), BreastSubtypeR$genes.sig50$EntrezGene.ID)]
  x_NC = data.frame(x_NC)
  colnames(x_NC) = samplenames
  
  ## prepare log expr
  x_NC.log = log2(x_NC+1 )
  
  ## get matrix for AIMS (entrezID as colnames)
  x_SSP = x[ na.omit(match(as.character( genes.signature[ genes.signature$SSP_based == "Yes",]$EntrezGene.ID) , rownames(x) ) ),]
  colnames(x_SSP) = samplenames
  
  result = list(x_NC = x_NC, x_NC.log = x_NC.log , x_SSP = x_SSP )
  
  return(result)
}


#' #' Function for consensus subtype
#' @noRd
get_consensus_subtype <- function(patient_row) {
  patient_row = unlist(patient_row, use.names = FALSE)
  counts <- table(patient_row)
  max_subtype <- names(counts)[which.max(counts)]
  return(max_subtype)
}

#' Function for entropy calculation
#' @noRd
get_entropy = function(patient_row) {
  freq = table(patient_row)
  prob = freq / sum(freq)
  entropy = -sum(prob * log2(prob))
  return(entropy)
}


#' Function to get the average correlation and ROR
#' @noRd 
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

 
  ## distances.Subtype for ROR 
  sum_cols_list.Subtype = mapply(function(res_ihc){
    
    #res_ihc = res_ihc_iterative[[1]]
    
    res_ihc$distances.Subtype = as.data.frame(res_ihc$distances.Subtype )
    
    ## if FALSE, make the cell as NULL
    keep = res_ihc$predictions == consensus_subtypes
    res_ihc$distances.Subtype[ !keep, ] = as.list(rep(NA, 4 ))
    
    res = mutate_at(res_ihc$distances.Subtype, vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
    
    return(res )
    
  }, res_ihc_iterative , SIMPLIFY = FALSE, USE.NAMES = FALSE )
  
  ## count_na
  count_weight_save.Subtype <- Reduce(`+`, lapply(sum_cols_list.Subtype, function(x) {
    x[!is.na(x)] <- 1
    x[is.na(x)] <- 0
    return(x)
  }))
  
  
  ## sum all for each cell
  sum_cols_save.Subtype = Reduce(`+`, lapply(sum_cols_list.Subtype, function(x) {
    
    ## change NA cell to 0 cell
    x[is.na(x)] = 0
    
    return(x)
  }))
  
  ## get the mean for each cell
  ## only when subtype is supported by consensus_subtypes for each iteration and each patient
  sum_cols_save.Subtype = sum_cols_save.Subtype / count_weight_save.Subtype
  
  
  res = list( mean_distance = mean_cols_save, mean_distance.Subtype = sum_cols_save.Subtype, testdata =  mean_cols_save.testdata) 
  

  return(res)
}

#### Visualization  ####

#' Boxplot of Correlation per Subtype
#'
#' @name Vis_boxplot
#' @description
#' This function generates a boxplot to visualize the correlation distribution between different subtypes of breast cancer, based on the provided correlation table and subtype information.
#'
#' @param out A data frame containing the columns `"PatientID"` and `"Subtype"`. The `"PatientID"` column should have unique identifiers for each patient, and the `"Subtype"` column should specify the assigned subtype for each patient.
#' @param correlations A data frame or matrix containing the correlation values computed from NC-based methods. 
#'
#' @return A `ggplot` object representing the boxplot visualization of the correlation distributions across the different subtypes.
#'
#' @examples
#' data("OSLO2EMITOobj")
#' 
#' # Prepare data: Subtype information and correlation matrix
#' out = data.frame(
#'   PatientID = res$results$parker.original$BS.all$PatientID,
#'   Subtype = res$results$parker.original$BS.all$BS
#' )
#' correlations = res$results$parker.original$outList$distances
#' 
#' # Generate the boxplot
#' p = Vis_boxplot(out, correlations)
#' plot(p)
#'
#' @export

Vis_boxpot = function(out, correlations ){
  
  # out= data.frame(PatientID = res$results$parker.original$BS.all$PatientID,
  #                 Subtype = res$results$parker.original$BS.all$BS )
  # correlations =res$results$parker.original$outList$distances
  
  df = data.frame( predictions = out$Subtype, cor = apply(correlations, 1, max))
  
  plot =  ggplot( df, aes( x = predictions, y = cor) ) +
    geom_boxplot()+
    labs(x = "", y = "Correlation")+
    theme_classic() +
    theme( 
      axis.text = element_text( size = 12)
      )
 
  return(plot)
  
}


#' Heatmap Visualization of Gene Expression by Subtype
#'
#' @name Vis_heatmap
#' @description
#' This function generates a heatmap to visualize gene expression patterns across breast cancer subtypes, based on the provided gene expression matrix and subtype information.
#'
#' @param x A gene expression matrix, where genes are rows and samples are columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and `"Subtype"`. The `"PatientID"` column should contain unique patient identifiers, and the `"Subtype"` column should specify the assigned subtype for each patient.
#'
#' @return A `ggplot` or `heatmap` object (depending on implementation) representing the heatmap of gene expression across different subtypes.
#'
#' @examples
#' data("OSLO2EMITOobj")
#' 
#' # Prepare data: Gene expression matrix and subtype information
#' x = data_input$x_NC
#' out = data.frame(
#'   PatientID = res$results$parker.original$BS.all$PatientID,
#'   Subtype = res$results$parker.original$BS.all$BS
#' )
#' 
#' # Generate the heatmap
#' p = Vis_heatmap(x, out)
#' plot(p)
#'
#' @export

Vis_heatmap = function(x, out){

  scaled_mat = t(scale(t(x[,out$PatientID])))

  ## color
  col_fun = colorRamp2(c(min(scaled_mat), 0 , max(scaled_mat)), c("green", "black", "red"))
  
  ## column annotation
  col_anno = data.frame( row.names = out$PatientID, Subtype = out$Subtype )
  anno_col = HeatmapAnnotation(df= col_anno, show_legend = FALSE,col = list(Subtype = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" ) ))
  
  heatmap = Heatmap(scaled_mat, name = "Expr",
          col = col_fun,
          ## annotation
          top_annotation = anno_col,
          
          ## clustering
          ## as original heatmap plot
          clustering_distance_rows = "spearman",
          clustering_method_rows = "average",
          show_row_dend = FALSE,
          
          cluster_column_slices = TRUE,
          column_split = col_anno$Subtype,
          clustering_distance_columns  = "spearman",
          clustering_method_columns = "average",

          ## general
          show_column_names = FALSE,
          show_heatmap_legend = TRUE,
          row_names_gp = grid::gpar(fontsize = 6))
  
  return(heatmap)
}
  



#' PCA Plot Visualization of Gene Expression by Subtype
#'
#' @name Vis_PCA
#' @description
#' This function generates a PCA plot to visualize the principal components of gene expression data, colored by the assigned subtypes. Optionally, it can display a scree plot of eigenvalues to evaluate the explained variance.
#'
#' @param x A gene expression matrix, where genes are rows and samples are columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and `"Subtype"`. The `"PatientID"` column should contain unique patient identifiers, and the `"Subtype"` column should specify the assigned subtype for each patient.
#' @param Eigen Logical. If `TRUE`, the function will display a scree plot showing the eigenvalues of the principal components.
#'
#' @return A `ggplot` object representing the PCA plot, colored by subtype. If `Eigen` is set to `TRUE`, a scree plot of the eigenvalues is also included.
#'
#' @examples
#' data("OSLO2EMITOobj")
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x = data_input$x_NC.log
#' out = data.frame(
#'   PatientID = res$results$parker.original$BS.all$PatientID,
#'   Subtype = res$results$parker.original$BS.all$BS
#' )
#'
#' # Generate the PCA plot
#' p = Vis_PCA(x = x, out = out)
#' plot(p)
#'
#' # Generate PCA plot with scree plot of eigenvalues
#' p_with_eigen = Vis_PCA(x = x, out = out, Eigen = TRUE)
#' plot(p_with_eigen)
#'
#' @export

Vis_PCA = function(x, out, Eigen = FALSE){

  
  Subtype.color = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" )

  
  pca = prcomp(t(x), center = T, scale. = T)

  if(Eigen){
    
    # Scree plot
    variance = pca$sdev^2/ sum(pca$sdev^2) *100
    scree_data = data.frame(PC =  seq_along(variance) ,Variance = variance)
    
    screeplot = ggplot(scree_data[1:10,], aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", width = 0.5, fill = "steelblue")+
      geom_line() +
      geom_point(size=2)+
      scale_x_continuous(breaks = seq(1,10, 1) )+
      geom_text( aes(x = PC, label= paste0( round( Variance,2), "%") ),  nudge_y = 1) +
      labs(x ="Principal Component", y = "Percentage of variance Explained") +
      theme_classic()+
      theme(axis.text = element_text(size =12),
            axis.title = element_text(size = 14)
      )

    return(screeplot)
  } else {
    
    ## PCA plot
    scores = as.data.frame(pca$x)
    scores$PatientID = rownames(scores)
    scores = left_join(scores, out, by ="PatientID" )
    rownames(scores) = scores$PatientID
    
    pcaplot = ggplot(data = scores, aes(x = PC1, y = PC2, color = Subtype)) +
      geom_point( size = 2) +
      geom_vline( xintercept = 0, linetype="dashed") +
      geom_hline( yintercept = 0, linetype="dashed")+
      scale_color_manual( name = "Subtype", values = Subtype.color )+
      labs(x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 2), "% variance)"),
           y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 2), "% variance)")) +
      theme_minimal() +
      theme(axis.text = element_text(size =12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text =  element_text(size = 12)
      )
    
    return(pcaplot)
  }

  
 }

#' Pie Chart Visualization of Subtype Distribution
#'
#' @name Vis_pie
#' @description
#' This function generates a pie chart to visualize the distribution of breast cancer subtypes in a cohort, based on the provided `Subtype` data.
#'
#' @param out A data frame containing two columns: `"PatientID"` and `"Subtype"`. The `"PatientID"` column should contain unique patient identifiers, and the `"Subtype"` column should specify the assigned subtype for each patient.
#'
#' @return A `ggplot` object representing a pie chart showing the proportion of each subtype in the dataset.
#'
#' @examples
#' data("OSLO2EMITOobj")
#'
#' # Prepare data: Subtype information
#' out = data.frame(
#'   PatientID = res$results$parker.original$BS.all$PatientID,
#'   Subtype = res$results$parker.original$BS.all$BS
#' )
#'
#' # Generate the pie chart
#' p = Vis_pie(out = out)
#' plot(p)
#'
#' @export

Vis_pie = function(out){
  
  # out= data.frame(PatientID = rownames(res$results$sspbc$BS.all),
  #                 Subtype = res$results$sspbc$BS.all$BS )
  
  data = data.frame( table(out$Subtype))
  data = data %>% 
    mutate(perc = round( `Freq` / sum(`Freq`) * 100, 2),
           csum = rev(cumsum(rev(perc))), 
           pos = perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), perc/2, pos)) 
    
    
  
  Subtype.color = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" )

 plot =  ggplot(data, aes(x = "", y = perc, fill = Var1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "Subtype", values = Subtype.color ) +
    geom_text_repel(data = data,
                     aes(y = pos, label = paste0(Freq, " (", perc, "%" ,")" ) ),
                     size = 4.5, nudge_x =0.6, show.legend = FALSE , segment.color = NA) +
    coord_polar("y", direction = -1) +
    theme_void() +
    theme(
          legend.title = element_text(size = 14),
          legend.text =  element_text(size = 12)
    )
  
  
  return(plot)
  
}


#' Multi-Method Subtype Heatmap Visualization
#'
#' @name Vis_Multi
#' @description
#' This function generates a heatmap visualization of breast cancer subtypes based on multiple subtyping methods. It allows users to visualize how different methods classify samples into subtypes.
#'
#' @param out A data frame or matrix that includes the subtypes predicted by different methods. The rows should represent individual samples (e.g., `PatientID`), and columns should correspond to the subtypes predicted by different methods.
#'
#' @return A heatmap visualization of multi-method subtype results.
#'
#' @examples
#' data("OSLO2MEITOobj")
#'
#' # Assuming `res$res_subtypes` contains multi-method subtype results
#' p = Vis_Multi(res$res_subtypes[,-ncol(res$res_subtypes)])  # Removing the last entropy column
#' plot(p)
#'
#' @export

Vis_Multi = function(data){

  Labels = unique(as.vector( as.matrix( data)))
  
  ## preset
  categories = data.frame(
    Category = rep( c("NC-based", "SSP-based", "consensus") , c(8,2,1)),
    row.names = c("parker.original","genefu.scale", "genefu.robust", 
                  "cIHC","cIHC.itr", "PCAPAM50", 
                  "ssBC", "ssBC.v2", 
                  "AIMS", "sspbc", "consensus")
  )
  categories$Category = factor(categories$Category, levels =  c("NC-based", "SSP-based", "consensus") )
  
  ## color
  Subtype.color = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" )
  Category.color = setNames( c("#fb9a99", "#a6cee3", "#b2df8a") , c("NC-based", "SSP-based", "consensus") )
  
  ## make row annotation
  row_anno = data.frame(Category = categories[colnames(data),], row.names = colnames(data) )
  row_anno = HeatmapAnnotation(df =row_anno, which = c("row"), col =list(Category = Category.color ),
                               annotation_legend_param = list(title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                                                              gap = unit(2, "points"),labels_gp = grid::gpar(fontsize= 12) , border = "white"))
  # ## column annotation
  # col_anno = HeatmapAnnotation(which = c("column"), Entropy = anno_barplot( res$res_subtypes[rownames(data) ,ncol(res$res_subtypes) ],bar_width = 1 ) )

  data = data[order(data[,ncol(data)]),]
  Subtype.color = Subtype.color[ names(Subtype.color) %in% Labels]
  p =  ComplexHeatmap::Heatmap(t( as.matrix(data)), name="Subtypes", col = Subtype.color,
                          row_names_gp = grid::gpar(fontsize = 12,fontface = "bold" ),
                          right_annotation = row_anno,show_column_names = FALSE,
                          heatmap_legend_param = list(title = "Intrinsic Subtype", labels = Labels[match( names(Subtype.color), Labels)],
                                                      title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                                                      gap = unit(2, "points"),labels_gp = grid::gpar(fontsize= 12) , border = "white")
                          )
  
 return(p)
}



