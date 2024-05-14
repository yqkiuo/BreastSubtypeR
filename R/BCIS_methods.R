#' Collection of breast cancer PAM50 subtyping methods
#' 
#' @title Collection of breast cancer PAM50 subtyping methods
#' @description PAM50subtyping is a R package that integrates the access to PAM50 subtyping methods and 
#' visualization ...
#' 
#' @name PAM50subtyping
#' 
#' @docType package
#' 
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import AIMS
#' @importFrom rlang dots_list
#' 
#' 
NULL


#' Mapping
#' @name Mapping
#' @description
#' function for mapping ID and supplementing missing data if necessary
#' @param x Gene expression matrix
#' @param y Feature data provided by user. The table should contain at least three column, which are probe(probeid or transcriptID), EntrezGene.ID and symbol. 
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select IQR for Affy and select mean for Agilent)
#' @param impute Logic. Please specify if there are NA data adn want keep them
#' @param verbose Logic. 
#' @export

Mapping = function(x ,y, method = "mean", mapping = TRUE,impute = TRUE, verbose = TRUE, ...){
  
  arguments = rlang::dots_list(
    x = x,
    y = y,
    method = method,
    mapping = mapping,
    impute = impute,
    verbose = verbose, ..., .homonyms = "last"
  )
  
  call = rlang::call2(domapping, !!!arguments)
  res = eval(call)
  
  return(res)
  
}


#' PAM50_parker
#' @name PAM50_parker
#' @description
#' This calls parker-based PAM50 subtyping methods, including None, medianCtr, meanCtr, qCtr, Given.mdns, or chosen platform to do calibration. 
#' 
#' 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column, log of normalized
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with ER information or IHC column (what input)
#' @param calibration How to do calibration, "None"(default) means no calibration for gene expression matrix. When setting calibration =None, you dont need to set internal and external parameters.  "Internal" means calibration for gene expression matrix by itself. "External" means calibration by external cohort. 
#' @param internal Specify the strategy for internal calibration, medianCtr(default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians calculated by train cohorts. When users want to use Medians prepared by user selves, this parameter should be "Given.mdns", not platform name. 
#' @param medians If you specify "external" parameter as "Given.mdns", you should input matrix/table, 50 signatures in the first column and "Given.mdns" values in the second column.
#' @param do.mapping Logical. If it is microarray data or expression matrix at transcript level, please specify TRUE (this will call dochecking() function and run mapping step )
#' @param hasClinical Logical. Please specify if you prepared clinical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The subtypes estimated by parker-based PAM50 subtyping
#' @export

PAM50_parker = function(gene_expression_matrix, phenodata, calibration = "None", internal = NA, external=NA, medians = NA, hasClinical = FALSE, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    calibration = calibration,
    internal = internal,
    external = external,
    medians = medians,
    hasClinical = hasClinical, ...,.homonyms = "last"
  )
  

  call = rlang::call2(makeCalls.parker, !!!arguments)
  res_parker = eval(call)

  return(res_parker)
  
}

#' PAM50_IHC
#' @name PAM50_IHC
#' @description
#' This predicts intrinsic breast cancer subtypes with ER balanced subset for gene centering
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with IHC column
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The intrinsic subtype, confidential level, percentages.
#' @export

PAM50_IHC = function(gene_expression_matrix, phenodata, hasClinical = FALSE, ...){

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    hasClinical = hasClinical,..., .homonyms = "last"
  )
  
  
  call = rlang::call2(makeCalls.ihc, !!!arguments)
  res_IHC= eval(call)
  
  
  return(res_IHC)
}


#' PAM50_IHC.itr
#' @name PAM50_IHC.itr
#' @description
#' This predicts intrinsic breast cancer subtypes with ER balanced subset for gene centering
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with IHC column
#' @param iterative times to do iterative ER balanced procedure with certain ratio
#' @param ratio The options are either 1:1 or 54(ER+):64(ER-) (change it as default). The latter is ER ratio used for UNC230 train cohort
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The intrinsic subtype, confidential level, percentages.
#' @export

PAM50_IHC.itr = function(gene_expression_matrix, phenodata, iterative = 100, ratio = 54/64, hasClinical = FALSE, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    iterative = iterative,
    ratio = ratio, 
    hasClinical = hasClinical, .homonyms = "last"
    )
  
  
  call = rlang::call2(makeCalls.ihc.iterative, !!!arguments)
  res_IHC.itr = eval(call)

  return(res_IHC.itr)
}



## PCA-PAM50

#' PAM50_PCA
#' @name PAM50_PCA-PAM50
#' @description
#' This calls PCA-PAM50 (Raj-Kumar, PK.) to do PAM50 subtyping. 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with ER information or IHC column
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The intrinsic subtypes estimated by PCA-PAM50
#' @export

PAM50_PCA_PAM50 = function(gene_expression_matrix, phenodata, hasClinical =FALSE){
  
  # ## test data
  # gene_expression_matrix = data_input$x_parker
  # phenodata = FL.clinic
  # hasClinical =T

  ## first step
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    hasClinical = hasClinical, .homonyms = "last"
    )
  
  call = rlang::call2(makeCalls.PC1ihc, !!!arguments)
  res_PC1IHC = eval(call)
  
  
  ## second step
  df.pc1pam = data.frame(PatientID=res_PC1IHC$Int.sbs$PatientID, PAM50=res_PC1IHC$Int.sbs$Int.SBS.PC1ihc,
                         IHC=res_PC1IHC$Int.sbs$IHC,ER = phenodata[res_PC1IHC$Int.sbs$PatientID,]$ER, 
                         T = phenodata[res_PC1IHC$Int.sbs$PatientID,]$T, NODE= phenodata[res_PC1IHC$Int.sbs$PatientID,]$NODE,
                         stringsAsFactors=F)
  
  #View(df.pc1pam)
  
  arguments2 = rlang::dots_list(
    mat = gene_expression_matrix,
    df.pam = df.pc1pam,
    hasClinical = hasClinical, .homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.v1PAM, !!!arguments2)
  res_PCA_PAM50 = eval(call)

  
  return(res_PCA_PAM50)
  
}


#' 
#' PAM50_ssBC
#' @name PAM50_ssBC
#' @description
#' This calls ssBC to do PAM50 subtyping. 
#' It is better to group cohorts by ER+/- or TN/nTN. 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, like TN,ER and HER2
#' @param s Specify "ER" or "TN", which is included in phenodata
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The subtypes estimated by ssBC
#' @export

## only ER and TN

PAM50_ssBC = function(gene_expression_matrix, phenodata, s , hasClinical =FALSE, ...) {

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    s = s,
    hasClinical = hasClinical, ..., .homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.ssBC, !!!arguments)
  res_ssBC = eval(call)
  
  return(res_ssBC)
  
}


#' PAM50_AIMS
#' @name PAM50_AIMS
#' @description
#' This calls AIMS to do PAM50 subtyping. 
#' 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with ER information or IHC column
#' @return The subtypes estimated by PCA_AIMS
#' @export

## need to match EntrezID here

PAM50_AIMS = function(gene_expression_matrix,EntrezID ,...){
  
  require(AIMS)
  
  arguments = rlang::dots_list(
    eset =as.matrix(gene_expression_matrix),
    EntrezID = EntrezID ## need to check this part
  )
  
  call = rlang::call2(AIMS::applyAIMS, !!!arguments)
  res_PAM50_AIMS = eval(call)

}


#' PAM50_sspbc
#' @name PAM50_sspbc
#' @description
#' This calls sspbc to do PAM50 subtyping. 
#' 
#' @param gene_expression_matrix Gene expression matrix, gene in row sample in column
#' @return The subtypes estimated by sspbc
#' @export

PAM50_sspbc = function(gene_expression_matrix, ...){
  
  ## check dependencies
  ## must do
  if (!requireNamespace("sspbc", quietly = TRUE)) {
    sspbc_1.0.tar.gz = system.file("extdata", "sspbc","sspbc_1.0.tar.gz", package = "PAM50subtyping", mustWork = TRUE)
    install.packages(sspbc_1.0.tar.gz, repos = NULL, type="source")
    if (!requireNamespace("sspbc", quietly = TRUE)) {
      stop("Package 'sspbc' should be installed manually.")
    }
  }
  
  #data(sspbc.models)
  require(sspbc)
  
  arguments = rlang::dots_list(
    gex = gene_expression_matrix,
    id = rownames(gene_expression_matrix),
    id.type = "EntrezGene",
    ssp.name= "ssp.pam50", ..., .homonyms = "last"
  )
  
  call = rlang::call2(sspbc::applySSP, !!!arguments)
  res_PAM50_sspbc = eval(call)
  
}



#' PAM50_ensembl
#' @name PAM50_ensembl
#' @description
#' This calls Ensembl method to predict BC intrinsic subtype. 
#' 
#' @param data_input  list of two Gene expression matrix, output of Mapping function
#' @return The subtypes estimated by selected methods
#' @export

PAM50_Ensembl = function(data_input, phenodata, POP = TRUE, methods = NA, hasClinical = FALSE,... ){
  
  
  if(length(methods) < 2 ){
    stop("Please select proper methods")
  } 
  
  
  if(POP){ ## if it is population based cohort
  results <- sapply(methods, function(method) {
    
    ## try NC-based
    if( method == "parker.median") {
      print(paste0(method," is running!"))
      return(PAM50_parker(data_input$x_parker, phenodata, calibration = "Internal", internal = "medianCtr", hasClinical = hasClinical))
    }
    if( method == "parker.mean"){
      print(paste0(method," is running!"))
      return(PAM50_parker(data_input$x_parker, phenodata, calibration = "Internal", internal = "meanCtr", hasClinical = hasClinical))
    }
    if( method == "parker.quantile"){
      print(paste0(method," is running!"))
      return(PAM50_parker(data_input$x_parker, phenodata, calibration = "Internal", internal = "qCtr", hasClinical = hasClinical))
    }
    
    if(method == "cIHC"){
      ## try conventional IHC
      print(paste0(method," is running!"))
      return(PAM50_IHC(data_input$x_parker, phenodata, hasClinical = hasClinical))
    }
    
    if(method == "IHC.itr"){
      ## try iterative IHC
      print(paste0(method," is running!"))
      return( PAM50_IHC.itr(data_input$x_parker, phenodata, hasClinical = hasClinical))
    }
    
    if(method == "PCAPAM50"){
      print(paste0(method," is running!"))
      return(PAM50_PCA_PAM50(data_input$x_parker, phenodata, hasClinical = hasClinical))
    }
    
    if(method == "ssBC"){
      print(paste0(method," is running!"))
      return(PAM50_ssBC( data_input$x_parker, phenodata, s= "ER", hasClinical = hasClinical ))
    }
    
    if(method == "ssBC_JAMA"){
      print(paste0(method," is running!"))
      return(PAM50_ssBC( data_input$x_parker, phenodata, s= "ER_JAMA", hasClinical = hasClinical ))
    }
    
    if(method == "AIMS"){
      print(paste0(method," is running!"))
      data("genes.signature")
      ## loading library first or model
      genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])
      return(PAM50_AIMS(data_input$x_AIMS[ genes,], genes ))
    }
    
    if(method == "sspbc"){
      print(paste0(method," is running!"))
      return(PAM50_sspbc( gene_expression_matrix = as.matrix(data_input$x_AIMS) ))
    }
  
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  }
  names(results) = paste0("res_", names(results))
  
  return(results)

}


