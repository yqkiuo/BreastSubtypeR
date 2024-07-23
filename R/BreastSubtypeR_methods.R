#' Collection of breast cancer intrinsic subtyping methods
#' 
#' @title Collection of breast cancer intrinsic subtyping methods
#' @description BreastSubtypeR is an R package that integrates the access to intrinsic subtyping methods. 
#' 
#' @name BreastSubtypeR
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
#' function for mapping gene ID 
#' @param gene_expression_matrix Gene expression matrix
#' @param featuredata Feature data provided by user. The table should contain at least three column, which are probe(probeid or transcriptID), EntrezGene.ID and symbol. 
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select "IQR" for Affy,  "mean" for Agilent and "max" for RNAseq)
#' @param impute Logic. Please specify if there are NA data adn want keep them
#' @param verbose Logic. 
#' @export

Mapping = function(gene_expression_matrix ,featuredata, method = "mean", mapping = TRUE,impute = TRUE, verbose = TRUE, ...){
  
  arguments = rlang::dots_list(
    x = gene_expression_matrix,
    y = featuredata,
    method = method,
    mapping = mapping,
    impute = impute,
    verbose = verbose, ..., .homonyms = "last"
  )
  
  call = rlang::call2(domapping, !!!arguments)
  res = eval(call)
  
  return(res)
  
}


#' BS_parker
#' @name BS_parker
#' @description
#' This calls parker-based intrinsic subtyping methods. The methods for gene calibration include "None", "medianCtr", "meanCtr", "qCtr", "Given.mdns", or chosen platform. 
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
#' @return The subtypes estimated by parker-based subtyping
#' @export

BS_parker = function(gene_expression_matrix, phenodata, calibration = "None", internal = NA, external=NA, medians = NA, hasClinical = FALSE, ...){
  
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

#' BS_IHC
#' @name BS_IHC
#' @description
#' This predicts breast cancer intrinsic subtypes with ER balanced subset for gene centering.
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with IHC column
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The intrinsic subtype, confidential level, percentages.
#' @export

BS_IHC = function(gene_expression_matrix, phenodata, hasClinical = FALSE, ...){

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    hasClinical = hasClinical,..., .homonyms = "last"
  )
  
  
  call = rlang::call2(makeCalls.ihc, !!!arguments)
  res_IHC= eval(call)
  
  
  return(res_IHC)
}


#' BS_IHC.itr
#' @name BS_IHC.itr
#' @description
#' This predicts breast cancer intrinsic subtypes with ER subset for gene centering. It supports the selection of ER subset ratio. 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with IHC column
#' @param iterative Times to do iterative ER balanced procedure with certain ratio
#' @param ratio The options are either 1:1 or 54(ER+):64(ER-) (default). The latter wass ER ratio used for UNC230 train cohort.
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as "T" column, lymphatic node status as "NODE" column. 
#' @return The intrinsic subtypes, confidential level, percentages.
#' @export

BS_IHC.itr = function(gene_expression_matrix, phenodata, iterative = 100, ratio = 54/64, hasClinical = FALSE, ...){
  
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

#' BS_PCAPAM50
#' @name BS_PCAPAM50
#' @description
#' This calls PCA-PAM50 (Raj-Kumar, PK.) to do intrinsic subtyping. 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with ER information or IHC column
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The intrinsic subtypes estimated by PCA-PAM50
#' @export

BS_PCAPAM50 = function(gene_expression_matrix, phenodata, hasClinical =FALSE){
  
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
  df.pc1pam = data.frame(PatientID=res_PC1IHC$BS.all$PatientID, PAM50=res_PC1IHC$BS.all$BS.PC1ihc,
                         IHC=res_PC1IHC$BS.all$IHC,ER = phenodata[res_PC1IHC$BS.all$PatientID,]$ER, 
                         T = phenodata[res_PC1IHC$BS.all$PatientID,]$T, NODE= phenodata[res_PC1IHC$BS.all$PatientID,]$NODE,
                         stringsAsFactors=F)
  
  #View(df.pc1pam)
  
  arguments2 = rlang::dots_list(
    mat = gene_expression_matrix,
    df.pam = df.pc1pam,
    hasClinical = hasClinical, .homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.v1PAM, !!!arguments2)
  res_PCAPAM50 = eval(call)

  
  return(res_PCAPAM50)
  
}


#' 
#' BS_ssBC
#' @name BS_ssBC
#' @description
#' This calls ssBC to do intrinsic subtyping. 
#' It is better to group cohorts by ER+/- or TN/nTN. 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, like TN,ER and HER2
#' @param s Specify "ER" or "TN", which is included in phenodata
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The subtypes estimated by ssBC
#' @export

## only ER and TN

BS_ssBC = function(gene_expression_matrix, phenodata, s , hasClinical =FALSE, ...) {

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


#' BS_AIMS
#' @name BS_AIMS
#' @description
#' This calls AIMS to do intrinsic subtyping. 
#' 
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with ER information or IHC column
#' @return The subtypes estimated by PCA_AIMS
#' @export

## need to match EntrezID here

BS_AIMS = function(gene_expression_matrix,EntrezID ,...){
  
  require(AIMS)
  
  arguments = rlang::dots_list(
    eset =as.matrix(gene_expression_matrix),
    EntrezID = EntrezID ## need to check this part
  )
  
  call = rlang::call2(AIMS::applyAIMS, !!!arguments)
  res_AIMS = eval(call)

}


#' BS_sspbc
#' @name BS_sspbc
#' @description
#' This calls sspbc to do intrinsic subtyping. 
#' 
#' @param gene_expression_matrix Gene expression matrix, gene in row sample in column
#' @return The subtypes estimated by sspbc
#' @export

BS_sspbc = function(gene_expression_matrix, ...){
  
  ## check dependencies
  ## must do
  if (!requireNamespace("sspbc", quietly = TRUE)) {
    sspbc_1.0.tar.gz = system.file("extdata", "sspbc","sspbc_1.0.tar.gz", package = "BreastSubtypeR", mustWork = TRUE)
    install.packages(sspbc_1.0.tar.gz, repos = NULL, type="source")
    if (!requireNamespace("sspbc", quietly = TRUE)) {
      stop("Package 'sspbc' should be installed manually.")
    }
  }
  
  #data(sspbc.models)
  require(sspbc)
  
  ## adding subtype - prosigna ???
  arguments = rlang::dots_list(
    gex = gene_expression_matrix,
    id = rownames(gene_expression_matrix),
    id.type = "EntrezGene",
    ssp.name= "ssp.pam50", ..., .homonyms = "last"
  )
  
  call = rlang::call2(sspbc::applySSP, !!!arguments)
  res_sspbc = eval(call)
  
}



#' BS_Check
#' @name BS_Check
#' @description
#' This calls Ensembl method to predict BC intrinsic subtype. 
#' 
#' @param data_input  list of two Gene expression matrix, output of Mapping function; "x_parker" object should be log2 transformed manually. 
#' @return The subtypes estimated by selected methods
#' @export

BS_Check = function(data_input, phenodata, POP = TRUE, methods = NA, hasClinical = FALSE,... ){
  
  
  if(length(methods) < 2 ){
    stop("Please select proper methods")
  } 
  
  
  if(POP){ ## if it is population based cohort
  results <- sapply(methods, function(method) {
    
    ## try NC-based
    if( method == "parker.median") {
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_parker, phenodata, calibration = "Internal", internal = "medianCtr", hasClinical = hasClinical))
    }
    if( method == "parker.mean"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_parker, phenodata, calibration = "Internal", internal = "meanCtr", hasClinical = hasClinical))
    }
    if( method == "parker.quantile"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_parker, phenodata, calibration = "Internal", internal = "qCtr", hasClinical = hasClinical))
    }
    
    if(method == "cIHC"){
      ## try conventional IHC
      print(paste0(method," is running!"))
      return(BS_IHC(data_input$x_parker, phenodata, hasClinical = hasClinical))
    }
    
    if(method == "IHC.itr"){
      ## try iterative IHC
      print(paste0(method," is running!"))
      return( PAM50_IHC.itr(data_input$x_parker, phenodata, hasClinical = hasClinical))
    }
    
    if(method == "PCAPAM50"){
      print(paste0(method," is running!"))
      return(BS_PCAPAM50(data_input$x_parker, phenodata, hasClinical = hasClinical))
    }
    
    if(method == "ssBC"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_parker, phenodata, s= "ER", hasClinical = hasClinical ))
    }
    
    if(method == "ssBC_JAMA"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_parker, phenodata, s= "ER_JAMA", hasClinical = hasClinical ))
    }
    
    if(method == "AIMS"){
      print(paste0(method," is running!"))
      data("genes.signature")
      ## loading library first or model
      genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])
      return(BS_AIMS(data_input$x_AIMS[genes,], genes ))
    }
    
    if(method == "sspbc"){
      print(paste0(method," is running!"))
      return(BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_AIMS) ))
    }
  
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  }
  names(results) = paste0("res_", names(results))
  
  return(results)

}


