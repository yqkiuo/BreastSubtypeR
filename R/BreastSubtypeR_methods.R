#' Collection of breast cancer intrinsic subtyping methods
#' 
#' @title Collection of breast cancer intrinsic subtyping methods
#' @description BreastSubtypeR is an R package that integrates the access to intrinsic subtyping methods. 
#' 
#' @name BreastSubtypeR
#' 
#' @docType package
#' 
#' @import AIMS
#' @importFrom data.table data.table
#' @importFrom data.table set
#' @importFrom rlang call2
#' @importFrom rlang dots_list
#' @noRd 
#' 
NULL


#' Mapping
#' @name Mapping
#' @description
#' function for mapping gene ID 
#' @param gene_expression_matrix Gene expression matrix. Probe or gene symbol in row and sample in column.
#' @param featuredata Feature data provided by user. The table should contain at least three column, which are probe(probeid or transcriptID), EntrezGene.ID and symbol. 
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select "IQR" for Affy, "mean" for Agilent and "max" for RNAseq. 
#' @param impute Logic. Please specify whether to perform impute.knn on NA data
#' @param verbose Logic. 
#' @export

Mapping = function(gene_expression_matrix, featuredata, method = "mean", mapping = TRUE,impute = TRUE, verbose = TRUE, ...){
  
  arguments = rlang::dots_list(
    x = gene_expression_matrix,
    y = featuredata,
    method = method,
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

BS_parker = function(gene_expression_matrix, phenodata, calibration = "None", internal = NA, external=NA, medians = NA, Prosigna = FALSE, hasClinical = FALSE, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    calibration = calibration,
    internal = internal,
    external = external,
    medians = medians,
    Prosigna = Prosigna,
    hasClinical = hasClinical, ...,.homonyms = "last"
  )
  

  call = rlang::call2(makeCalls.parker, !!!arguments)
  res_parker = eval(call)

  return(res_parker)
  
}

#' BS_cIHC
#' @name BS_cIHC
#' @description
#' This predicts breast cancer intrinsic subtypes with ER balanced subset for gene centering.
#' @param gene_expression_matrix Gene expression matrix, gene in row and sample in column
#' @param featuredata Annotate genes
#' @param phenodata Clinical information table, with IHC column
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The intrinsic subtype, confidential level, percentages.
#' @export

BS_cIHC = function(gene_expression_matrix, phenodata,Prosigna = FALSE ,  hasClinical = FALSE, ...){

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    Prosigna = Prosigna, 
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

BS_cIHC.itr = function(gene_expression_matrix, phenodata, iteration = 100, ratio = 54/64, Prosigna = FALSE, hasClinical = FALSE,seed=118, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    iteration = iteration,
    ratio = ratio, 
    Prosigna = Prosigna, 
    hasClinical = hasClinical,
    seed = seed, .homonyms = "last"
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

BS_PCAPAM50 = function(gene_expression_matrix, phenodata, Prosigna = FALSE, hasClinical =FALSE,seed=118){
  
  # ## test data
  # gene_expression_matrix = data_input$x_NC
  # phenodata = FL.clinic
  # hasClinical =T

  ## first step
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    Prosigna = Prosigna,
    hasClinical = hasClinical, .homonyms = "last"
    )
  
  call = rlang::call2(makeCalls.PC1ihc, !!!arguments)
  res_PC1IHC = eval(call)
  
  
  ## second step
  df.pc1pam = data.frame(PatientID=res_PC1IHC$BS.all$PatientID, PAM50=res_PC1IHC$BS.all$BS,
                         IHC=res_PC1IHC$BS.all$IHC,ER = phenodata[res_PC1IHC$BS.all$PatientID,]$ER, 
                         T = phenodata[res_PC1IHC$BS.all$PatientID,]$T, NODE= phenodata[res_PC1IHC$BS.all$PatientID,]$NODE,
                         stringsAsFactors=F)
  
  #View(df.pc1pam)
  
  arguments2 = rlang::dots_list(
    mat = gene_expression_matrix,
    df.pam = df.pc1pam,
    Prosigna = Prosigna,
    hasClinical = hasClinical,
    seed=118, .homonyms = "last"
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
#' @param s Specify "ER" or "TN" or "ER_JAMA" or "HER2+" or "TNBC". The original quantile is "ER" and "TN" of TNBC-BreastCancerRes2015.  If you choose "ER_JAMA" or "HER2+" or "TNBC", it means you choose quantile from TNBC-JAMAOncol2024. 
#' @param hasClinical Logical. Please specify if you prepared clincical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @return The subtypes estimated by ssBC
#' @export

## only ER and TN

BS_ssBC = function(gene_expression_matrix, phenodata, s , Prosigna = FALSE, hasClinical =FALSE, ...) {

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    s = s,
    Prosigna = Prosigna,
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
#' @param EntrezID A list of Entrez gene ids
#' @return The subtypes estimated by PCA_AIMS
#' @export

## need to match EntrezID here

BS_AIMS = function(gene_expression_matrix, EntrezID, ...){
  
  require(AIMS)
  # 
  # data("genes.signature")
  # ## loading library first or model
  # genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])
  # return(BS_AIMS(data_input$x_SSP[genes,], genes ))
  # 
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

BS_sspbc = function(gene_expression_matrix, ssp.name= "ssp.pam50" ,...){
  
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
    ssp.name= ssp.name, ..., .homonyms = "last"
  )
  
  call = rlang::call2(sspbc::applySSP, !!!arguments)
  res_sspbc = eval(call)

}



#' BS_Check
#' @name BS_Check
#' @description
#' This calls consensus subtype . 
#' 
#' @param data_input  list of two Gene expression matrix, output of Mapping function; "x_NC" object should be log2 transformed manually. 
#' @return The subtypes estimated by selected methods
#' @export

## return Subtype table; and every object
## return entropy level


BS_Check = function(data_input, phenodata, methods = NA, Prosigna = FALSE, hasClinical = FALSE,... ){
  # 
  # data_input = data_input
  # phenodata = clinic.oslo
  # methods = methods
  # POP = TRUE
  # Prosigna = TRUE
  # hasClinical = TRUE
  
  if(length(methods) < 2 ){
    stop("Please select two methods at least")
  } 
  
  if (length(methods[str_detect(methods, pattern =  "parker.median|parker.mean|parker.quantile|ssBC|ssBC_JAMA|cIHC|cIHC.itr|AIMS|sspbc")] ) < length(methods)){
    stop("Please provide right method names")
  }

  ## check ER and if NC-based methods are feasible
  if(  !("ER" %in% colnames(phenodata)) & ( length(methods[str_detect(methods, pattern =  "ssBC|ssBC_JAMA|cIHC|cIHC.itr")] ) > 0 ) ) {
    stop("Please do not select any of ssBC, ssBC_JAMA, cIHC and cIHC.itr")
  } 
  
  if(  !("HER2" %in% colnames(phenodata)) & ( "ssBC_JAMA" %in% methods  ) ) {
    stop("ssBC_JAMA is not supported")
  }

  ## run each method
  results = sapply(methods, function(method) {
    
    ## try NC-based
    if( method == "parker.median") {
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "medianCtr", Prosigna = Prosigna ,  hasClinical = hasClinical))
    }
    if( method == "parker.mean"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "meanCtr", Prosigna = Prosigna , hasClinical = hasClinical))
    }
    if( method == "parker.quantile"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "qCtr", Prosigna = Prosigna , hasClinical = hasClinical))
    }
    if(method == "cIHC"){
      ## try conventional IHC
      print(paste0(method," is running!"))
      return(BS_cIHC(data_input$x_NC.log, phenodata,Prosigna = Prosigna ,  hasClinical = hasClinical))
    }
    
    if(method == "cIHC.itr"){
      ## try iterative IHC
      print(paste0(method," is running!"))
      return( BS_cIHC.itr(data_input$x_NC.log, phenodata,Prosigna = Prosigna , hasClinical = hasClinical))
    }
    
    if(method == "PCAPAM50"){
      print(paste0(method," is running!"))
      return(BS_PCAPAM50(data_input$x_NC.log, phenodata,Prosigna = Prosigna , hasClinical = hasClinical))
    }
    
    if(method == "ssBC"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_NC.log, phenodata, s= "ER",Prosigna = Prosigna , hasClinical = hasClinical ))
    }
    
    if(method == "ssBC_JAMA"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_NC.log, phenodata, s= "ER_JAMA", Prosigna = Prosigna ,hasClinical = hasClinical ))
    }
    
    if(method == "AIMS"){
      print(paste0(method," is running!"))
      data("genes.signature")
      ## loading library first or model
      genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])
      res_AIMS = BS_AIMS(data_input$x_SSP[genes,], genes )
      res_AIMS$BS.all = as.data.frame( res_AIMS$cl)
      colnames(res_AIMS$BS.all) = "BS"
      
      if( Prosigna){
        res_AIMS$BS.all$BS.prosigna = res_AIMS$BS.all$BS
      }
      return(res_AIMS)
    }
    
    
    if(method == "sspbc" ){
      print(paste0(method," is running!"))
      
      res_sspbc = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.pam50"  )
      
      BS.all = as.data.frame(res_sspbc)
      colnames(BS.all) = "BS"
      
      if(Prosigna) {
        res_sspbc.prosigna = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.subtype"  )
        BS.all$BS.prosigna = res_sspbc.prosigna[,1]
      } 
      
      return(list(BS.all = BS.all ) )
    }
  
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  names(results) = paste0(names(results))
  

  res_subtypes = data.table( row_id = colnames(data_input$x_NC))
  if (Prosigna) {
    res_subtypes.prosigna = data.table(row_id = colnames(data_input$x_NC))
  }

  for (method in methods) {
    if (!is.null(results[[method]])) {
      set(res_subtypes, j = method, value = results[[method]]$BS.all$BS)
      if (Prosigna) {
        set(res_subtypes.prosigna, j = method, value = results[[method]]$BS.all$BS.prosigna)
      }
    }
  }

  ## adding consensus
  res_subtypes = as.data.frame(res_subtypes)
  rownames(res_subtypes) = colnames(data_input$x_NC); res_subtypes[,1]= NULL

  consensus.subtype = apply(res_subtypes, 1, get_consensus_subtype)
  res_subtypes$consensus.subtype = consensus.subtype


  if (Prosigna) {
    res_subtypes.prosigna = as.data.frame(res_subtypes.prosigna)
    rownames(res_subtypes.prosigna) = colnames(data_input$x_NC); res_subtypes.prosigna[,1]= NULL

    consensus.subtype_prosigna = apply(res_subtypes.prosigna, 1, get_consensus_subtype)
    res_subtypes.prosigna$consensus.subtype = consensus.subtype_prosigna
  }
  
  if (Prosigna) {
    res = list(res_subtypes = res_subtypes, res_subtypes.prosigna = res_subtypes.prosigna, results = results)
  } else {
    res = list(res_subtypes = res_subtypes, results = results)
  }

}



