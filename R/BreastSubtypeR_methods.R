#' Collection of breast cancer intrinsic subtyping methods
#' 
#' @title Collection of breast cancer intrinsic subtyping methods
#' @description BreastSubtypeR is an R package that grant the access to intrinsic subtyping methods. 
#' 
#' @name BreastSubtypeR
#' 
#' @docType package
#' 
#' @import AIMS
#' @import stringr
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
#' @param gene_expression_matrix Gene expression matrix. Probe in row and sample in column.
#' @param featuredata Feature data provided by user. The table should contain at least three column, which are probe (probeid or transcriptID), SYMBOL and ENTREZID.  
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select "IQR" for Affy, "mean" for Agilent and "max" for RNAseq. 
#' @param impute Logic. Please specify whether to perform impute.knn on NA data
#' @param verbose Logic. 
#' @return lists of three gene expression matrix. "x_NC.log" is prepared for NC-based methods and "x_SSP" for SSP-based methods. 
#' @details
#' The row can be gene symbol names in gene expression matrix/table, but you need to add one extra SYMBOL column and rename it as probe in feature table. 
#' 
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' data = OSLO2EMIT0.103.genematrix_noNeg[,clinic.oslo$PatientID]
#' data_input = Mapping(gene_expression_matrix = data, featuredata = anno_feature, impute = TRUE, verbose = TRUE )
#' 
#' @export

Mapping = function(gene_expression_matrix, featuredata, method = "max", impute = TRUE, verbose = TRUE, ...){
  
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
#' This calls parker-based intrinsic subtyping methods.  
#' 
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param calibration The calibration method to use. Options are "None", "Internal", or "External". If "Internal" is selected, see the "internal" parameter for further details. If "External" is selected, see the "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are median-centered ("medianCtr", default), mean-centered ("meanCtr"), or quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for external medians, which are calculated by the training cohort. If you want to use user-provided medians, set this parameter to "Given.mdns" and provide the medians via the "medians" parameter. 
#' @param medians If "Given.mdns" is specified for the "external" parameter, input a matrix/table where the first column contains 50 genes and the second column contains the corresponding "Given.mdns" values.
#' @param Subtype Logic. Specify whether to predict four subtypes by removing the Normal-like subtype. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @return The intrinsic subtypes estimated using the Parker-based method. 
#' 
#' @references 
#' Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. Journal of Clinical Oncology. 2009;27(8). https://doi.org/10.1200/JCO.2008.18.1370
#' Gendoo DMA, Ratanasirigulchai N, Schr"der MS, Par' L, Parker JS, Prat A, et al. Genefu: An R/Bioconductor package for computation of gene expression-based signatures in breast cancer. Bioinformatics. 2016;32(7). https://doi.org/10.1093/bioinformatics/btv693
#' 
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' res = BS_parker(data_input$x_NC.log, phenodata = NA, calibration = "Internal", internal = "medianCtr",  Subtype = FALSE, hasClinical =FALSE)
#' 
#' @export

BS_parker = function(gene_expression_matrix, phenodata = NA, calibration = "None", internal = NA, external=NA, medians = NA, Subtype = FALSE, hasClinical = FALSE, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    calibration = calibration,
    internal = internal,
    external = external,
    medians = medians,
    Subtype = Subtype,
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
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param Subtype Logic. Specify whether to predict four subtypes by removing the Normal-like subtype. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @return The intrinsic subtypes estimated using the conventinal IHC (cIHC) method.
#' 
#' @references 
#' Ciriello G, Gatza ML, Beck AH, Wilkerson MD, Rhie SK, Pastore A, et al. Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer. Cell. 2015;163(2). https://doi.org/10.1016/j.cell.2015.09.033
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' res = BS_cIHC(data_input$x_NC.log, phenodata= clinic.oslo,  Subtype = FALSE, hasClinical =FALSE)
#' 
#' @export

BS_cIHC = function(gene_expression_matrix, phenodata, Subtype = FALSE ,  hasClinical = FALSE, ...){

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    Subtype = Subtype, 
    hasClinical = hasClinical,..., .homonyms = "last"
  )
  
  
  call = rlang::call2(makeCalls.ihc, !!!arguments)
  res_IHC= eval(call)
  
  
  return(res_IHC)
}


#' BS_cIHC.itr
#' @name BS_cIHC.itr
#' @description
#' This predicts breast cancer intrinsic subtypes with ER subset for gene centering. It supports the selection of ER subset ratio. 
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param iterative Times to do iterative ER balanced procedure with certain ratio. 
#' @param ratio The options are either 1:1 or 54 (ER+) : 64 (ER-) (default). The latter was ER ratio used for UNC230 train cohort.
#' @param Subtype Logic. Specify whether to predict four subtypes by removing the Normal-like subtype. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' 
#' @return The intrinsic subtypes, confidential level, percentages.
#' 
#' @references 
#' Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature. 2012;486(7403). https://doi.org/10.1038/nature10983
#' 
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' res =  BS_cIHC.itr(data_input$x_NC.log, phenodata = clinic.oslo,  Subtype = FALSE, hasClinical =FALSE)
#' 
#' @export

BS_cIHC.itr = function(gene_expression_matrix, phenodata, iteration = 100, ratio = 54/64, Subtype = FALSE, hasClinical = FALSE,seed=118, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    iteration = iteration,
    ratio = ratio, 
    Subtype = Subtype, 
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
#' This calls PCA-PAM50 to do intrinsic subtyping. 
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param Subtype Logic. Specify whether to predict four subtypes by removing the Normal-like subtype. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' 
#' @return The intrinsic subtypes estimated by PCA-PAM50 method. 
#' 
#' @references 
#' Raj-Kumar PK, Liu J, Hooke JA, Kovatich AJ, Kvecher L, Shriver CD, et al. PCA-PAM50 improves consistency between breast cancer intrinsic and clinical subtyping reclassifying a subset of luminal A tumors as luminal B. Sci Rep. 2019;9(1). https://doi.org/10.1038/s41598-019-44339-4
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' res = BS_PCAPAM50(data_input$x_NC.log, phenodata = clinic.oslo, Subtype = FALSE, hasClinical =FALSE, seed=118)
#' 
#' @export

BS_PCAPAM50 = function(gene_expression_matrix, phenodata, Subtype = FALSE, hasClinical =FALSE,seed=118){

  samples = phenodata$PatientID

  if ( "ER"  %in% colnames(phenodata) ) {
    phenodata$ER_status = NA
    
    phenodata$ER_status[which(phenodata$ER %in% c("ER+"))] = "pos"
    phenodata$ER_status[which(phenodata$ER %in% c("ER-"))] = "neg"
    phenodata = phenodata[order(phenodata$ER_status,decreasing=T),]
  } else {
    stop("Please prepare ER status in clinical table")
  }
  
  gene_expression_matrix = gene_expression_matrix[,phenodata$PatientID]

  ## first step
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    Subtype = Subtype,
    hasClinical = hasClinical, .homonyms = "last"
    )
  
  call = rlang::call2(makeCalls.PC1ihc, !!!arguments)
  res_PC1IHC = eval(call)
  
  
  ## second step
  
  if (hasClinical) {
    df.pc1pam = data.frame(PatientID=res_PC1IHC$BS.all$PatientID, PAM50=res_PC1IHC$BS.all$BS,
                           T = phenodata[res_PC1IHC$BS.all$PatientID,]$T, NODE= phenodata[res_PC1IHC$BS.all$PatientID,]$NODE,
                           stringsAsFactors=F)
  } else {
    df.pc1pam = data.frame(PatientID=res_PC1IHC$BS.all$PatientID, PAM50=res_PC1IHC$BS.all$BS,
                           stringsAsFactors=F)
  }
  

  arguments2 = rlang::dots_list(
    mat = gene_expression_matrix,
    df.pam = df.pc1pam,
    Subtype = Subtype,
    hasClinical = hasClinical,
    seed=118, .homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.v1PAM, !!!arguments2)
  res_PCAPAM50 = eval(call)
  
  ## reorder
  res_PCAPAM50$BS.all = res_PCAPAM50$BS.all[match(samples,res_PCAPAM50$BS.all$PatientID ),]

  return(res_PCAPAM50)
  
}


#' 
#' BS_ssBC
#' @name BS_ssBC
#' @description
#' This calls ssBC to do intrinsic subtyping. 
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param s Options are "ER" or "TN" or "ER.v2" or "HER2+" or "TNBC". Specify the medians you want. The original quantile is "ER" and "TN" of TNBC-BreastCancerRes2015.  If you choose "ER.v2" or "HER2+" or "TNBC", it means you choose quantile from TNBC-JAMAOncol2024. 
#' @param Subtype Logic. Specify whether to predict four subtypes by removing the Normal-like subtype. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' 
#' @return The intrinsic subtypes estimated by ssBC method
#' 
#' @references 
#' Zhao X, Rodland EA, Tibshirani R, Plevritis S. Molecular subtyping for clinically defined breast cancer subgroups. Breast Cancer Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
#' Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L, et al. Survival, pathologic response, and genomics in CALGB 40601 (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or without lapatinib in HER2-positive breast cancer. In: Journal of Clinical Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#' 
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' res = BS_ssBC(data_input$x_NC.log, phenodata = clinic.oslo, s = "ER.v2", Subtype = FALSE, hasClinical =FALSE)
#' 
#' @export


BS_ssBC = function(gene_expression_matrix, phenodata, s , Subtype = FALSE, hasClinical =FALSE, ...) {

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    s = s,
    Subtype = Subtype,
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
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param EntrezID A list of Entrez gene ids
#' 
#' @return The subtypes estimated by AIMS method
#' 
#' @references 
#' Paquet ER, Hallett MT. Absolute assignment of breast cancer intrinsic molecular subtype. J Natl Cancer Inst. 2015;107(1). https://doi.org/10.1093/jnci/dju357
#' 
#' @examples
#' 
#' data("BreastSubtypeR")
#' data("OSLO2MEITOobj")
#' 
#' genes = as.character( BreastSubtypeR$genes.signature$EntrezGene.ID[which( BreastSubtypeR$genes.signature$AIMS == "Yes" )])
#' res = BS_AIMS(data_input$x_SSP[genes,], rownames(data_input$x_SSP[genes,]) )
#' 
#' @export

## need to match EntrezID here

BS_AIMS = function(gene_expression_matrix, EntrezID, ...){
  
  ## check dependencies
  suppressMessages(suppressWarnings({
    require(AIMS,quietly = TRUE)
  }))
  
  arguments = rlang::dots_list(
    eset =as.matrix(gene_expression_matrix),
    EntrezID = EntrezID 
  )
  
  call = rlang::call2(AIMS::applyAIMS, !!!arguments)
  
  res_AIMS = eval(call)

}


#' BS_sspbc
#' @name BS_sspbc
#' @description
#' This calls sspbc to do intrinsic subtyping. 
#' 
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param ssp.name Specify model names. Option is either "ssp.pam50" or "ssp.subtype". The latter one was used to predict four subtypes by removing the Normal-like subtype. 
#' 
#' @return The subtypes estimated by sspbc method
#' 
#' @references 
#' Staaf J, H"kkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al. RNA sequencing-based single sample predictors of molecular subtype and risk of recurrence for clinical assessment of early-stage breast cancer. NPJ Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#' 
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' res = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.pam50" )
#' 
#' @export

BS_sspbc = function(gene_expression_matrix, ssp.name= "ssp.pam50" ,...){
  
  ## check dependencies
  if (!requireNamespace("sspbc", quietly = TRUE)) {
    sspbc_1.0.tar.gz = system.file("extdata", "sspbc","sspbc_1.0.tar.gz", package = "BreastSubtypeR", mustWork = TRUE)
    install.packages(sspbc_1.0.tar.gz, repos = NULL, type="source")
    if (!requireNamespace("sspbc", quietly = TRUE)) {
      stop("Package 'sspbc' should be installed manually.")
    }
  }
  
  #data(sspbc.models)
  suppressMessages(suppressWarnings({
    require(sspbc,quietly = TRUE)
  }))
  
  

  arguments = rlang::dots_list(
    gex = gene_expression_matrix,
    id = rownames(gene_expression_matrix),
    id.type = "EntrezGene",
    ssp.name= ssp.name, ..., .homonyms = "last"
  )
  
  call = rlang::call2(sspbc::applySSP, !!!arguments)
  res_sspbc = eval(call)

}



#' BS_Multi
#' @name BS_Multi
#' @description
#' This calls consensus subtyping and users can specify subtyping methods. 
#' 
#' @param data_input The output of Mapping() function. 
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param methods Specify methods. 
#' @param Subtype Logic. Specify whether to predict four subtypes by removing the Normal-like subtype. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @return The subtypes estimated by selected methods
#' 
#' @examples
#' 
#' data("OSLO2MEITOobj")
#' methods = c( "parker.original", "genefu.scale", "genefu.robust")
#' res.test = BS_Multi(data_input = data_input, phenodata = clinic.oslo, methods = methods, Subtype = FALSE, hasClinical = FALSE)
#' 
#' 
#' @export

BS_Multi = function(data_input, phenodata, methods = NA, Subtype = FALSE, hasClinical = FALSE,... ){

  if(length(methods) < 2 ){
    stop("Please select two methods at least")
  } 
  
  if (length(methods[str_detect(methods, pattern =  "parker.original|genefu.scale|genefu.robust|ssBC|ssBC.v2|cIHC|cIHC.itr|PCAPAM50|AIMS|sspbc")] ) < length(methods)){
    stop("Please provide right method names")
  }

  ## check ER and if methods are feasible
  if(  !("ER" %in% colnames(phenodata)) & ( length(methods[str_detect(methods, pattern =  "ssBC|ssBC.v2|cIHC|cIHC.itr")] ) > 0 ) ) {
    stop("Please do not select any of ssBC, ssBC.v2, cIHC and cIHC.itr")
  } 
  
  if(  !("HER2" %in% colnames(phenodata)) & ( "ssBC.v2" %in% methods  ) ) {
    stop("ssBC.v2 is not supported")
  }

  ## run each method
  results = sapply(methods, function(method) {
    
    ## try NC-based
    if( method == "parker.original") {
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "medianCtr", Subtype = Subtype ,  hasClinical = hasClinical))
    }
    if( method == "genefu.scale"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "meanCtr", Subtype = Subtype , hasClinical = hasClinical))
    }
    if( method == "genefu.robust"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "qCtr", Subtype = Subtype , hasClinical = hasClinical))
    }
    if(method == "cIHC"){
      ## try conventional IHC
      print(paste0(method," is running!"))
      return(BS_cIHC(data_input$x_NC.log, phenodata,Subtype = Subtype ,  hasClinical = hasClinical))
    }
    
    if(method == "cIHC.itr"){
      ## try iterative IHC
      print(paste0(method," is running!"))
      return( BS_cIHC.itr(data_input$x_NC.log, phenodata,Subtype = Subtype , hasClinical = hasClinical))
    }
    
    if(method == "PCAPAM50"){
      print(paste0(method," is running!"))
      return(BS_PCAPAM50(data_input$x_NC.log, phenodata, Subtype = Subtype , hasClinical = hasClinical))
    }
    
    if(method == "ssBC"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_NC.log, phenodata, s= "ER",Subtype = Subtype, hasClinical = hasClinical ))
    }
    
    if(method == "ssBC.v2"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_NC.log, phenodata, s= "ER.v2", Subtype = Subtype, hasClinical = hasClinical ))
    }
    
    if(method == "AIMS"){
      print(paste0(method," is running!"))
      data("BreastSubtypeR")
      genes.signature = BreastSubtypeR$genes.signature
      
      ## loading library first or model
      genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])
      res_AIMS = BS_AIMS(data_input$x_SSP[genes,], rownames(data_input$x_SSP[genes,]) )
      res_AIMS$BS.all = data.frame( PatientID = rownames(res_AIMS$cl) ,
                                    BS = res_AIMS$cl[,1])
      
      if( Subtype){
        res_AIMS$BS.all$BS.Subtype = res_AIMS$BS.all$BS
      }
      return(res_AIMS)
    }
    
    
    if(method == "sspbc" ){
      print(paste0(method," is running!"))
      
      res_sspbc = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.pam50"  )
      
      BS.all = data.frame( PatinetID = rownames(res_sspbc),
                           BS = res_sspbc[,1],
                           row.names = rownames(res_sspbc) )

      if(Subtype) {
        res_sspbc.Subtype = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.subtype"  )
        BS.all$BS.Subtype = res_sspbc.Subtype[,1]
      } 
      
      return(list(BS.all = BS.all ) )
    }
  
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  names(results) = paste0(names(results))
  

  res_subtypes = data.table( row_id = colnames(data_input$x_NC))
  if (Subtype) {
    res_subtypes.Subtype = data.table(row_id = colnames(data_input$x_NC))
  }

  for (method in methods) {
    if (!is.null(results[[method]])) {
      set(res_subtypes, j = method, value = results[[method]]$BS.all$BS)
      if (Subtype) {
        set(res_subtypes.Subtype, j = method, value = results[[method]]$BS.all$BS.Subtype)
      }
    }
  }

  ## adding consensus
  res_subtypes = as.data.frame(res_subtypes)
  rownames(res_subtypes) = colnames(data_input$x_NC); res_subtypes[,1]= NULL

  consensus = apply(res_subtypes, 1, get_consensus_subtype)
  res_subtypes$consensus = consensus


  if (Subtype) {
    res_subtypes.Subtype = as.data.frame(res_subtypes.Subtype)
    rownames(res_subtypes.Subtype) = colnames(data_input$x_NC); res_subtypes.Subtype[,1]= NULL

    consensus_Subtype = apply(res_subtypes.Subtype, 1, get_consensus_subtype)
    res_subtypes.Subtype$consensus = consensus_Subtype
  }
  
  if (Subtype) {
    res = list(res_subtypes = res_subtypes, res_subtypes.Subtype = res_subtypes.Subtype, results = results)
  } else {
    res = list(res_subtypes = res_subtypes, results = results)
  }

}



