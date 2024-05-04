#'@title MolecularPermutationClassifier high level constructor    
#'
#'@description High level constructor for MolecularPermutationClassifier   
#'subclasses using available Bioconductor's Breast Cancer example datasets. 
#'
#'@param Class name of MolecularPermutationClassifier child class to use.    
#'@param libname lowercase character with the name of the breastCancerXXX   
#'    database to be loaded. At present, XXX can be "upp", "nki", "vdx", 
#'    "mainz", "transbig" or "unt". See reference for available breast cancer   
#'    citations.
#'@param verbose should the user feedback be displayed? By default value is 
#'    "verbose" global option parameter, if present, or FALSE otherwise.
#'
#'@return MolecularPermutationClassifier subclass object with exprs,   
#'annotation and targets slots taken from the libname used.
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include MolecularPermutationClassifierClass.R   
#'@importMethodsFrom Biobase fData pData 
#'@importFrom methods as new validObject
#'@exportMethod loadBCDataset   
#'@docType methods   
#'@name loadBCDataset   
#'@rdname MolecularPermutationClassifierConstructor   
#'@aliases loadBCDataset-methods   
#'@family MolecularPermutationClassifier PAM50  
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez    
#'    \email{efernandez@@bdmg.com.ar}
#'@references
#'\itemize{
#'    \item Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G and 
#'         Quackenbush J (2011). breastCancerUPP: Gene expression dataset 
#'         published by Miller et al. [2005] (UPP).. R package version 1.3.1, 
#'         \url{http://compbio.dfci.harvard.edu/}.   
#'    \item Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G and 
#'         Quackenbush J (2011). breastCancerNKI: Genexpression dataset  
#'         published by van't Veer et al. [2002] and van de Vijver et al.   
#'         [2002] (NKI).. R package version 1.3.1,  
#'         \url{http://compbio.dfci.harvard.edu/}.   
#'    \item Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G and 
#'         Quackenbush J (2011). breastCancerVDX: Gene expression datasets 
#'         published by Wang et al. [2005] and Minn et al. [2007] (VDX). R   
#'         package version 1.3.1, \url{http://compbio.dfci.harvard.edu/}.
#'    \item Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G and 
#'         Quackenbush J (2011). breastCancerTRANSBIG: Gene expression dataset 
#'         published by Desmedt et al. [2007] (TRANSBIG).. R package version  
#'         1.3.1, \url{http://compbio.dfci.harvard.edu/}.  
#'    \item Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G and 
#'         Quackenbush J (2011). breastCancerMAINZ: Gene expression dataset 
#'         published by Schmidt et al. [2008] (MAINZ).. R package version  
#'         1.3.1, \url{http://compbio.dfci.harvard.edu/}.  
#'    \item Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G and 
#'         Quackenbush J (2011). breastCancerUNT: Gene expression dataset 
#'         published by Sotiriou et al. [2007] (UNT).. R package version  
#'         1.3.1, \url{http://compbio.dfci.harvard.edu/}.  
#'}
#'
#'@examples
#'##Using Breast Cancer NKI database, if available, to create a PAM50 class. 
#'if(requireNamespace("breastCancerNKI")){
#'     object<-loadBCDataset(Class=PAM50, libname="nki", verbose=TRUE) 
#'     object   
#'
#'     ##Now we can inspect the object  
#'     head(exprs(object))      ##The gene expression   
#'     head(annotation(object)) ##The available annotation
#'     head(targets(object))    ##The clinical data present in the package 
#'}
#'
setGeneric(name="loadBCDataset", def=function(Class, libname=c("upp", "nki", 
    "vdx", "mainz", "transbig", "unt"), verbose=getOption("verbose",    
    default=FALSE)){
    standardGeneric("loadBCDataset")
})
#'
#'@name loadBCDataset   
#'@rdname MolecularPermutationClassifierConstructor   
#'@inheritParams loadBCDataset   
#'@aliases loadBCDataset,classGeneratorFunction-method   
setMethod(f="loadBCDataset", signature="classGeneratorFunction",   
definition=function(Class, libname=c("upp", "nki", "vdx", "mainz",
    "transbig", "unt"), verbose=getOption("verbose", default=FALSE)){ 

    ##Check database   
    stopifnot(!missing(libname))
    stopifnot(length(libname)==1)
    stopifnot(libname %in% c("upp", "nki", "vdx", "mainz", "transbig", "unt")) 

    ##Load the appropriate dataset and get its dataset name
    datasetName<-switch(libname,
        upp={requireNamespace("breastCancerUPP")    
            c(datum="upp", package="breastCancerUPP")},   
        nki={requireNamespace("breastCancerNKI")
            c(datum="nki", package="breastCancerNKI")},   
        vdx={requireNamespace("breastCancerVDX")
            c(datum="vdx", package="breastCancerVDX")},   
        transbig={requireNamespace("breastCancerTRANSBIG")
            c(datum="transbig",  package="breastCancerTRANSBIG")},  
        mainz={requireNamespace("breastCancerMAINZ")
            c(datum="mainz", package="breastCancerMAINZ")},   
        unt={requireNamespace("breastCancerUNT")
            c(datum="unt", package="breastCancerUNT")}   
    )

    ##Data unification: annotation, clinical and exprs.   
    ##Load into memory and remove original variable  
    if(verbose){message("Loading dataset...")}   
    data(list=datasetName["datum"],    
        package=datasetName["package"],    
        envir=environment())
    dataset<-eval(parse(text=datasetName["datum"]))
    rm(list=datasetName["datum"], envir=environment())   
    
    ##Get annotation: probe, entrez and symbol   
    if(verbose){message("Building a ", Class@className,  " object")}  
    annotData <-fData(dataset)   
    ##Check ProbeName field  
    probe<-c("probe", "ProbeName")   
    stopifnot(any(probe %in% names(annotData)))  
    probe<-probe[probe %in% names(annotData)][1]  
    ##Check EntrezID field  
    entrez<-c("EntrezGene.ID", "EntrezID")   
    stopifnot(any(entrez %in% names(annotData)))  
    entrez<-entrez[entrez %in% names(annotData)][1]  
    ##Check NCBI.gene.symbol field  
    symbol<-c("Gene.symbol", "HUGO.gene.symbol")   
    stopifnot(any(symbol %in% names(annotData)))  
    symbol<-symbol[symbol %in% names(annotData)][1]  

    ##Finally assign the fields 
    annotData$probe <- as.character(annotData[,probe])  
    annotData$EntrezGene.ID  <- as.integer(as.character(annotData[,entrez])) 
    annotData$NCBI.gene.symbol <- as.character(annotData[, symbol]) 

    ##Get clinical and experimental data
    clinData <- pData(dataset)  
    expData <- exprs(dataset)  
    
    ##Create PAM50 object  
    object<-Class()
    object@exprs<-expData
    object@annotation<-annotData
    object@targets<-clinData
    validObject(object)

    return(object)
})
