#'\code{filtrate} centroid genes from PAM50 classification   
#'
#'Remove exprs rows not required by MolecularPermutationClassifier subclasses 
#'to classify samples, in this case PAM50. This means to only keep genes with   
#'valid EntrezGeneID, i. e., not NA and present in PAM50 signature centroids. 
#'In addition, annotation slot will only keep "probe", "EntrezGene.ID" and   
#'"NCBI.gene.symbol" fields required by genefu's intrinsic.cluster.predict   
#'function.
#'
#'@param object a PAM50 object.
#'@param verbose should the user feedback be displayed? By default value is 
#'"verbose" global option parameter, if present, or FALSE otherwise.
#'
#'@return  MolecularPermutationClassifier subclass with updated slots:  
#'\item{@@exprs}{only rows required by the classifier.}   
#'\item{@@annotation}{consistent with exprs rows and only "probe",  
#'    "EntrezGene.ID" and "NCBI.gene.symbol" annotation fields.}
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include PAM50Class.R   
#'@import genefu   
#'@rdname PAM50Filtrate   
#'@family PAM50   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez
#'    \email{efernandez@@bdmg.com.ar}
#'@examples
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'pam50centroids
#'pam50centroids<-filtrate(pam50centroids, verbose=TRUE)   
#'pam50centroids
#'
#'##Using Breast Cancer NKI database, if available.  
#'if(requireNamespace("breastCancerNKI")){
#'     object<-loadBCDataset(Class=PAM50, libname="nki", verbose=TRUE) 
#'     object   
#'     object<-filtrate(object, verbose=TRUE)  
#'     object   
#'}
#'
setMethod(f="filtrate", signature="PAM50", definition=function(object,  
    verbose=getOption("verbose", default=FALSE)){   

    ##Check for data  
    stopifnot(!missing(object))
    stopifnot(nrow(exprs(object))>0)

    ##Keep probes with valid in EntrezGene.ID field.  
    if(verbose){message("Keeping only annotated EntrezGene.ID genes")}
    validIDs<-!is.na(annotation(object)$EntrezGene.ID)
    #Cannot use the setters because they check integrity, the only way is @
    object@exprs<-exprs(object)[validIDs, , drop=FALSE]      
    object@annotation<-annotation(object)[validIDs, , drop=FALSE]      
    stopifnot(validObject(object))#Check if the subset is valid   
    

    ##Look for the signature genes in the array 
    ##message("Loading PAM50 data and stored in pam50 object") 
    if(verbose){message("Keeping only PAM50 available genes.")}
    validIDs<-annotation(object)$EntrezGene.ID %in%   
    genefu::pam50$centroids.map$EntrezGene.ID
    object@exprs<-exprs(object)[validIDs, , drop=FALSE]      
    object@annotation<-annotation(object)[validIDs, , drop=FALSE]      

    ##Configuration of the annotation table for "genefu" functions. 
    annotation(object)<-annotation(object)[, c("probe", "EntrezGene.ID",  
        "NCBI.gene.symbol")]

    ##Check if the object is still valid  
    validObject(object)
    return(object)
})
