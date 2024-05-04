#'\code{PAM50} high level coerce functions
#'
#'These functions (setAs and as.PAM50) are intended to be used with limma 
#'\code{\link{MAList-class}} in order to coerce its structure into a    
#'compatible PAM50 class.  
#'
#'Basically the $M and $genes items are copied into a   
#'MolecularPermutationClassifier's exprs and annotation slots respectively.   
#'In addition, if present, $targets content is also copied to the same named
#'slot.
#'
#'@param object MAList object with at least M and genes items, optionally 
#'targets.    
#'@param Class character with the name of class "PAM50" to be coerced. 
#'@param strict,ext see \code{\link{as}} function.
#'
#'@return a PAM50 object with the respective copied data.
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include PAM50Class.R   
#'@importClassesFrom limma MAList  
#'@importFrom limma show  
#'@importMethodsFrom methods coerce  
#'@name as   
#'@rdname PAM50Constructor   
#'@usage as(object,Class,strict=TRUE,ext=possibleExtends(thisClass,Class))   
#'@family PAM50   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez
#'    \email{efernandez@@bdmg.com.ar}
#'@export
#'@examples
#'##Example 1: Create a PAM50 object -----------------------------------------  
#'##1) Just an empty object
#'object<-PAM50()
#'object
#'
#'##2) Using Breast Cancer NKI database, if available. 
#'if(requireNamespace("breastCancerNKI")){
#'     object<-loadBCDataset(Class=PAM50, libname="nki", verbose=TRUE) 
#'     object   
#'     ##Now we can inspect the object  
#'     head(exprs(object))      ##The gene expression   
#'     head(annotation(object)) ##The available annotation
#'     head(targets(object))    ##The clinical data present in the package 
#'}
#'
#'##Example 2: Build a PAM50 object with user data --------------------------   
#'##Option 1: using PAM50 constructor. The user will only need:   
#'##a) The M gene expression object, i. e., gene in rows and sample in columns  
#'##b) The annotation data.frame which must include the compulsory fields   
#'## "probe", "NCBI.gene.symbol" and "EntrezGene.ID"
#'M<-pam50$centroids
#'genes<-pam50$centroids.map
#'names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")  
#'object<-PAM50(exprs=M, annotation=genes)   
#'object
#'
#'##Option 2: Two ways to build it from a MAList (as or as.PAM50)-------------
#'##Let's use PAM50 classifier's centroids toy example, i. e., the five subject 
#'##subtypes, which must correctly classify all the subject. 
#'M<-pam50$centroids
#'genes<-pam50$centroids.map
#'names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")  
#'maux<-new("MAList", list(M=M, genes=genes))  
#'##calling as function  
#'object<-as(maux, "PAM50")   
#'object
#'##same result with as.PAM50 function
#'object<-as.PAM50(maux)
#'object
#'
setAs(from="MAList", to="PAM50", def=function(from){  
    ##Create the new class 
    to<-PAM50(exprs=from$M, annotation=from$genes)   

    ##Check for targets (optional slot)
    if("targets" %in% names(from)){  
        to@targets<-from$targets
    }

    ##Check if it is valid and return  
    validObject(to)
    return(to)
})
#'
#'@exportMethod as.PAM50   
#'@docType methods   
#'@name as.PAM50   
#'@rdname PAM50Constructor   
#'@inheritParams setAs   
#'@aliases as.PAM50-methods   
setGeneric(name="as.PAM50", def=function(object){   
    standardGeneric("as.PAM50")
})
#'
#'@name as.PAM50   
#'@rdname PAM50Constructor   
#'@inheritParams setAs   
#'@aliases as.PAM50,MAList-method   
setMethod(f="as.PAM50", signature="MAList",   
    definition=function(object){
    as(object, "PAM50")   
})
