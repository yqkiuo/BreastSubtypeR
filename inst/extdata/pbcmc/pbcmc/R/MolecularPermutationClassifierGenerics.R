#'Virtual functions for MolecularPermutationClassifier hierarchy
#'
#'The following functions establish an organized framework for 
#'MolecularPermutationClassifier subclasses data processing. In this context,  
#'the later are supposed to be implemented with respective responsibilities.   
#'In particular, once the class is created the user has to:  
#'\describe{
#'    \item{filtrate:}{Removes, from the exprs matrix, subjects not required 
#'         by the classification algorithm.}
#'    \item{classify:}{Generates subject classification according to
#'         subclass implementations (PAM50, etc.).}
#'    \item{permute:}{Obtains subject classification based on the null  
#'         correlation distribution by means permutation simulation.}  
#'    \item{subtype:}{Obtaind the new classification using permutation   
#'         results.}   
#'    \item{subjectReport:}{A friendly report for physician treatment   
#'         decision support.}  
#'    \item{databaseReport:}{A pdf with all subjectReports, if a database is
#'         available.}   
#'}
#'
#'@param object MolecularPermutationClassifier child class object   
#'@param verbose should the user feedback be displayed? By default value  
#'    is "verbose" global option parameter, if present, or FALSE otherwise.   
#'@param nPerm integer with number of permutations. Default: 1e4L.
#'@param pCutoff numeric with p-value or fdr cutoff used, i.e.,   
#'    variable<pCutoff. Default: 0.01.  
#'@param where character with significant value used. Default value is "fdr".  
#'@param keep should null distribution simulation values be kept?.
#'    Default: FALSE   
#'@param seed integer to use as random seed. Default: 1234567890.   
#'@param BPPARAM an optional BiocParallelParam instance determining the 
#'parallel back-end to be used during evaluation, or a list of  
#'BiocParallelParam instances, to be applied in sequence for nested calls to  
#'bplapply. Default=bpparam().   
#'@param subject integer to select the appropriate subject to report.   
#'@param fileName character with the name of the pdf report file to save.
#'@param ... additional parameters for future implementations.  
#'
#'@return A MolecularPermutationClassifier child according to the actual 
#'object class.   
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include MolecularPermutationClassifierClass.R   
#'@exportMethod filtrate   
#'@docType methods   
#'@name filtrate   
#'@rdname MolecularPermutationClassifierGenerics   
#'@aliases filtrate-methods   
#'@family MolecularPermutationClassifier PAM50  
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez    
#'    \email{efernandez@@bdmg.com.ar}
#'
#'@examples
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'pam50centroids
#'pam50centroids<-filtrate(pam50centroids, verbose=TRUE)   
#'pam50centroids<-classify(pam50centroids, std="none", verbose=TRUE)  
#'##Let's run a quick example with 100 permutations. It is recommended at 
#'##least 10.000   
#'pam50centroids<-permutate(pam50centroids, nPerm=100, pCutoff=0.01,  
#'corCutoff=0.1, verbose=TRUE)   
#'pam50centroids
#'
setGeneric(name="filtrate", def=function(object, verbose=getOption("verbose",  
default=FALSE)){
    standardGeneric("filtrate")
})
#'@name classify   
#'@exportMethod classify   
#'@rdname MolecularPermutationClassifierGenerics   
#'@inheritParams filtrate   
#'@aliases classify-methods   
setGeneric(name="classify", def=function(object, ...,  
verbose=getOption("verbose", default=FALSE)){   
    standardGeneric("classify")
})
#'@name permutate   
#'@exportMethod permutate   
#'@rdname MolecularPermutationClassifierGenerics   
#'@inheritParams filtrate   
#'@importFrom BiocParallel bplapply bpparam bpprogressbar<-
#'@aliases permutate-methods   
setGeneric(name="permutate", def=function(object, nPerm=1e4L, pCutoff=0.01, 
where="fdr", keep=FALSE, ... , seed=1234567890, BPPARAM=bpparam(),   
verbose=getOption("verbose", default=TRUE)){   
    standardGeneric("permutate")
})
#'@name subtypes   
#'@exportMethod subtypes   
#'@rdname MolecularPermutationClassifierGenerics   
#'@inheritParams filtrate   
#'@aliases subtypes-methods   
setGeneric(name="subtypes", def=function(object, pCutoff=0.01, ..., 
where=c("fdr", "pvalue")[1]){   
    standardGeneric("subtypes")
})
#'@name subjectReport   
#'@exportMethod subjectReport   
#'@rdname MolecularPermutationClassifierGenerics   
#'@inheritParams filtrate   
#'@aliases subjectReport-methods   
setGeneric(name="subjectReport", def=function(object, subject){  
    standardGeneric("subjectReport")
})
#'@name databaseReport   
#'@exportMethod databaseReport   
#'@rdname MolecularPermutationClassifierGenerics   
#'@inheritParams filtrate   
#'@aliases databaseReport-methods   
setGeneric(name="databaseReport", def=function(object, fileName, ..., 
verbose=getOption("verbose", default=TRUE)){   
    standardGeneric("databaseReport")
})
