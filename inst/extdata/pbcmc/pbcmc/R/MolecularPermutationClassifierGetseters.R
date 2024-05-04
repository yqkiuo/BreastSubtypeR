#'Accessors for MolecularPermutationClassifier child class slots   
#'
#'Slot setters/getters for MolecularPermutationClassifier hierarchy classes   
#'    
#'@param object MolecularPermutationClassifier subclass object
#'@param value according to the function call:  
#'\itemize{
#'    \item parameters: named list with at least the following fields:   
#'        \describe{
#'            \item{$nPerm}{integer with number of permutations. Default:   
#'                1e4L}
#'            \item{$where}{character with significant value used. Default   
#'                value is "fdr".}  
#'            \item{$pCutoff}{numeric with p-value or fdr cutoff used, i.e., 
#'                variable<pCutoff. Default: 0.01}  
#'            \item{$corCutoff}{numeric with correlation difference between
#'                classes cutoff used, i.e., 
#'            \eqn{|\rho(profile,class_A)-\rho(profile,class_B)|>corCutoff}}
#'            \item{$keep}{should null distribution simulation values be   
#'                kept?. Default: FALSE}  
#'        }
#'    \item annotation: data.frame with individual annotations (genes, etc). 
#'         Minimal compulsory fields are:
#'    \describe{
#'        \item{$probe}{same characters as in row.names(M).}
#'        \item{$EntrezGene.ID}{integer with NCBI Entrez Data Base.}   
#'        \item{$NCBI.gene.symbol}{character with gene mnemonic, 
#'        a.k.a. gene symbol.}  
#'    }
#'    \item exprs: matrix with gene exprs profile, where genes are in rows 
#'         and subjects as columns, a.k.a., \strong{M matrix}. 
#'    \item targets: data.frame with additional subject data.  
#'}
#'@param ... additional parameters according to function call. 
#'
#'@return according to function call one of the following objects:   
#'\item{parameters}{named list see value parameter}
#'\item{exprs}{matrix with gene exprs profile, where genes are in rows and  
#'     subjects as columns, a.k.a., \strong{M matrix}.}  
#'\item{annotation}{data.frame see value parameter} 
#'\item{classification}{named list with at least the following fields:} 
#'    \describe{
#'        \item{$class}{factor with with all possible class levels.}  
#'    }
#'\item{permutation}{named list with at least the following fields:} 
#'    \describe{
#'        \item{pvalues}{numeric matrix with subjects in row and classes 
#'             in columns.}  
#'        \item{$fdr}{numeric matrix with False Discovery Rate correction  
#'             of pvalues by row.}
#'    }
#'\item{parameters<-}{MolecularPermutationClassifier object with parameters 
#'     updated slot.}  
#'\item{exprs<-}{MolecularPermutationClassifier object with exprs updated
#'     slot.}   
#'\item{annotation<-}{MolecularPermutationClassifier object with annotation 
#'     updated slot.}  
#'\item{targets<-}{MolecularPermutationClassifier object with targets 
#'     updated slot.}  
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include MolecularPermutationClassifierClass.R   
#'@exportMethod parameters   
#'@docType methods   
#'@name parameters   
#'@rdname MolecularPermutationClassifierGetseters   
#'@aliases parameters-methods   
#'@importMethodsFrom BiocGenerics annotation "annotation<-" as.data.frame
#'cbind colnames "colnames<-" do.call eval lapply ncol nrow order paste   
#'rbind rownames "rownames<-" sort subset table unique unlist 
#'@family MolecularPermutationClassifier   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez
#'    \email{efernandez@@bdmg.com.ar}
#'@examples
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'
#'##Now we can inspect pam50centroids object   
#'head(exprs(pam50centroids))      ##The gene expression
#'head(annotation(pam50centroids)) ##The available annotation 
#'head(targets(pam50centroids))    ##The clinical data present in the package  
#'
#'##Work with the parameters 
#'parameters(pam50centroids)       ##Display them
#'aux<-parameters(pam50centroids)    
#'aux$keep<-TRUE                   ##Set keep to FALSE  
#'parameters(pam50centroids)<-aux
#'parameters(pam50centroids)        
#'
#'##Also exprs<-, annotation<- and targets<- available functions to update    
#'##the respective slots  

setGeneric(name="parameters", def=function(object){   
    standardGeneric("parameters")
})
#'
#'@name parameters   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases parameters,MolecularPermutationClassifier-method   
setMethod(f="parameters", signature="MolecularPermutationClassifier",   
    definition=function(object){
    return(object@parameters)
})
#'
#'@exportMethod parameters<-   
#'@docType methods   
#'@name parameters<-   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases parameters<--methods   
setGeneric(name="parameters<-", def=function(object, value){  
    standardGeneric("parameters<-")
})
#'
#'@name parameters<-   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases parameters<-,MolecularPermutationClassifier-method   
setReplaceMethod(f="parameters", signature="MolecularPermutationClassifier",   
    definition=function(object, value){   
    object@parameters<-value
    validObject(object)
    return(object)
})
#'
#'@importFrom Biobase exprs  
#'@exportMethod exprs   
#'@name exprs   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases exprs,MolecularPermutationClassifier-method   
setMethod(f="exprs", signature="MolecularPermutationClassifier",   
    definition=function(object){
    return(object@exprs)
})
#'
#'@importFrom Biobase exprs<-  
#'@name exprs<-   
#'@exportMethod exprs<-   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases exprs<-,MolecularPermutationClassifier,ANY-method   
setReplaceMethod(f="exprs", signature="MolecularPermutationClassifier",   
    definition=function(object, value){   
    object@exprs<-value
    validObject(object)
    return(object)
})
#'
#'@importFrom BiocGenerics annotation  
#'@name annotation   
#'@exportMethod annotation   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases annotation,MolecularPermutationClassifier-method   
setMethod(f="annotation", signature="MolecularPermutationClassifier",   
    definition=function(object, ...){   
    return(object@annotation)
})
#'
#'@importFrom BiocGenerics annotation<-  
#'@name annotation<-   
#'@exportMethod annotation<-   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases annotation<-,MolecularPermutationClassifier,ANY-method   
setReplaceMethod(f="annotation", signature="MolecularPermutationClassifier",   
    definition=function(object, value){   
    object@annotation<-value
    validObject(object)
    return(object)
})
#'@name targets   
#'@exportMethod targets   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases targets-methods   
setGeneric(name="targets", def=function(object){   
    standardGeneric("targets")
})
#'
#'@name targets   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases targets,MolecularPermutationClassifier-method   
setMethod(f="targets", signature="MolecularPermutationClassifier",   
    definition=function(object){
    return(object@targets)
})
#'
#'@name targets<-   
#'@exportMethod targets<-   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases targets<--methods   
setGeneric(name="targets<-", def=function(object, value){  
    standardGeneric("targets<-")
})
#'
#'@name targets<-   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases targets<-,MolecularPermutationClassifier-method   
setReplaceMethod(f="targets", signature="MolecularPermutationClassifier",   
    definition=function(object, value){   
    object@targets<-value
    validObject(object)
    return(object)
})
#'
#'@name classification   
#'@exportMethod classification   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases classification-methods   
setGeneric(name="classification", def=function(object){   
    standardGeneric("classification")
})
#'
#'@name classification   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases classification,MolecularPermutationClassifier-method   
setMethod(f="classification", signature="MolecularPermutationClassifier",   
    definition=function(object){
    return(object@classification)
})
#'
#'@name permutation   
#'@exportMethod permutation   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases permutation-methods   
setGeneric(name="permutation", def=function(object){   
    standardGeneric("permutation")
})
#'
#'@name permutation   
#'@rdname MolecularPermutationClassifierGetseters   
#'@inheritParams parameters   
#'@aliases permutation,MolecularPermutationClassifier-method   
setMethod(f="permutation", signature="MolecularPermutationClassifier",   
    definition=function(object){
    return(object@permutation)
})
