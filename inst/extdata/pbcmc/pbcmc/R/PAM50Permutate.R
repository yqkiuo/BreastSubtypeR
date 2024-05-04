#'\code{permutate} subject gene-expression for PAM50 confidence   
#'
#'Calculate the null Spearman's \eqn{\rho} distribution of each subtype by   
#'means of gene label permutation, in order to evaluate if the observed 
#'values could be obtained by random change.  
#'
#'@param object a MolecularPermutationClassifier subclass object.   
#'@param nPerm integer with number of permutations. Default: 1e4L
#'@param pCutoff numeric with p-value or fdr cutoff used, i.e.,   
#'    variable<pCutoff. Default: 0.01  
#'@param where character with significant value used. Default value is "fdr".  
#'@param keep should null distribution simulation values be kept?.
#'    Default: FALSE   
#'@param corCutoff numeric with correlation difference between classes 
#'    cutoff used, i.e.,  
#'    \eqn{|\rho(profile, class_A)-\rho(profile, class_B)|>corCutoff}.  
#'    Default 0.1   
#'@param seed integer to use as random seed. Default: 1234567890.   
#'@param BPPARAM an optional BiocParallelParam instance determining the 
#'parallel back-end to be used during evaluation, or a list of  
#'BiocParallelParam instances, to be applied in sequence for nested calls to  
#'bplapply. Default=bpparam().   
#'@param verbose should the user feedback be displayed? By default value is 
#'   "verbose" global option parameter, if present, or FALSE otherwise. 
#'
#'@return a PAM50 object with the following updated slots:
#'\item{@@permutation}{
#'    \describe{
#'        \item{$pvalues}{numeric matrix with subtype pvalues obtained as the 
#'             number of times the permuted correlation is greater or equal  
#'             the observed correlation divided the number of permutations.}
#'        \item{$fdr}{subtype adjusted pvalues for each subject with False 
#'             Discovery Rate.}  
#'        \item{$correlations}{list with subject matrix correlation of each  
#'             permutation simulation.}  
#'        \item{$subtype}{data.frame with classification results 
#'              obtained by subtype function.}   
#'    }
#'}
#'\item{@@parameters}{$nPerm, $pCutoff, $where and $keep updated accordingly.}  
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include PAM50Class.R   
#'@import genefu parallel  
#'@rdname PAM50Permutate   
#'@family PAM50   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez
#'    \email{efernandez@@bdmg.com.ar}
#'@references    
#'\enumerate{
#'    \item Haibe-Kains B, Schroeder M, Bontempi G, Sotiriou C and   
#'         Quackenbush J, 2014, genefu: Relevant Functions for Gene
#'         Expression Analysis, Especially in Breast Cancer. R package
#'         version 1.16.0, \url{www.pmgenomics.ca/bhklab/} 
#'    \item Perou CM, Sorlie T, Eisen MB, et al., 2000, Molecular portraits 
#'         of human breast tumors. Nature 406:747-752.  
#'    \item Perou CM, Parker JS, Prat A, Ellis MJ, Bernard PB., 2010, 
#'         Clinical implementation of the intrinsic subtypes of 
#'         breast cancer, The Lancet Oncology 11(8):718-719.  
#'}
#'@examples
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'pam50centroids
#'pam50centroids<-filtrate(pam50centroids, verbose=TRUE)   
#'pam50centroids<-classify(pam50centroids, std="none", verbose=TRUE)  
#'
#'##Let's run a quick example with 100 permutations. It is recommended at 
#'##least 10.000   
#'pam50centroids<-permutate(pam50centroids, nPerm=100, pCutoff=0.01,  
#'corCutoff=0.1, verbose=TRUE)   
#'pam50centroids
setMethod(f="permutate", signature="PAM50", definition=function(object,  
    nPerm=1e4, pCutoff=0.01, where="fdr", keep=FALSE, corCutoff=0.1,
    seed=1234567890, BPPARAM=bpparam(),   
    verbose=getOption("verbose",default=TRUE)){

    ##Check object   
    stopifnot(nrow(exprs(object))>0)
    if(length(classification(object)$subtype)==0){
        stop(paste("No call to classify(object, ...) has been made yet.",
            "Please run classify(object)."))  
    }
    
    ##PAM50 parameters   
    pam50.aux<-genefu::pam50
    pam50.aux$std<-"none"
    set.seed(seed)

    ##Parallel permutations   
    if(verbose){
        message("Obtaining ", nPerm, " permutations for ",ncol(exprs(object)),  
            " subjects..." )  
        bpprogressbar(BPPARAM)<-TRUE
    }

    ##For each subject obtains its permutations   
    out<-bplapply(1:ncol(exprs(object)), function(subject, dataset, annot, 
        model){
        ##Obtain the permutation  
        genePermutations<-do.call(rbind, lapply(1:nPerm,   
            function(permutation){
            return(dataset[subject, sample(ncol(dataset))] )  
        }))
        colnames(genePermutations)<-colnames(dataset)

        ##Get the correlation  
        corPermutations<-intrinsic.cluster.predict(sbt.model=model,    
            data=genePermutations, annot=annot, do.mapping=TRUE,  
            do.prediction.strength=FALSE, verbose=FALSE)$cor   

        ##Obtain the p-value  
        pvals<-rowSums(t(corPermutations)>=    
        classification(object)$correlation[subject, , drop=TRUE])/nPerm  

        ##Output
        if(keep){
            return(list(pvals=pvals, permutedCor=corPermutations))   
        }else{ 
            return(pvals)
        }
    }, dataset=t(exprs(object)), annot=annotation(object), model=pam50.aux, 
    BPPARAM=BPPARAM)

    ##For each subject  
    ##Permutations completed feedback  
    if(verbose){
        message("Obtaining ", nPerm, " permutations for ",  
            ncol(exprs(object)), " subjects... done." )
    }

    ##Format data   
    if(!keep){
        pvals<-do.call(rbind, out)   
    }else{
        pvals<-do.call(rbind, lapply(out,function(x){x$pvals}))   
        permutedCor<-lapply(out,function(x){x$permutedCor})
        names(permutedCor)<-colnames(exprs(object))
        ##Update object   
        object@permutation$correlation<-permutedCor
    }
    row.names(pvals)<-colnames(exprs(object))

    ##Update object   
    object@permutation$pvalues<-pvals
    object@permutation$fdr<-t(apply(pvals, 1, p.adjust, method="fdr")) 
    object@parameters$pCutoff<-pCutoff
    object@parameters$nPerm<-nPerm
    object@parameters$keep<-keep
    object@parameters$where<-where
    ##Get permuted subtypes calls 
    object <- subtypes(object, pCutoff=pCutoff, corCutoff=corCutoff,    
        where=where)

    validObject(object)
    return(object)
})
