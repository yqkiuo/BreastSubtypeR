#'\code{Show} a MolecularPermutationClassifier subclass object
#'
#'Basic MolecularPermutationClassifier class information display function   
#'(slots, dimensions, etc).  
#'
#'@param object an object of MolecularPermutationClassifier class hierarchy 
#'
#'@return console messages displaying the class content  
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include MolecularPermutationClassifierClass.R   
#'@exportMethod parameters   
#'@docType methods   
#'@name show   
#'@rdname MolecularPermutationClassifierShow   
#'@aliases show,MolecularPermutationClassifier-method   
#'@family MolecularPermutationClassifier   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez    
#'    \email{efernandez@@bdmg.com.ar}
#'@examples
#'##For an empty object 
#'object<-PAM50()
#'object
#'
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'pam50centroids
#'
setMethod(f="show", signature="MolecularPermutationClassifier",   
    definition=function(object){
    ##Check is not empty 
    if(nrow(exprs(object))==0){
        cat("An empty", class(object),  
            "molecular permutation classifier object\n") 
    }else{
        ##Loaded object   
        cat("A", class(object), "molecular permutation classifier object\n")   
        cat("Dimensions:\n")
        aux<-do.call(rbind, lapply(list(exprs=exprs(object),   
        annotation=annotation(object), targets=targets(object)), dim))  
        colnames(aux)<-c("nrow", "ncol")   
        show(aux)

        ##Classification already ran  
        
        if(length(classification(object)$subtype)>0){
            cat("Classification: \n")   
            ##Array items   
            aux<-as.data.frame(do.call(rbind, lapply(classification(object),   
                dim)))
            colnames(aux)<-c("nrow", "ncol")   
            show(aux)

            ##Factor items   
            classes<-unlist(lapply(classification(object), class))   
            if(any(classes=="factor")){
                aux<-lapply(classification(object)[names(classes)[classes==
                    "factor"]], table)   
                show(aux)
            }
        }

        ##Permutation already run  
        if(nrow(permutation(object)$pvalues)>0){
            cat("Permutations test ran with following parameters:\n")   
            cat(" Permutations=", parameters(object)$nPerm, ", ",
            parameters(object)$where, "<", parameters(object)$pCutoff,  
                ", corCutoff>", parameters(object)$corCutoff,  
                ", keep=", parameters(object)$keep, sep="") 

            ##Permutation
            cat("\nPermutation: \n")   
            cat("correlation available:" ,"correlation" %in% 
                names(permutation(object)),"\n")
            aux<-as.data.frame(do.call(rbind,lapply(permutation(object),
                dim)))
            colnames(aux)<-c("nrow", "ncol")   
            show(aux)
        }
    }
})
