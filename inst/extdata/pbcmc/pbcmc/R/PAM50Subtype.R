#'Subject \code{subtypes} for PAM50 adaptation with permuted results. 
#'
#'PAM50 subtypes are obtained using permuted test results. The idea is to 
#'give confidence in PAM50 subtype assessment (Perou et al. 2000 & 2010). 
#'In this context, the observed Spearman's \eqn{\rho} correlation is tested   
#'against the null distribution obtained for each subtype. Then, only   
#'significant correlations are used in according to the following scheme:   
#'\describe{
#'    \item{\strong{Not assigned}}{all subtype have fdr > pcutoff. Hence, 
#'         there is evidence that the observed \eqn{\rho} can be obtained by 
#'         random chance.}  
#'    \item{\strong{Assigned}}{only one fdr <= pcutoff. There is not enough
#'         evidence to say that the observed \eqn{\rho} does not belong to 
#'         the null distribution.} 
#'    \item{\strong{Ambiguous}}{more than one have fdr <= pcutoff. Then, 
#'         one of the following alternatives holds given the result of  
#'         \eqn{|\rho(profile, class_A)-\rho(profile, class_B)|>corCutoff}. 
#'         \describe{   
#'            \item{Assigned}{If the statement is TRUE.}
#'            \item{Ambiguous}{If the statement is FALSE.}    
#'        }
#'    }
#'}
#'Under the above scheme, the physician has an objective measurement to  
#'support the patient treatment decision. Both, with the given permuted   
#'subtype and by interpreting the p-value or fdr of each subtype null 
#'distribution test.   
#'
#'@param object a MolecularPermutationClassifier subclass object.   
#'@param pCutoff numeric with p-value/fdr cutoff used depending on "where"   
#'    selection. Default: 0.01.  
#'@param corCutoff numeric with correlation difference between classes 
#'    cutoff used, i.e.,  
#'    \eqn{|\rho(profile, class_A)-\rho(profile, class_B)|>corCutoff}.  
#'    Default 0.1   
#'@param where character with significant value used. Default value is "fdr".  
#'
#'@return a PAM50 object with the updated slots: 
#'\item{@@permutation}{
#'    \describe{
#'        \item{$subtype}{data.frame with the following fields}
#'        \item{$PAM50}{the original PAM50 subtype} 
#'        \item{$Permuted}{factor with the following levels:
#'            \itemize{
#'                \item "Not assigned":  all subtype have fdr > pcutoff   
#'                \item "Assigned": only one fdr <= pcutoff  
#'                \item "Ambiguous": more than one fdr <= pcutoff 
#'        }
#'     }   
#'     \item{$Classes}{a character according to "Permuted" field:  
#'        \itemize{
#'            \item the unique PAM50 subtype if "Assigned"  
#'            \item a combination for "Ambiguous" or   
#'            \item NA if "Not assigned".
#'        }
#'     }   
#'     \item{$Class}{idem as Classes but "Ambiguous" is set to PAM50 calls}  
#'     \item{$Subtype}{Classes but "Ambiguous" is kept as "Ambiguous" string}.
#'    }
#'}
#'\item{@@parameters}{$pCutoff, $corCutoff and $where are updated   
#' accordingly.}   
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include PAM50Class.R   
#'@import genefu   
#'@rdname PAM50Subtype   
#'@family PAM50   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez
#'    \email{efernandez@@bdmg.com.ar}
#'@references    
#'\enumerate{
#'    \item Perou CM, Sorlie T, Eisen MB, et al., 2000, Molecular portraits 
#'         of human breast tumors. Nature 406:747-752.  
#'    \item Perou CM, Parker JS, Prat A, Ellis MJ, Bernard PB., 2010, 
#'         Clinical implementation of the intrinsic subtypes of 
#'         breast cancer, The Lancet Oncology 11(8):718-719.  
#'}
#'@examples
#'##Using pam50centroids package example data, which already had been
#'##filtrated, classified and permutated. 
#'data(pam50centroids)
#'summary(pam50centroids)
#'
#'##Now, let's change pCutoff and corCutoff without the need to run pemutate 
#'##again
#'pam50centroids<-subtypes(pam50centroids, pCutoff=0.01, corCutoff=Inf,  
#'     where="fdr")   
#'pam50centroids    
#'summary(pam50centroids)##Note that only Basal is not Ambiguos  

setMethod(f="subtypes", signature="PAM50", definition=function(object,  
pCutoff=0.01, corCutoff=0.1,   where=c("fdr", "pvalue")[1]){   

    ##Check for permuted.pvalues  
    stopifnot("pvalues" %in% names(permutation(object)))  
    stopifnot("fdr" %in% names(permutation(object)))  
    stopifnot("correlation" %in% names(classification(object)))  
    stopifnot(where %in% c("pvalue", "fdr")) 

    ##Get the pvals to work with   
    pvals<-switch(where,    
        fdr=permutation(object)$fdr,
        pvalues=permutation(object)$pvalues
    )

    ##Get the subtypes  
    permutedSubtype<-lapply(1:nrow(pvals),function(subject){
        possibles<- pvals[subject, ] <= pCutoff    
        out<-"Not Assigned"   
        if(sum(possibles)>1){
            ##Check if they are close in correlation space: 
            ##If they are close it means that they are Ambiguous,   
            ##if not Assigned.  
            corOrder<-sort(
                classification(object)$correlation[subject, ][possibles],   
                decreasing=TRUE)
            out<-as.character(ifelse(-diff(corOrder)[1] >= corCutoff,  
                "Assigned", "Ambiguous"))   
        }else{
            if(sum(possibles)==1){
                out<-"Assigned"
            }
        }

        return(data.frame(
            PAM50=classification(object)$subtype[subject],
            Permuted=out,
            Classes=switch(out,
                "Not Assigned"=NA,   
                "Ambiguous"=toString(names(possibles)[possibles]),
                "Assigned"=colnames(pvals)[
                which.max(classification(object)$correlation[subject, ])])   
            )
        )    
    })

    ##Update the object  
    permutedSubtype<-do.call(rbind, permutedSubtype)   
    permutedSubtype$Classes<-as.character(permutedSubtype$Classes)
    permutedSubtype$Class<-permutedSubtype$Classes
    idAmbiguous<-permutedSubtype$Permuted=="Ambiguous"
    permutedSubtype$Subtype<-permutedSubtype$Classes
    ##Check for Ambiguous  
    if(any(idAmbiguous, na.rm=TRUE)){   
        permutedSubtype$Class[idAmbiguous]<-as.character(
        permutedSubtype$PAM50[idAmbiguous])
        permutedSubtype$Subtype[idAmbiguous]<-"Ambiguous"
    }
    row.names(permutedSubtype)<-row.names(permutation(object)$pvalues)    

    ##Update object slots  
    ##@permutation
    object@permutation$subtype<-permutedSubtype
    ##@parameters
    object@parameters$pCutoff<-pCutoff
    object@parameters$corCutoff<-corCutoff
    return(object)
})
