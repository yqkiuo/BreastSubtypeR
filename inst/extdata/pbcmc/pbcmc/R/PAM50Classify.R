#'\code{classify} subjects with PAM50 molecular signature   
#'
#'Obtain PAM50 subtype using genefu centroid Spearman's correlation 
#'implementation. If std=="median" probes with the same mapping are averaged.   
#'Then, the complete database is center normalized using gene median   
#'expression. This is done in order to assure selecting the same "gene" to
#'those in "genefu" library, instead of the most variant probe (default in 
#'geneid.map), when more than one probe match the same gene. This selection 
#'is based on probe population variance that could depends on the number of
#'accounted genes.   
#'
#'@param object a MolecularPermutationClassifier subclass object.   
#'@param std character to select standardization alternative "none" 
#'    (default), "scale" and "robust" as in genefu original implementation,
#'    plus the suggested "median" if many subjects are available.
#'@param verbose should the user feedback be displayed? By default value is 
#'    "verbose" global option parameter, if present, or FALSE otherwise.
#'
#'@return a PAM50 object with the updated slots: 
#'\item{@@exprs}{updated matrix with the used std parameter.}  
#'\item{@@classification}{
#'    \describe{
#'        \item{$subtype}{subject named factor with all classifier possible  
#'             levels, i.e, "Basal", "Her2", "LumA", "LumB" and "Normal".}
#'        \item{$probability}{numeric matrix with subtype class probability   
#'             for each subject, as in genefu, obtained as the positive  
#'             proportion of correlation explained by each subtype.} 
#'        \item{$correlation}{numeric matrix with Spearman's rho
#'             correlation of each subject to the corresponding PAM50
#'             subtypes.}   
#'    }
#'}
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include PAM50Class.R   
#'@importFrom limma avereps  
#'@import genefu   
#'@importFrom stats median na.omit p.adjust
#'@rdname PAM50Classify   
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
#'
#'@examples
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'
#'##Get the original PAM50 calls using genefu implementation 
#'pam50centroids<-classify(pam50centroids, std="none", verbose=TRUE)  
#'classification(pam50centroids)

setMethod(f="classify", signature="PAM50", definition=function(object,  
    std=c("none", "scale", "robust", "median")[1], verbose=getOption("verbose",
    default=FALSE)){

    ##Check for standardization method 
    stopifnot(std %in% c("median", "none", "scale", "robust"))   
    stopifnot(length(std)==1)
    ##Force filtrate if not already performed   
    if(verbose){
        message("Enforcing filtrate(object)...")   
    }
    object<-filtrate(object, verbose=verbose)   

    ##auxiliary pam   
    pam50.aux <- genefu::pam50  
    pam50.aux$std<-std

    ##dataset center scale using median estimation using all subjects
    ##present.
    if(std=="median"){
        ##Check for single subject 
        if(ncol(exprs(object))==1){
            stop("Not posible to use std=\"median\" for a single subject.")}
        ##Check for appropiate annotation i. e. EntrezGene.ID-Symbol tuples 
        ##for the same EntrezGene.ID have the same Symbol 
        if(verbose){
        message("Annotation check over identical EntrezGeneID probes...")}   
        if(nrow(unique(annotation(object)[,-1])) !=   
            length(unique(annotation(object)$EntrezGene.ID))){
            ##Search for the wrong tuples
            ids<-unique(annotation(object)$EntrezGene.ID[
                duplicated(annotation(object)$EntrezGene.ID)])
            check<-sapply(ids, function(posible){   
                nrow(unique(annotation(object)[
                    annotation(object)$EntrezGene.ID==posible, -1]))!=1   
            })
            stop(paste("Inconsistent annotation for EntrezGene.ID:", 
                toString(ids[check]),
                ". Different symbols are present for the same EntrezGene.ID.",
                " Please, check and correct the annotation."))  
        }
                
        ##Average probes with the same EntrezGene.ID   
        if(verbose){
            message("Averaging over identical EntrezGeneID probes...")}
        object@exprs<-avereps(exprs(object),    
            ID=as.character(annotation(object)$EntrezGene.ID))
        ##Keep class consistency  
        object@annotation<-unique(annotation(object)[, c("EntrezGene.ID",   
                "NCBI.gene.symbol")])
        object@annotation$probe<-row.names(annotation(object))
        object@annotation<-annotation(object)[, c("probe",   
                "EntrezGene.ID", "NCBI.gene.symbol")]   
        mask<-annotation(object)$probe
        names(mask)<-as.character(annotation(object)$EntrezGene.ID)
        row.names(object@exprs)<-mask[row.names(object@exprs)]
        stopifnot(validObject(object))

        ##Median scaling   
        if(verbose){message("Center scaling using gene median...")}
        med.total<-apply(exprs(object), 1, median, na.rm=TRUE) 
        exprs(object)<-t(scale(t(exprs(object)), center=med.total,   
            scale=FALSE))
        row.names(exprs(object)) <- annotation(object)$probe  
        dataset<-t(exprs(object))
        pam50.aux$std<-"none"
    }else{
        dataset<-t(exprs(object))
    }

    ##Get the PAM50 calls 
    if(verbose){message("Getting PAM50 subtypes...")}  
    subtype.norm <- intrinsic.cluster.predict(sbt.model=pam50.aux,  
        data=dataset, annot=annotation(object), do.mapping=TRUE,  
        do.prediction.strength=FALSE, verbose=verbose)   

    ##Update object   
    if(std != "median"){  
        object@exprs<-t(subtype.norm$profiles)
        #keep the same individuals and order   
        object@annotation<-annotation(object)[
            row.names(annotation(object)) %in% row.names(exprs(object)), ]     
        object@exprs<-exprs(object)[order(row.names(exprs(object))),
            ,drop=FALSE]    
        object@annotation<-annotation(object)[order(
            row.names(annotation(object))), ]   
        stopifnot(validObject(object))
    }

    ##Obtained classification   
    object@classification$subtype <- factor(subtype.norm$subtype,  
        levels=c("Basal", "Her2", "LumA", "LumB", "Normal"))
    names(object@classification$subtype) <- colnames(exprs(object))  
    object@classification$probability <-subtype.norm$subtype.proba   
    object@classification$correlation <- subtype.norm$cor  

    validObject(object)
    return(object)
})
