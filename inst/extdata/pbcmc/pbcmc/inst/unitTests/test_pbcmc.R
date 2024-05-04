# rm(list=ls())   
# library("RUnit")   
# library("pbcmc")   
# source("./unitTests/test_pbcmc.R")   
##-----------------------------------------------------------------------------
##MolecularPermutationClassifier-class Tests   
##-----------------------------------------------------------------------------
##Virtual test: MolecularPermutationClassifier cannot be created   
test_MolecularPermutationClassifier<-function(){
    checkException(new("MolecularPermutationClassifier"),
        msg="MolecularPermutationClassifier object cannot be created: OK.",   
    silent=TRUE)
}

##-----------------------------------------------------------------------------
##PAM50-class Tests   
##-----------------------------------------------------------------------------
##Empty object test: can PAM50 empty object be created?
test_PAM50<-function(){
    checkTrue(validObject(PAM50()),    
    msg="PAM50 empty object created: OK.")    
}

##Build PAM50 object with user data   
test_PAM50UserData<-function(){
    ##Installed Library   
    checkTrue(requireNamespace("genefu"),    
        msg="genefu library is available: OK.")
    
    ##Auxiliary function   
    test_PAM50Centroids<-function(object){
        checkTrue(!missing(object))
        checkTrue(validObject(object), msg="PAM50 using user data: OK.")   
    
        ##exprs slot   
        checkTrue(is.matrix(exprs(object)),    
            msg="PAM50Centroids exprs slot type: OK.")
        checkEquals(dim(exprs(object)), c(50, 5),  
            msg="PAM50Centroids exprs slot dimension: OK.")
    
        ##annotation slot   
        checkTrue(is.data.frame(annotation(object)),    
            msg="PAM50Centroids annotation slot type: OK.")
        checkEquals(dim(annotation(object)), c(50, 3),  
            msg="PAM50Centroids annotation slot dimension: OK.")
    
        ##target slot   
        checkTrue(is.data.frame(targets(object)),    
            msg="PAM50Centroids targets slot type: OK.")
        checkEquals(dim(targets(object)), c(0, 0),  
            msg="PAM50Centroids targets slot dimension: OK.")
    }
    
    ##Option 1: PAM50 constructor 
    M<-pam50$centroids
    genes<-pam50$centroids.map
    names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")  
    object<-PAM50(exprs=M, annotation=genes)   
    test_PAM50Centroids(object)

    ##Option 2: Two ways to build it from a MAList   
    M<-pam50$centroids
    genes<-pam50$centroids.map
    names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")  
    maux<-new("MAList", list(M=M, genes=genes))  
    test_PAM50Centroids(as(maux, "PAM50"))   
    test_PAM50Centroids(as.PAM50(maux))
}

##Getters/setters: check parameters  
test_parameters<-function(){
    data(pam50centroids)
    checkEquals(names(parameters(pam50centroids)), c("nPerm", "where",  
        "pCutoff", "corCutoff", "keep"), msg="parameters getter names: OK.")  
    checkEquals(unlist(parameters(pam50centroids)),    
        c(nPerm=10000, where="fdr", pCutoff=0.01, corCutoff=0.1, keep=TRUE),    
        msg="parameters getter: OK.")  
    aux<-parameters(pam50centroids)
    aux$keep<-FALSE
    parameters(pam50centroids)<-aux
    checkEquals(unlist(parameters(pam50centroids)),    
        c(nPerm=10000, where="fdr", pCutoff=0.01, corCutoff=0.1, 
        keep=FALSE), msg="parameters setter: OK.") 
}

##Getters/setters: check exprs  
test_exprs<-function(){
    data(pam50centroids)
    checkTrue(is.matrix(exprs(pam50centroids)), msg="exprs getter: OK.") 
    aux<-exprs(pam50centroids)
    aux[1, 1]<-NA   
    exprs(pam50centroids)<-aux
    checkTrue(is.matrix(exprs(pam50centroids)), msg="exprs setter type: OK.")
    checkTrue(is.na(exprs(pam50centroids)[1, 1]),   
        msg="exprs setter values: OK.") 
}

##Getters/setters: check annotation  
test_annotation<-function(){
    data(pam50centroids)
    checkTrue(is.data.frame(annotation(pam50centroids)),    
        msg="annotation getter: OK.")  
    aux<-annotation(pam50centroids)
    aux$probe[1]<-"test"
    annotation(pam50centroids)<-aux
    checkEquals(annotation(pam50centroids)$probe[1], "test",   
        msg="annotation setter value: OK.") 
}

##Getters/setters: check targets  
test_targets<-function(){
    data(pam50centroids)
    checkTrue(is.data.frame(targets(pam50centroids)),    
        msg="targets getter: OK.")  
    aux<-targets(pam50centroids)
    checkEquals(dim(targets(pam50centroids)), c(0, 0),  
        msg="targets dimension: OK.")  
    targets(pam50centroids)<-data.frame(foo=1:3)
    checkEquals(targets(pam50centroids), data.frame(foo=1:3),   
        msg="targets setter: OK.")  
}

##Bioconductor's dataset loading:  
##NKI example   
test_loadBCDataset<-function(){
    ##Installed Library   
    checkTrue(requireNamespace("breastCancerNKI"),    
        msg="breastCancerNKI library is available: OK.")

    ##Loading NKI dataset  
    object<-loadBCDataset(Class=PAM50, libname="nki", verbose=FALSE)  
    checkTrue(validObject(object), msg="NKI loaded PAM50 object: OK.")   

    ##exprs slot   
    checkTrue(is.matrix(exprs(object)), msg="NKI exprs slot type: OK.")   
    checkEquals(dim(exprs(object)), c(24481, 337),  
        msg="NKI exprs slot dimension: OK.")

    ##annotation slot   
    checkTrue(is.data.frame(annotation(object)),    
        msg="NKI annotation slot type: OK.")
    checkEquals(dim(annotation(object)), c(24481, 10),  
        msg="NKI annotation slot dimension: OK.")

    ##target slot   
    checkTrue(is.data.frame(targets(object)),
        msg="NKI targets slot type: OK.")
    checkEquals(dim(targets(object)), c(337, 21),  
        msg="NKI targets slot dimension: OK.")
}

##Workflow filtrate: are genes filtrated?    
test_filtrate<-function(){
    ##Loading NKI dataset  
    object<-loadBCDataset(Class=PAM50, libname="nki", verbose=FALSE)  
    object<-filtrate(object, verbose=FALSE)   
    checkTrue(validObject(object),
        msg="NKI is still valid after filtrate: OK.")  

    ##Are slots updated?  
    checkEquals(dim(exprs(object)), c(57, 337),  
        msg="NKI genes are filtrated: OK.")
    checkEquals(dim(annotation(object)), c(57, 3),  
        msg="NKI annotation are filtrated: OK.")
    checkEquals(dim(targets(object)), c(337, 21),  
        msg="NKI targets are not afected by filtrate: OK.") 
}

##Workflow classify: are subject classified?    
test_classify<-function(){
    ##Loading NKI dataset  
    object<-loadBCDataset(Class=PAM50, libname="nki", verbose=FALSE)  
    object<-filtrate(object, verbose=FALSE)   
    object<-classify(object, std="median", verbose=FALSE)  
    checkTrue(validObject(object),    
        msg="NKI is still valid after classify: OK.")  

    ##Are slots updated?  
    checkEquals(dim(exprs(object)), c(48, 337),  
        msg="NKI genes were median by classify: OK.")  
    checkEquals(dim(annotation(object)), c(48, 3),  
        msg="NKI annotation updated accordingly to median genes: OK.") 
    checkEquals(dim(targets(object)), c(337, 21),  
        msg="NKI targets are not afected by classify: OK.") 

    ##Check PAM50 results  
    ##classification slot   
    checkEquals(names(classification(object)), c("subtype", "probability",  
        "correlation"), msg="NKI classification results: OK.")
    ##classification$subtype
    aux<-as.table(c(Basal=70, Her2=55, LumA=105, LumB=76, Normal=31))
    names(attr(aux, "dimnames"))<-""   
    checkEquals(table(classification(object)$subtype), aux,   
        msg="NKI classification results: OK.") 
    ##classification$probability
    checkTrue(all(classification(object)$probability>=0) &   
        all(classification(object)$probability<=1),
        msg="NKI classification probability: OK.") 
    checkEquals(dim(classification(object)$probability), c(337, 5),  
        msg="NKI classification probability dimension: OK.")
    ##classification$correlation
    checkTrue(all(classification(object)$correlation>=-1) &   
        all(classification(object)$correlation<=1),
        msg="NKI classification correlation: OK.")     
    checkEquals(dim(classification(object)$correlation), c(337, 5),  
        msg="NKI classification correlation dimension: OK.")
}

##Workflow permutate: are genes permutated?    
test_permutate<-function(){
    ##Loading NKI dataset  
    object<-loadBCDataset(Class=PAM50, libname="nki", verbose=FALSE)  
    object<-filtrate(object, verbose=FALSE)   
    object<-classify(object, std="median", verbose=FALSE)  
    object<-permutate(object, nPerm=100, pCutoff=0.1, corCutoff=0.1, 
        keep=TRUE, verbose=FALSE)   
    checkTrue(validObject(object),    
        msg="NKI is still valid after permutate: OK.")  

    ##Check PAM50 results  
    ##permutation slot   
    checkEquals(names(permutation(object)), c("correlation", "pvalues", "fdr", 
        "subtype"), msg="NKI permutation results.") 
    ##permutation$correlation
    checkEquals(class(permutation(object)$correlation), "list",   
        msg="NKI permutation correlation class: OK.")
    checkEquals(length(permutation(object)$correlation), 337,   
        msg="NKI permutation correlation length: OK.")
    checkEquals(unique(unlist(lapply(permutation(object)$correlation,
        class))), "matrix",   
        msg="NKI permutation correlation matrix structure: OK.")   
    checkEquals(unique(unlist(lapply(permutation(object)$correlation, dim))),   
        c(100, 5), msg="NKI permutation correlation matrix dimension: OK.") 
    ##permutation$pvalues
    checkEquals(class(permutation(object)$pvalues), "matrix",   
        msg="NKI permutation pvalues class: OK.")
    checkEquals(dim(permutation(object)$pvalues), c(337, 5),  
    msg="NKI permutation pvalues dimmension: OK.")
    checkTrue(all(permutation(object)$pvalues>=0) &   
        all(permutation(object)$pvalues<=1),
        msg="NKI permutation pvalues values: OK.")
    ##permutation$fdr
    checkEquals(class(permutation(object)$fdr), "matrix",   
        msg="NKI permutation fdr class: OK.")
    checkEquals(dim(permutation(object)$fdr), c(337, 5),  
        msg="NKI permutation fdr dimmension: OK.")
    ##permutation$subtype
    checkEquals(class(permutation(object)$subtype), "data.frame",   
        msg="NKI permutation subtype class: OK.")
    aux<-classification(object)$subtype
    names(aux)<-NULL
    checkEquals(permutation(object)$subtype$PAM50, aux,   
        msg="NKI permutation subtype PAM50 field: OK.")   
    checkEquals(names(permutation(object)$subtype), c("PAM50", "Permuted",  
        "Classes", "Class",  "Subtype"), 
        msg="NKI permutation subtype names: OK.")
    checkTrue(all(c("Not Assigned", "Assigned", "Ambiguous") %in%    
        levels(permutation(object)$subtype$Permuted)),    
        msg="NKI permutation subtype levels: OK.")
    checkEquals(dim(permutation(object)$subtype), c(337, 5),  
        msg="NKI permutation subtype dimension: OK.")
}    

##Workflow summary: contingency table ok?    
test_summary<-function(){
    data(pam50centroids)
    checkEquals(as.character(class(pam50centroids)), "PAM50",   
        msg="PAM50Centroids is a PAM50 object: OK.")   
    checkTrue(all(diag(suppressMessages(summary(pam50centroids))) == c(rep(1,  
        5), 0)), msg="PAM50Centroids are classified: OK")   
}    

##Workflow subjectReport: is the plot created at least? 
test_subjectReport<-function(){
    data(pam50centroids)
    report<-subjectReport(pam50centroids, subject=1)   
    checkEquals(names(report), c("data", "layout", "plot"), 
        msg="subjectReport names check: OK.") 
    checkEquals(class(report$data), "list",   
        msg="subjectReport data class: OK.") 
    checkEquals(class(report$layout), c("Layout", "ggproto"),  
        msg="subjectReport panel class: OK.") 
    checkEquals(class(report$plot), c("gg", "ggplot"),  
        msg="subjectReport plot class: OK.") 
}

##Workflow databaseReport: is the file created?   
test_databaseReport<-function(){
    data(pam50centroids)
    fileName<-paste(tempfile(), ".pdf", sep="")  
    ##Check if a subjectReport can be obtained  
    test_subjectReport()

    ##Generate the report  
    databaseReport(pam50centroids, fileName=fileName, verbose=FALSE)  
    checkTrue(file.exists(fileName),    
        msg="databaseReport check for the pdf file: OK.")  
}

##-----------------------------------------------------------------------------
##Test functions   
##-----------------------------------------------------------------------------
##MolecularPermutationClassifier-class Tests   
# test_MolecularPermutationClassifier()   

##PAM50-class Tests   
# test_PAM50()   
# test_parameters()   
# test_exprs()   
# test_annotation()   
# test_targets()   
# test_loadBCDataset()   
# test_filtrate()   
# test_classify()   
# test_permutate()   
# test_summary()   
# test_subjectReport()   
# test_databaseReport()   
