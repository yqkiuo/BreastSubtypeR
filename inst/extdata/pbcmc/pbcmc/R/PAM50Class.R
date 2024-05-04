#'PAM50 S4 implementation in R
#'
#'This is a concrete MolecularPermutationClassifier based on Perou et al.   
#'(2000 & 2010) PAM50 molecular signature, using genefu package
#'implementation (Haibe-Kains et al. 2014).
#'
#'@section Superclassses:   
#'Direct descendant from \code{\link{MolecularPermutationClassifier-class}}. 
#'
#'@section Subclasses:   
#'None declared.   
#'
#'@slot parameters named list with at least the following fields:   
#'    \describe{
#'        \item{$nPerm}{integer with number of permutations. Default:   
#'            1e4L}
#'        \item{$where}{character with significant value used. Default   
#'             value is "fdr".} 
#'        \item{$pCutoff}{numeric with p-value or fdr cutoff used, i.e., 
#'             variable<pCutoff. Default: 0.01} 
#'        \item{$keep}{should null distribution simulation values be kept?.  
#'             Default: FALSE}  
#'        \item{corCutoff}{PAM50 additional numeric parameter with the   
#'             correlation difference between classes cutoff used, i.e., 
#'            \eqn{|\rho(profile,class_A)-\rho(profile,classB)|>corCutoff}}
#'}
#'@slot exprs matrix with gene exprs profile, where genes are in  
#'rows and subjects as columns, a.k.a., M matrix. 
#'@slot annotation data.frame with individual annotations (genes, etc). 
#'Minimal compulsory fields are: 
#'    \describe{
#'        \item{$probe}{same characters as in row.names(M).}
#'        \item{$EntrezGene.ID}{integer with NCBI Entrez Data Base.}   
#'        \item{$NCBI.gene.symbol}{character with gene mnemonic, a.k.a.    
#'             gene symbol.}  
#'}
#'@slot targets data.frame with additional subject data (optional). 
#'@slot classification named list with at least the following fields:   
#'\describe{
#'  \item{$subtype}{factor with PAM50 subtype of each sample.}
#'  \item{$probability}{matrix with the subtype probability  
#'      of each subtype per sample, as in genefu library.}  
#'  \item{$correlation}{matrix with the observed correlation  
#'      of each subtype per sample.}  
#'}            
#'@slot permutation named list with at least the following fields:   
#'\describe{
#'  \item{$correlation}{Only if keep==TRUE is a list of the   
#'      five subtypes containing a matrix with the permuted   
#'      null distribution correlations.}
#'  \item{$pvalues}{matrix with the subject's p-values of the
#'      permutation test per subject.}   
#'  \item{$fdr}{matrix with the corresponding adjusted p-values.} 
#'  \item{$subtype}{data.frame where each subject has the 
#'      reported "PAM50" subtype, the "Permuted" test result i.e.   
#'      "Assigned", "Not Assigned" or "Ambiguous"; "Classes" 
#'      whether is a single PAM50 subtype or more than one if
#'      Ambiguous case; "Class" if it is needed to assign just one
#'      i.e., a single PAM50 subtype or Not Assigned.}   
#'}
#'
#'@section Function:   
#'Redefinition from MolecularPermutationClassifier: filtrate, classify,
#'permutate, subjectReporta and databaseReport. 
#'
#'@include MolecularPermutationClassifierClass.R   
#'@docType class   
#'@name PAM50-class   
#'@rdname PAM50-class   
#'@exportClass PAM50   
#'@export PAM50   
#'@family MolecularPermutationClassifier PAM50  
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
#'         breast cancer, The Lancet Oncology 11(8):718-719  
#'}
#'
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
#'##Example 2: Build a PAM50 object with user data -------------------------   
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
#'##Example3: Work with PAM50 object: filtrate, classify and permutate--------
#'##1)Keep only annotated genes presentes in PAM50 centroids 
#'object<-filtrate(object, verbose=TRUE)   
#'
#'##2)Get PAM50 subtypes without any normalization   
#'object<-classify(object, std="none", verbose=TRUE)  
#'##Now we can inspect the how the calssification went
#'head(classification(object))
#'
#'##3)Obtain the permutation subtype 
#'##Let's run a quick example with 100 permutations. It is recommended at 
#'##least 10.000   
#'object<-permutate(object, nPerm=100, pCutoff=0.01, corCutoff=0.1, 
#'     keep=TRUE, seed=1234567890, verbose=TRUE) 
#'object
#'##Now we can inspect the how the permutation went
#'head(permutation(object))
#'##Which parameters were used? 
#'parameters(object)
#'
#'##Example 4: Obtain summary statistics and reports--------------------------  
#'##1) Let's check if we have a diagonal contigency matrix, i. e., no mistake   
#'##is made in subtype assesment.
#'summary(object)
#'
#'##2)Let's take a look at the how the patient genes behave according 
#'## to PAM50  
#'subjectReport(object, subject=1)   
#'##3)Just get a pdf with all the used subjects (PAM50 centroids in this    
#'##example).
#'#databaseReport(object, fileName="PAM50.pdf", verbose=TRUE)  

PAM50<-setClass(Class="PAM50", contains="MolecularPermutationClassifier",   
    validity=function(object){
        ##Check the integrity according to MolecularPermutationClassifier   
        ##i.e., annotation, exprs and parameters
        validSuper<-getValidity(getClassDef("MolecularPermutationClassifier"))
        stopifnot(validSuper(object))
        
        ##Parameters added: correlation cut-off 
        stopifnot(all(c("corCutoff")%in%names(parameters(object))))
        
        ##Check slot customization for PAM50 class   
        ##Global classification check  
        stopifnot(all(c("subtype","probability","correlation")%in%
        names(classification(object))))
        ##classification$subtype check   
        stopifnot(all(c("Basal","Her2","LumA","LumB","Normal")%in%
        levels(classification(object)$subtype)))
        ##classification$probability check   
        stopifnot(all(c("Basal","Her2","LumA","LumB","Normal")%in%
        colnames(classification(object)$probability)))
        ##classification$correlation correlation   
        stopifnot(all(c("Basal","Her2","LumA","LumB","Normal")%in%
        colnames(classification(object)$correlation)))
        ##Global dimensions check: all PAM50 subtypes are present and
        ##all individuals were classified if runned   
        stopifnot(ncol(classification(object)$probability)==5)
        stopifnot(ncol(classification(object)$correlation)==5)
        if(length(classification(object)$subtype)!=0){
            stopifnot(length(classification(object)$subtype)==
                ncol(exprs(object)))
            stopifnot(nrow(classification(object)$probability)==
                ncol(exprs(object)))
            stopifnot(nrow(classification(object)$correlation)==
                ncol(exprs(object)))
        }

        ##Global permutation check  
        stopifnot(all(c("correlation","pvalues","fdr","subtype")%in%
        names(permutation(object))))
        ##permutation$correlation
        ##Check only if keep==TRUE 
        if(parameters(object)$keep){
            ##Are the 5 subtypes permuted according to the parameters set?   
            stopifnot(length(permutation(object)$correlation)==
                ncol(exprs(object)))
            stopifnot(all(colnames(exprs(object))%in%
                names(permutation(object)$correlation)))
            stopifnot(unlist(unique(lapply(permutation(object)$correlation,
            dim)))==c(parameters(object)$nPerm, 5))   
        }
        ##permutation$pvalues
        stopifnot(all(c("Basal","Her2","LumA","LumB","Normal")%in%
            colnames(permutation(object)$pvalues)))
        ##permutation$fdr
        stopifnot(all(c("Basal","Her2","LumA","LumB","Normal")%in%
            colnames(permutation(object)$fdr)))
        ##permutation$subtype
        stopifnot(c("PAM50","Permuted","Classes","Class","Subtype")%in%
            names(permutation(object)$subtype))    
        ##If permutations were obtained, then check them  
        if(nrow(permutation(object)$pvalues)!=0){
            stopifnot(dim(permutation(object)$pvalues)==
                c(ncol(exprs(object)),5))
            stopifnot(dim(permutation(object)$fdr)==c(ncol(exprs(object)),5))
            stopifnot(dim(permutation(object)$subtype)==
                c(ncol(exprs(object)),5))
        }    

        return(TRUE)
    },
    prototype=list(
        parameters=list(nPerm=1e4L, where="fdr", pCutoff=0.01, keep=FALSE, 
            corCutoff=0.1),
        exprs=matrix(ncol=0,nrow=0),    
        annotation=data.frame(probe=character(), EntrezGene.ID=integer(),   
            NCBI.gene.symbol=character()),
        classification=list(
            subtype=factor(levels=c("Basal","Her2","LumA","LumB","Normal")),
            probability=matrix(nrow=0,ncol=5,dimnames=list(character(0),
                c("Basal","Her2","LumA","LumB","Normal"))),
            correlation=matrix(nrow=0,ncol=5,dimnames=list(character(0),
                c("Basal","Her2","LumA","LumB","Normal")))),
        permutation=list(correlation=list(),
            pvalues=matrix(nrow=0,ncol=5,dimnames=list(character(0),
                c("Basal","Her2","LumA","LumB","Normal"))),
            fdr=matrix(nrow=0,ncol=5,dimnames=list(character(0),
                c("Basal","Her2","LumA","LumB","Normal"))),
            subtype=data.frame(
                PAM50=factor(levels=c("Basal","Her2","LumA","LumB","Normal")),
                Permuted=factor(levels=c("Assigned",
                    "Not Assigned","Ambiguous")),   
                Classes=character(),Class=character(),Subtype=character()))
    )
)    
