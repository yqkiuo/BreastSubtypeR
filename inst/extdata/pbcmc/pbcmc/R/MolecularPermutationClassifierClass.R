#'Class \code{MolecularPermutationClassifier} S4 implementation in R   
#'
#'Virtual class to represent gene-based molecular signature classification 
#'by means of permutation test.
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
#'    \describe{
#'        \item{$class}{factor with all possible class levels.}   
#'    }
#'@slot permutation named list with at least the following fields:   
#'    \describe{
#'        \item{$pvalues}{numeric matrix with subjects in row and classes 
#'            in columns.}   
#'        \item{$fdr}{numeric matrix with False Discovery Rate correction  
#'            of p-values by row.} 
#'    }
#'
#'@section Superclasses:   
#'None declared.   
#'
#'@section Subclasses:   
#'    \describe{
#'        \item{PAM50}{Peruo et al. (2000 and 2010) breast cancer subtypes,
#'             i. e., Luminal A, Luminal B, Basal, Her2 or Normal-like  
#'             subtypes as implemented in genefu library (Haibe-Kains et    
#'             al. 2014).}  
#'    }
#'
#'@section Functions:   
#'MolecularPermutationClassifier S4 class includes the following functions:  
#'    \itemize{
#'        \item Integrity check:  
#'            \describe{
#'                \item{validity}{will check appropriate annotation 
#'                    data.frame minimal required columns, all named   
#'                    parameters and if exprs and annotation dimension  
#'                    matches.}
#'                \item{prototype}{just for an empty class with default  
#'                    values: nPerm=1e4L, where="fdr", pCutoff=0.01, 
#'                    corCutoff=0.1 and keep=FALSE}.  
#'            }
#'    \item Generics:   
#'        \describe{
#'            \item{\code{\link{show}},\code{\link{print}}}{basic
#'             class display wrappers.} 
#'            \item{\code{\link{summary}}}{classifier statistics.}   
#'        }
#'    \item Constructors (as this class is virtual see subclass'    
#'        'documentation).
#'        \describe{
#'            \item{\code{\link{setAs}}}{MAList to \code{\link{PAM50}}}  
#'            \item{\code{\link{as.PAM50}}}{wrapper for \code{\link{PAM50}}  
#'                setAs from MAList.}  
#'            \item{\code{\link{loadBCDataset}}}{wrapper to load  
#'                BreastCancerXX data (Class, exprs, annotation,
#'                clinical data).}   
#'        }
#'    \item Getters for the corresponding slots (\code{\link{parameters}},  
#'         \code{\link{exprs}}, \code{\link{annotation}},  
#'         \code{\link{targets}}, \code{\link{classification}}  
#'         and \code{\link{permutation}}).  
#'    \item Setters for the corresponding slots (\code{\link{parameters<-}},  
#'         \code{\link{annotation<-}} and \code{\link{targets<-}}). 
#'    \item Particular (virtual) functions: 
#'        \describe{
#'            \item{\code{\link{filtrate}}}{remove from the  
#'                \code{\link{exprs}} matrix subjects not required by the  
#'                classification algorithm.}   
#'            \item{\code{\link{classify}}}{generate subject classification  
#'                according to subclasses implementation (PAM50, etc.).}   
#'            \item{\code{\link{permutate}}}{obtain subject classification  
#'                based on the null correlation distribution by means 
#'                permutation simulation.}   
#'            \item{\code{\link{subtypes}}}{obtain the new classification 
#'                using permutation results.}  
#'            \item{\code{\link{subjectReport}}}{a friendly report for 
#'                Physician treatment decision support.} 
#'            \item{\code{\link{databaseReport}}}{a pdf with all 
#'                subjectReports, if a database is available.}   
#'        }
#'    }
#'
#'@seealso \code{\link{PAM50}} for a complete example,   
#'\code{\link{loadBCDataset}} to load BreastCancerXX dataset,
#'\code{\link{filtrate}}, \code{\link{classify}} and  
#'\code{\link{permutate}} to get corresponding Breast Cancer subtype.  
#'Getters/Setters for this class are \code{\link{parameters}},   
#'\code{\link{exprs}}, \code{\link{annotation}}, \code{\link{targets}},  
#'\code{\link{classification}} and \code{\link{permutation}}.  
#'
#'@include pbcmcPackage.R   
#'@docType class   
#'@importClassesFrom methods ANY classGeneratorFunction data.frame
#'@importClassesFrom methods list matrix 
#'@name MolecularPermutationClassifier-class   
#'@rdname MolecularPermutationClassifier-class   
#'@exportClass MolecularPermutationClassifier   
#'@family MolecularPermutationClassifier   
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez    
#'    \email{efernandez@@bdmg.com.ar}
#'@references    
#'    \enumerate{
#'        \item Haibe-Kains B, Schroeder M, Bontempi G, Sotiriou C and   
#'             Quackenbush J, 2014, genefu: Relevant Functions for Gene    
#'             Expression Analysis, Especially in Breast Cancer. R package
#'             version 1.16.0, \url{www.pmgenomics.ca/bhklab/} 
#'        \item Perou CM, Sorlie T, Eisen MB, et al., 2000, Molecular  
#'             portraits of human breast tumors. Nature 406:747-752 
#'        \item Perou CM, Parker JS, Prat A, Ellis MJ, Bernard PB., 2010, 
#'             Clinical implementation of the intrinsic subtypes of breast
#'             cancer, The Lancet Oncology 11(8):718-719   
#'    }
#'
setClass(
    Class="MolecularPermutationClassifier", contains="VIRTUAL",   
    slots=c(
        parameters="list",
        exprs="matrix",
        annotation="data.frame",
        targets="data.frame",
        classification="list",
        permutation="list"
    ),
    validity=function(object){
        ##Check the integrity  
        ##annotation slot   
        stopifnot(all(c("probe", "EntrezGene.ID", "NCBI.gene.symbol") %in% 
            names(annotation(object))))
        stopifnot(nrow(exprs(object))==nrow(annotation(object)))
        ##parameters
        stopifnot(all(c("nPerm", "where", "pCutoff", "keep") %in%    
            names(parameters(object))))
        ##exprs and annotation consistency 
        if(nrow(exprs(object))>0){
            stopifnot(all(row.names(annotation(object)) ==   
                row.names(exprs(object))))
        }
        return(TRUE)
    },
    prototype=list(
        parameters=list(nPerm=1e4L, where="fdr", pCutoff=0.01, keep=FALSE), 
        exprs=matrix(ncol=0,nrow=0),    
        annotation=data.frame(probe=character(), EntrezGene.ID=integer(),   
            NCBI.gene.symbol=character()),
        classification=list(),
        permutation=list()
    )
)
