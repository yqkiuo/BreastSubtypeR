#'Example PAM50 objects for pbcmc package   
#'
#'The dataset corresponds to the Permutation-Based Confidence for Molecular
#'Classification package \code{\link{PAM50}} example objects, that was  
#'\code{\link{filtrate}}d, \code{\link{classify}}ed and  
#'\code{\link{permutate}}d using the following parameters:
#'\describe{
#'    \item{Permutations}{10000}
#'    \item{fdr}{0.01}
#'    \item{corCutoff}{0.1}
#'    \item{keep}{TRUE}
#'}
#'
#'@format pam50centroids corresponds with \strong{pam50$centroids} dataset   
#'available in genefu package. 
#'
#'@return a PAM50 object with the results obtained for pam50centroids   
#'simulations under the given parameters (see Detail section.) 
#'
#'@references    
#'\enumerate{
#'    \item Bioscience data mining group \url{http://www.bdmg.com.ar}   
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
#'@author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, German A. Gonzalez  
#'    \email{ggonzalez@@bdmg.com.ar}, Andrea S. Llera 
#'    \email{allera@@leloir.org.ar} and Elmer Andres Fernandez
#'    \email{efernandez@@bdmg.com.ar}
#'
#'@docType data   
#'@include pbcmcPackage.R   
#'@name pam50centroids   
#'@usage data(pam50centroids)   
#'@rdname pbcmcData   
#'@keywords datasets pam50 centroids 
#'@family PAM50   
NULL
