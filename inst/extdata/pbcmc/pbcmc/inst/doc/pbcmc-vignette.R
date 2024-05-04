### R code from vignette source 'pbcmc-vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Speedup
###################################################
texec<-matrix(ncol=2, nrow=7, byrow=TRUE,  
c(7,3.35,
    6,3.07,
    5,2.97,
    4,2.79,
    3,2.71,
    2,1.89,
    1,1.00))

colnames(texec)<-c("Cores", "SpeedUp")   
texec<-as.data.frame(texec)

library("ggplot2")
p<-ggplot(data=texec, aes(x=Cores, y=SpeedUp))+  
        geom_point(size=2)+geom_abline(slope=1,intercept=0)+
        # geom_smooth(se=FALSE)+   
        geom_line(color="blue", linetype="dashed")+   
        ylab("Speed Up")+xlab("Number of cores")+ 
        theme_bw()
p            


###################################################
### code chunk number 2: General options for R
###################################################
options(prompt="R> ", continue="+  ", width=70, useFancyQuotes=FALSE, digits=4) 
suppressMessages(library("pbcmc"))


###################################################
### code chunk number 3: Loading datasets
###################################################
library("pbcmc")
library("BiocParallel")
object<-loadBCDataset(Class=PAM50, libname="nki", verbose=TRUE)  
object


###################################################
### code chunk number 4: PAM50 with microarray data
###################################################
library("breastCancerNKI")
data("nki")


###################################################
### code chunk number 5: PAM50 with microarray data2
###################################################
##The expression   
M<-exprs(nki)[, 1:5, drop=FALSE]  
head(M)

##The annotation   
genes<-fData(nki)[, c("probe", "NCBI.gene.symbol", "EntrezGene.ID")] 
head(genes)

##Additional information (optional)  
targets<-pData(nki)[1:5, ,drop=FALSE]   
head(targets)


###################################################
### code chunk number 6: PAM50 with own data
###################################################
M<-pam50$centroids
genes<-pam50$centroids.map
names(genes)<-c("probe", "NCBI.gene.symbol", "EntrezGene.ID")  
object<-PAM50(exprs=M, annotation=genes)   
object


###################################################
### code chunk number 7: Exploring the slots
###################################################
head(exprs(object))      ##The gene expression values for each subject


###################################################
### code chunk number 8: Exploring the slots2
###################################################
head(annotation(object)) ##The compulsory annotation fields
head(targets(object))    ##The clinical data, if available.


###################################################
### code chunk number 9: Workflow
###################################################
object<-filtrate(object, verbose=TRUE)   
object<-classify(object, std="none", verbose=TRUE)  
object<-permutate(object, nPerm=10000, pCutoff=0.01, where="fdr", 
    corCutoff=0.1, keep=TRUE, seed=1234567890, verbose=TRUE, 
    BPPARAM=bpparam())


###################################################
### code chunk number 10: Permutation results
###################################################
object


###################################################
### code chunk number 11: Summary
###################################################
summary(object)


###################################################
### code chunk number 12: SubjectReport
###################################################
subjectReport(object, subject=1)   


###################################################
### code chunk number 13: Session Info
###################################################
sessionInfo()


