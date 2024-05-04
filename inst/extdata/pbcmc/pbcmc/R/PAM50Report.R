#'PAM50 permutation test results reports
#'
#'\code{subjectReport} is basically a grid.arrange object which basically 
#'consists of three main parts: a summary table, a two row ggplot2 
#'facet_wrap with scatter ggplots (Wickham 2009) of subject expression and   
#'PAM50 centroids (Perou et al. 2000 & 2010) and a textGrob with the
#'simulation parameter used. Particularly: 
#'\describe{
#'    \item{tableGrob}{ with the following fields:
#'        \describe{
#'            \item{$Summary}{subject name, PAM50 and Permuted subtype}   
#'            \item{$Fields}{for the five PAM50 subtypes:    
#'                \itemize{
#'                    \item Correlation: PAM50 centroid correlation with   
#'                        observed subject exprs.  
#'                    \item p-value: permutation p-value obtained using   
#'                        the simulation.   
#'                    \item FDR: adjusted p-value using False   
#'                         Discovery Rate.  
#'                }
#'            }
#'        }
#'    }
#'    \item{ggplot facet_wrap}{two rows to display scatter subject exprs vs
#'         PAM50 centroids, in addition to a the linear regression fix. If 
#'         subject, has an unique subtype, then the graph is in red. In
#'         addition, if simulated permutations were run with keep=TRUE
#'         option, then null distribution boxplots are plotted with
#'         observed correlations as a big round point.}     
#'    \item{textGrob}{the permutation @@parameter slot used in the  
#'         simulation.}   
#'}
#'
#'\code{summary} it basically prints descriptive data of PAM50 dataset,
#'the test parameters used, a frequency table of PAM50 Subtypes and a 
#'contingency table with Classes vs PAM50 Subtypes.  
#'
#'\code{databaseReport} basically is a pdf report where the first page is a 
#'global summary of the database, i.e., a \code{summary} contingency table   
#'of permutation test classes against original PAM50 subtypes results. The   
#'following pages are the database respective \code{subjectReport} outputs. 
#'
#'@param object a PAM50 object.
#'@param subject integer to select the appropriate subject to report.   
#'@param fileName character with the name of the pdf report file to save.
#'@param ... additional parameters for pdf function call. 
#'@param verbose should the user feedback be displayed? By default value is 
#'    "verbose" global option parameter, if present, or FALSE otherwise.
#'
#'@return depending on function call:
#'\item{subjectReport}{a grid.arrange object.}  
#'\item{databaseReport}{a pdf file with database summary and  
#'    \code{subjectReport}s.}
#'\item{summary}{Console summary statistics plus a data.frame}   
#'
#'@seealso \code{\link{PAM50}} for a complete example.   
#'
#'@include PAM50Class.R   
#'@importFrom ggplot2 aes element_blank element_text facet_wrap   
#'geom_blank geom_boxplot geom_point geom_text ggplot scale_colour_discrete   
#'scale_colour_manual stat_smooth theme xlab xlim ylab ylim  
#'@importFrom gridExtra tableGrob grid.arrange 
#'@importFrom reshape2 melt  
#'@importFrom ggplot2 geom_blank xlim ylim element_blank element_text geom_text 
#'@importFrom cowplot plot_grid  
#'@importFrom utils capture.output data setTxtProgressBar txtProgressBar   
#'@importFrom grDevices dev.off pdf 
#'@importFrom grid unit viewport 
#'@rdname PAM50SubjectReport   
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
#'         Clinical implementation of the intrinsic subtypes of breast
#'         cancer, The Lancet Oncology 11(8):718-719.   
#'    \item Wickham H, ggplot2: elegant graphics for data analysis.
#'          Springer New York, 2009.   
#'}
#'@examples
#'##Using pam50centroids package example data
#'data(pam50centroids)
#'pam50centroids
#'
#'##This object has already run filtrate, classify and permutate. So, now  
#'##we can obtain some reports:
#'##1) database summary  
#'summary(pam50centroids)
#'
#'##2)Individual subject report. If keep=FALSE boxplot panel is not available   
#'subjectReport(pam50centroids, subject=1)##Basal subtype  
#'subjectReport(pam50centroids, subject=1)##Her2 subtype  
#'
#'##3) complete database report 
#'#databaseReport(pam50centroids, fileName="PAM50.pdf", verbose=TRUE)  

setMethod(f="subjectReport", signature="PAM50", definition=function(object,  
    subject){
    ##Check data   
    stopifnot(!missing(subject))
    stopifnot(length(subject) == 1)  
    stopifnot(subject <= ncol(exprs(object)))  

    ##Subject summary   
    subjectSummary<-c(
        Subject=colnames(exprs(object))[subject],
        PAM50.Subtype=as.character(classification(object)$subtype[subject]),
        Permuted.Subtype=
        switch(as.character(permutation(object)$subtype$Permuted[subject]),
            "Assigned"=permutation(object)$subtype$Subtype[subject],
            "Not Assigned"="Not Assigned",  
            "Ambiguous"=paste("Ambiguous:\n", gsub(pattern=", ",  
                replacement=",\n",    
            permutation(object)$subtype$Classes[subject]), sep="")   
        )
    )
    subjectSummary<-as.data.frame(subjectSummary)
    names(subjectSummary)<-"Summary"
    subjectSummary$Field<-c("Correlation", "p-value", "FDR")  
    subjectSummary<-cbind(subjectSummary, rbind(   
        classification(object)$correlation[subject,],    
        permutation(object)$pvalue[subject,],
        permutation(object)$fdr[subject,]))

    ##Match subject centroids genes with PAM50 ones  
    ##PAM50 exprs   
    pam50exprs<-genefu::pam50$centroids
    row.names(pam50exprs)<-genefu::pam50$centroids.map$EntrezGene
    pam50exprs<-pam50exprs[row.names(pam50exprs) %in%   
        annotation(object)$EntrezGene.ID,]
    pam50exprs<-pam50exprs[order(row.names(pam50exprs)), ]   

    ##Subject exprs   
    subjectexprs<-exprs(object)[, subject]   
    names(subjectexprs)<-annotation(object)$EntrezGene.ID
    subjectexprs<-subjectexprs[row.names(pam50exprs)]

    ##Data reshape   
    pam50exprs<-cbind(pam50exprs, Subject=subjectexprs)   
    pam50exprs<-as.data.frame(pam50exprs)
    pam50exprs<-na.omit(pam50exprs)
    pam50exprs<-melt(pam50exprs, id.vars="Subject")   
    names(pam50exprs)<-c("Subject", "Subtype", "exprs")  

    ##Panel 6x6 subject exprs vs PAM50 centroids  
    ##golbal variables   
    Subject<-Subtype<-Color<-Rho<-NULL    
    p<-ggplot()
    p<-p+geom_point(data=pam50exprs, aes(x=exprs, y=Subject))  
    p<-p+facet_wrap(~Subtype, nrow=2)   
    p<-p+stat_smooth(data=pam50exprs, aes(x=exprs, y=Subject),  
        method="lm", se=FALSE)   
    p<-p+xlab("PAM50 Centroids") + ylab("Subject Expression")
    ##Change point colour if Assigned
    if(permutation(object)$subtype$Permuted[subject]=="Assigned"){
        p<-p+geom_point(aes(x=exprs, y=Subject, colour="red"),  
        data=subset(pam50exprs,    
            pam50exprs$Subtype==classification(object)$subtype[subject]))    
        p<-p+scale_colour_discrete(guide=FALSE)
    }

    ##Check for permuted.cor  
    if("correlation" %in% names(permutation(object))){  
        ##Format the data  
        ##Observed
        correlation<-melt(classification(object)$correlation[subject, ],   
            value.name="Rho")
        correlation$Subtype<-row.names(correlation)
        correlation$Color<-rep("black", nrow(correlation))   
        ##Permuted
        pCorrelations<-melt(permutation(object)$correlation[[subject]])
        names(pCorrelations)<-c("n.perm", "Subtype", "Rho")  
        pCorrelations$Color<-rep("black", nrow(pCorrelations))   
        ##Color the Assigned subtype, if available   
        if(permutation(object)$subtype$Permuted[subject]=="Assigned"){
            ##Observed
            id<-correlation$Subtype ==   
                permutation(object)$subtype$Subtype[subject]
            correlation$Color[id]<-"red"
            ##Permuted
            id<-pCorrelations$Subtype ==   
                permutation(object)$subtype$Subtype[subject]
            pCorrelations$Color[id]<-"red"    
        }
        correlation$Color<-factor(correlation$Color)
        pCorrelations$Color<-factor(pCorrelations$Color,    
            levels=levels(correlation$Color))

        ##Make the plot  
        p2<-ggplot()
        p2<-p2+geom_boxplot(data=pCorrelations, aes(y=Rho, x=Subtype,  
            color=Color))
        p2<-p2+geom_point(data=correlation, aes(y=Rho, x=Subtype,  
            color=Color), size=5)   
        p2<-p2+scale_colour_manual(values=levels(correlation$Color),    
            guide=FALSE)
        p2<-p2+theme(axis.text.x=element_text(angle=90))    
    }

    ##The report output  
    pblank<-ggplot()+geom_blank()
    pblank<-pblank+xlim(-2,2)+ylim(0,2)
    pblank<-pblank+xlab("")+ylab("")
    pblank<-pblank+theme(axis.ticks=element_blank(),    
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        strip.text.x = element_blank(),  
        axis.line.x=element_blank(), axis.line.y=element_blank(),   
        plot.margin = unit(c(0,0,0,0), "cm")) 
    t1<-tableGrob(format(subjectSummary, digits=3))   
    t2<-pblank+geom_text(aes(x=0, y=1, label=paste("Permutations:",  
        parameters(object)$nPerm, ", pcutoff<",  
        parameters(object)$pCutoff, " & corCutoff>", 
        parameters(object)$corCutoff, sep="")))   
    ##Check for permuted.cor  
    if("correlation" %in% names(permutation(object))){  
#         plot_grid(plotlist=list(t1, p,  
#             t2+theme(plot.margin = unit(c(0, 0,  0, 0), "lines"))),    
#             ncol=1, rel_heights=c(1, 2, 0.3))
        grid.arrange(t1, p, t2, ncol=1, nrow=3, heights=c(1, 2, 0.2))     
        ##fill the 6x6 hole with the permutation boxplot 
        vp<-viewport(width=1/3-0.009, height=0.3, x=1-1/3+0.009, y=0.25, 
            just="left")
        print(p2, vp = vp) 
#         vp<-viewport(width=1/3-0.009, height=1/2.5,  
#             x=1-1/3+0.009, y=0.2, just="left") 
#         return(print(p2, vp = vp))
    }else{
        grid.arrange(t1, p, t2, ncol=1, nrow=3, heights=c(1, 2, 0.2)) 
#         return(plot_grid(plotlist=list(t1, p,  
#             t2+theme(plot.margin = unit(c(0, 0,  0, 0), "lines"))),    
#             ncol=1, rel_heights=c(1, 2, 0.3)))
        
    }
})
#'
#'@inheritParams subjectReport   
#'@rdname PAM50SubjectReport   
#'@family PAM50   
setMethod(f="databaseReport", signature="PAM50", definition=function(object,  
    fileName, ..., verbose=getOption("verbose",default=TRUE)){  
    ##Check data   
    stopifnot(!missing(fileName))

    ##The report output  
    databaseSummary<-suppressMessages(summary(object))
    pblank<-ggplot()+geom_blank()
    pblank<-pblank+xlim(-2,2)+ylim(0,2)
    pblank<-pblank+xlab("")+ylab("")
    pblank<-pblank+theme(axis.ticks=element_blank(),    
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        strip.text.x = element_blank(),  
        axis.line.x=element_blank(), axis.line.y=element_blank(),   
        plot.margin = unit(c(0,0,0,0), "cm")) 
    t1<-pblank+geom_text(aes(x=0, y=1, label=paste("Permutations:",  
        parameters(object)$nPerm, ", p-value/FDR cutoff <",    
        parameters(object)$pCutoff, "& Correlation cutoff >",    
        parameters(object)$corCutoff, sep="")))   
    t2<-tableGrob(databaseSummary)

    if(verbose){
        message("Generating database PAM50 permutation report..." )   
        pb <- txtProgressBar(min=0, max=ncol(exprs(object)),initial=0,style=3) 
    }
    ##File report generation  
    pdf(file=fileName, ...)   
    grid.arrange(t2, t1, ncol=1, nrow=2, heights=c(1, 0.1))   
    out<-lapply(1:ncol(exprs(object)), function(subject){   
        ##Permutation progress...   
        if(verbose){setTxtProgressBar(pb, subject)}   
        subjectReport(object=object, subject=subject)   
    })
    dev.off()
    ##Report completed feedback  
    if(verbose){
        setTxtProgressBar(pb, ncol(exprs(object)))   
        close(pb)
        message("Generating database PAM50 permutation report... done." )  
    }
})
#'
#'@inheritParams subjectReport   
#'@rdname PAM50SubjectReport   
#'@export summary   
#'@family PAM50   
setMethod(f="summary", signature="PAM50", definition=function(object, ...){ 
    ##Check if permutation test has been run  
    if(!"subtype" %in% names(permutation(object))){  
        stop("No PAM50 permutation results found!!!")
    }else{
        ##Summary head   
        message("PAM50 Permutation Test results!!! \n")
        message("Permutations: ", parameters(object)$nPerm)  
        message("pcutoff < ", parameters(object)$pCutoff) 
        message("corCutoff > ", parameters(object)$corCutoff, "\n")
        message("Global results:")   
        aux<-capture.output(show(table(permutation(object)$subtype$Permuted)))
        message(paste(aux[2], "\n", aux[3], sep="")) 
        message("Particular results:")   
        permSummary<-table(Classes=permutation(object)$subtype$Classes,    
        Subtype=permutation(object)$subtype$Subtype, useNA="always")   
        colnames(permSummary)[is.na(colnames(permSummary))]<-"Not Assigned"   
        rownames(permSummary)[is.na(rownames(permSummary))]<-"Not Assigned"   
        return(permSummary)
    }
})
