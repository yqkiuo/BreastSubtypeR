#' @title Single-Sample Predictors for Breast Cancer (sspbc)
#'
#' @description
#' This function assign classes to breast cancer samples using a selection of provided models. It works on raw gene expression data.
#' Provided models were developed using gene expression data from mRNAseq generated using HiSat/StringTie.
#' For recommended usage see details and examples.
#'
#' @usage applySSP(tsv, ssp="",
#' plot=FALSE,txt=FALSE, report=FALSE, add.is.num=TRUE,
#' mylas=1,ssp.name="",gex=NULL, id=NULL, id.type="Gene.ID",
#' full.out=FALSE, output)
#'
#' @param tsv {a character string specifying file path to a gene.tsv file from StingTie that must contain column "Gene.ID" and column "FPKM".
#'		Gene.ID is expected to be GENCODE Human Release 27 gene identifiers.}
#' @param ssp {a character string specifying file name of a Rdata file equal to a ssp model (R object named aims.gs).}
#' @param plot {logical with default FALSE. Should a pdf file (all.probs.pdf) with plot of posterior probabilities be created?}
#' @param txt {logical with default FALSE. Should txt files (class.txt, all.probs.txt) with class and posterior probabilities be created?}
#' @param report {logical with default FALSE. Should reporting info be printed when using function?}
#' @param add.is.num {logical with default FALSE. Should $isnumeric be set? Option to fix a practical issue with $isnumeric variable that was not set for original models.}
#' @param mylas {numeric with default 1. Used in plot of posterior probabilities  (see plot)}
#'
#' @param ssp.name {One of "" (default), "ssp.er", "ssp.pr", "ssp.her2", "ssp.her2.ern", "ssp.her2.erp", "ssp.ki67", "ssp.nhg", "ssp.pam50", "ssp.subtype", "ssp.ror", "ssp.cc15".
#' Option to specify ssp model by short name (ssp models provided as part of package, sspbc.models list with models).
#'  If argument ssp.name is specified, i.e., not "" (default), argument ssp is ignored and sspbc.models list is used to load ssp model.}
#'
#' @param gex {NULL (default) or a vector or matrix with numeric gene expression data. Option to input gex data as vector or matrix.
#'		If gex is a matrix rows are expected to be genes and columns are expected to be samples.}
#'
#' @param id {NULL (default) or a character vector corresponding to ID for genes for gene expression data provided as vector or matrix.
#'		Expected default type for gene identifiers is Gene.ID (see id.type)
#'		If gex is provided, i.e., not NULL (default), argument tsv is ignored and argument id must be provided.
#'		If gex is provided, i.e., not NULL (default), the option to save results as txt and plot are disabled. Arguments 'plot', 'txt', and 'mylas' are ignored.}
#'
#' @param id.type {One of "Gene.ID" (default), "Gene.Name","HGNC","EntrezGene". Specify input gene identifier.
#'		Used by function translate_id2entrez_v1.2.R to translate gene identifiers before classification with provided ssp models.}
#'
#' @param full.out {logical with default FALSE. Should full classification resultslist be returned?
#'			Option to get output as resultslist (as by applyAIMS in the AIMS package)
#'			If full.out=FALSE only class will be returned.
#'			If full.out=TRUE the complete resultslist from applyAIMS will be returned.}
#'
#' @param output {a character string specifying file path to an output directory for txt files and plot.}
#'
#' @details
#' applySSP works by assigning classes by simple gene rules on the the form gene A < gene B defined by different
#' models as described by Paquet and Hallett (2014) J Natl Cancer Inst 107(1):357. The development of the models included
#' with applySSP is described in Staaf J. et al. medRxiv 2021.12.03.21267116; https://doi.org/10.1101/2021.12.03.21267116
#'
#' For most users the recommended way of using applySSP is by providing gene expression data as a loaded numeric matrix or vector
#' and by specifying ssp model by short name and by providing gene identifiers as a character vector (see examples).
#'
#' For other users the possibility to provide gene expression data in a gene.tsv file from StingTie and to specify
#' ssp model by full name from development is available.
#'
#' The included ssp models were developed on gene expression data using gene definitions from GENCODE Human Release 27
#' summarized on GENCODE gene and including all genes annotated with EntrezGene. When providing gene identifiers as a character
#' vector, applySSP can handle a number of different types of gene identifiers (see id.type) and will translate these to EntrezGene
#' as needed.
#'
#'
#'@return
#' \item{cl}{Molecular class identified by a given ssp model.}
#'
#' \item{all.probs}{Posterior probability values for all samples and all subtypes.}
#'
#'
#'@author
#' Johan Staaf (johan.staaf@@med.lu.se),
#' Johan Vallon-Christersson (johan.vallon-christersson@@med.lu.se)
#'
#' @seealso
#' \code{\link{sspbc.models}
#' \link{sspbc.models.fullname}
#' \link{Gene.ID.ann}
#' \link{testmatrix}
#' }
#'
#' @export
#'
#' @import e1071
#' @import methods
#' @import stats
#' @import grid
#' @noRd
#'
#' @references
#' Staaf J. et al. medRxiv 2021.12.03.21267116
#' (\href{https://doi.org/10.1101/2021.12.03.21267116}{medRxiv})
#' (\href{https://github.com/StaafLab/sspbc}{GitHub})
#'
#' Paquet and Hallett (2014) J Natl Cancer Inst 107(1):357
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/25479802}{PubMed})
#' (\href{https://github.com/meoyo/AIMS}{GitHub})
#'
#' @examples
#' ## Load the example dataset
#' data(testmatrix)
#'
#' ## Load the sspbc models
#' data(sspbc.models)
#' data(sspbc.models.fullname)
#'
#' ## Assign status to multiple samples in the testmatrix for the selected ssp model (recommended for most users)
#' myresults <- applySSP(gex=testmatrix, id=row.names(testmatrix), ssp.name="ssp.pam50")
#' myresults <- applySSP(gex=testmatrix, id=row.names(testmatrix),
#' ssp="Training_Run19081Genes_noNorm_SSP.subtypeMost.Fcc15_5x5foldCV.num.rules.50_24.selRules.AIMS.GS.RData")
#'

applySSP <- function(tsv, ssp="", plot=FALSE, txt=FALSE, report=FALSE, add.is.num=TRUE, mylas=1,
                     ssp.name="", gex=NULL, id=NULL, id.type="Gene.ID", full.out=FALSE,
                     output){
  
  data("Gene.ID.ann")
  
  if(report){
    # report on Gene.ID.ann
    cat("Gene.ID.ann nrow:",nrow(Gene.ID.ann),"\n")
  }
  
  # load gex data as matrix
  # if gex=NULL read data from tsv file
  if(is.null(gex)){
    # read tsv from single sample gene.tsv file
    mymatrix <- read_StringTie_tsv_FPKM(tsv, id=Gene.ID.ann$Gene.ID)
    # head(mymatrix)
  }else{
    # if data is a matrix make sure it is numeric
    if(is.matrix(gex)){
      if(is.numeric(gex)){
        mymatrix <- gex
        rownames(mymatrix) <- id
      }else{
        stop("provided gex matrix is not numeric")
      }
    }else{
      # if data is a vector convert it into matrix
      if( is( gex, "numeric")){
        mymatrix <- as.matrix(gex)
        rownames(mymatrix) <- id
      }
    }
  }
  
  # translate gene id
  myid <- translate_id2entrez(id=rownames(mymatrix), ann=Gene.ID.ann, id.type=id.type, e=TRUE)
  
  # load ssp model
  if((ssp.name=="") & (ssp %in% names(sspbc.models.fullname))){
    aims.gs <- sspbc.models.fullname[[ssp]] # loaded ssp object must be named aims.gs
  }else{
    # ssp.name <- "ssp.cc15"
    if(ssp.name %in% names(sspbc.models)){
      aims.gs <- sspbc.models[[ssp.name]]
    }else{
      stop(paste("Specified ssp.name not among names in sspbc.models: ", paste(names(sspbc.models), collapse=", ")))
    }
  } #
  
  # fix issue with $isnumeric
  if(add.is.num){
    # add logi (aims.gs$ all.pairs), i.e, one TRUE for each aims.gs$ all.pairs
    aims.gs$ one.vs.all.tsp[[aims.gs $ k]]$isnumeric[names(aims.gs$ one.vs.all.tsp[[aims.gs $ k]]$tables)] <- TRUE
  }#
  
  # apply ssp model
  resultslist <- applyAIMS(mymatrix, myid, aims.gs)
  
  names(resultslist)
  
  # handle results when a single sample is classified
  if(ncol(mymatrix)==1){
    cl <- resultslist$cl[[1]]
    all.probs <- c(resultslist$all.probs[[1]])
    names(all.probs) <- colnames(resultslist$all.probs[[1]])
    
    if(report){
      cat("assigned class:",cl,"\n")
    }
    
    # plot only available when a single sample is classified
    if(plot){
      pdf(file=paste(output,"all.probs.pdf",sep="/"), width = 2.53, height = 1.19, pointsize = 6)
      op <- par(oma=c(2,2,0,0), mar=c(2,2,1,0))
      barplot(all.probs, las=mylas)
      mtext("posterior probability", side=2, line=3, cex=0.9)
      # mtext("class", side=1, line=2, cex=0.9)
      par(op)
      dev.off()
    }# end if plot
    
    # text only available when a single sample is classified
    if(txt){
      write.table(data.frame(all.probs), file=paste(output,"all.probs.txt",sep="/"), col.names=FALSE, quote=FALSE, sep="\t")
      write.table(cl, file=paste(output,"class.txt",sep="/"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    }# end if txt
    
  }
  
  # handle results when a matrix with multiple samples are classified
  if(ncol(mymatrix)>1){
    cl <- resultslist$cl
    colnames(cl) <- "class"
    all.probs <- resultslist$all.probs[[1]]
    row.names(all.probs) <- row.names(cl)
  }
  
  if(full.out){
    return(resultslist)
  }else{
    return(cl)
  }
  
}# end function




###########################################
#	Johan Vallon-Christersson             #
#	johan.vallon-christersson@med.lu.se   #
###########################################


####
# function: read_StringTie_tsv_FPKM
#
#	INPUT:
#		- gene.tsv file from StingTie. Must contain column holding gene id with column name 'Gene ID' or 'Gene.ID' and column 'FPKM' with geneexpression data.
#		- vector of unique gene id to return data for (genes to be included in returned matrix). Must be same type of id found in column 'Gene.ID' in gene.tsv file.
#	OUTPUT:
#		- returns genematrix with one column (FPKM data sum on gene id) and rows for gene id (row names are id)

#	v1:	first implementation
# 		the gene.tsv is read using read.delim() and invalid characters in column names are translated to ".". That is, column name 'Gene ID' will be translated to 'Gene.ID'.



#	# examples and manual stuff
#			
#		# specify input gene.tsv 
#			gene.tsv <- "gene.tsv"
#				
#		# load gene annotation
#			load("Gene.ID.ann.Rdata")
#				
#			# some gene id to collect
#			some.gene.id <-  Gene.ID.ann$Gene.ID[1:5]
#			
#		# run function
#			mymatrix <- read_StringTie_tsv_FPKM(tsv=gene.tsv, id=some.gene.id, report=TRUE)
#
#				head(mymatrix)


#####################################
# function read_StringTie_tsv_FPKM
#####################################

# function
read_StringTie_tsv_FPKM <- function(tsv, id, report=FALSE){
  
  # number of columns and rows in genematrix 
  my.ncol <- 1
  my.nrow <- length(id)
  
  # test verify that id are unique  
  if(length(id)==length(unique(id))){
    if(report){
      cat("provided id(s) are unique, ")	
    }
  }else{
    stop("mfkr id")
  }			
  
  # create genematrix
  genematrix <- matrix(nrow=my.nrow, ncol=my.ncol)
  
  # rownames and colname
  row.names(genematrix) <- id
  colnames(genematrix) <- "FPKMsum"
  
  # head(genematrix)
  
  # read gene.tsv
  if(report){
    cat("creating matrix, ")
  }
  
  # read from file
  read.data <- read.delim(file=tsv, as.is=TRUE)
  # head(read.data)
  
  # aggregate fpkm on gene 	
  sum_fpkm_gene <- aggregate(FPKM ~ Gene.ID, data=read.data, sum)
  # head(sum_fpkm_gene)
  
  # find id index in sum_fpkm_gene$Gene.ID
  id.i <- match(id, sum_fpkm_gene$Gene.ID)
  
  # test verify that all provided id are found
  if(length(which(is.na(id.i)==TRUE))==0){
    # add id.i FPKM to genematrix
    genematrix[,1] <- sum_fpkm_gene$FPKM[id.i]
  }else{
    stop("provided id not found in tsv")
  }
  
  if(report){
    cat("done")	
  }
  return(genematrix)
  
}# end function

#####################################
#####################################


###########################################
#	Johan Vallon-Christersson             #
#	johan.vallon-christersson@med.lu.se   #
###########################################


####
# function: translate_id2entrez
#
#	INPUT:
#		- vector of gene id to translate (id of type used in our StringTie pipeline, i.e., 'Gene.ID')
#			- for v1.2, other identifiers supported. id is expected to be inputted as.character or will be transformed as.caharcter
#		- gene annotation object uncluding columns 'Gene.ID' and 'EntrezGene' (the gene annotation file created for our StringTie pipeline).
#			- for v1.2, to use other identifier as input these must also be included in the gene annotation object. 
#		- vector of unique gene id to return data for (genes to be included in returned matrix). Must be same type of id found in column 'Gene.ID' in gene.tsv file.
#	OUTPUT:
#		- vector of EntrezGene or vector with EntrezGeneMust with appended prefix 'e' (as required by our SSP models).
#
#	v1:	first implementation
#		Returned vector will have NA for id not found in gene annotation object.
#		Note that in gene annotation object all id (Gene.ID) without an EntrezGene have character "NA" in column EntrezGene and these will be returned
#
#	v1.2:
#		added parameter 'id.type' (default Gene.ID) to specify type for inputted ID
#			allowed vales: Gene.ID (default), Gene.Name, HGNC, or EntrezGene.
#		
#			
#		
#		




#	# examples and manual stuff
#				
#		# load gene annotation
#			load("Gene.ID.ann.Rdata")
#				
#			# some gene id to translate examples
#			some.gene.id <-  Gene.ID.ann$Gene.ID[c(1:5)]
#				# some.gene.id <-  Gene.ID.ann$Gene.ID[c(1:5, 19625)]
#				# some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
#			# type <- "HGNC"
#			some.gene.id <-  Gene.ID.ann$Gene.Name[c(1:5, 19625)]
#				# some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
#			# type <- "EntrezGene"
#			some.gene.id <-  Gene.ID.ann$EntrezGene[c(1:5, 19625)]
#				# some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
#			
#		# run function
#			myid <- translate_id2entrez(id=some.gene.id, ann=Gene.ID.ann, e=TRUE, report=TRUE)
#
#				head(myid)


#####################################
# function translate_id2entrez
#####################################

# function
translate_id2entrez <- function(id, ann, id.type="Gene.ID", e=FALSE, report=FALSE){
  
  if(id.type %in% c("Gene.ID","Gene.Name","HGNC","EntrezGene")){
    # find first index for id in ann
    id.i <- match(as.character(id), ann[,id.type])
  }else{
    stop("specified type not supported")
  }		
  
  if(report){
    # test verify that all provided id are found
    if(length(which(is.na(id.i)==TRUE))==0){
      cat("all id found in ann, ")				
    }else{
      cat(length(which(is.na(id.i)==TRUE)), "id not found in ann, ")
    }
  }
  
  # get EntrezGene and assign id as names 
  entrez <- ann$EntrezGene[id.i]
  names(entrez) <- id
  
  entrezNA <- length(which(entrez[is.na(entrez)==FALSE]=="NA"))
  
  if(report){
    if(entrezNA>0){
      cat(entrezNA, "found id have EntrezGene NA, ")
    }else{
      cat("all found id have an EntrezGene, ")
    }				
  }
  
  # add prefix	
  if(e){
    if(report){
      cat("adding prefix e, ")				
    }
    entrez[is.na(entrez)==FALSE & entrez!="NA"] <- paste("e", entrez[is.na(entrez)==FALSE & entrez!="NA"], sep="")
  }				
  
  if(report){
    cat("done")	
  }
  
  return(entrez)
  
}# end function

#####################################
#####################################




#### SSP_functions.R ####
.smaller <- function(x,y){
  x < y
}

## Functions necessary to run AIMS
.comp.sel.pairs <- function(dataset,sel.pairs,func=.smaller){
  to.ret <- matrix(NA,nrow=length(sel.pairs),ncol(dataset$exprs))
  for (spi in 1:length(sel.pairs)){
    ss.cur <- strsplit(sel.pairs[spi],"<")[[1]]
    gene1 = which(dataset$GeneName == ss.cur[1])
    gene2 = which(dataset$GeneName == ss.cur[2])
    #stopifnot(length(gene1) == 1 & length(gene2) == 1)
    if (length(gene1) == 1 & length(gene2) == 1){
      to.ret[spi,] <- func(dataset$exprs[gene1,],dataset$exprs[gene2,])
    }
    else{
      message(paste("You are missing the pair or have more than one",sel.pairs[spi],"in",dataset$name))
    }
  }
  to.ret <- apply(to.ret,2,as.numeric)
  rownames(to.ret) <- sel.pairs
  to.ret
}

.one.vs.all.tsp <- function(D,GeneName,one.vs.all.tsp){
  ## Need to add some QC
  ## First compute
  train.pairs <- .comp.sel.pairs(list(exprs=D,GeneName=GeneName),one.vs.all.tsp$all.pairs)

  classes <- matrix("",ncol=length(one.vs.all.tsp$k),nrow=ncol(D))
  prob <- matrix(0,ncol=length(one.vs.all.tsp$k),nrow=ncol(D))
  colnames(classes) <- colnames(prob) <- as.character(one.vs.all.tsp$k)
  rownames(classes) <- rownames(prob) <-colnames(D)
  nb.d <- data.frame(t(train.pairs))
  all.probs <- list()
  for (ki in one.vs.all.tsp$k){
    message(sprintf("Current k = %d",ki))
    prob.train <- predict(one.vs.all.tsp$one.vs.all.tsp[[ki]],nb.d,type="raw")
    cur.cl <- apply(prob.train,1,function(prob.cur){colnames(prob.train)[which.max(prob.cur)]})
    cur.prob <- apply(prob.train,1,function(prob.cur){max(prob.cur)})
    prob[,as.character(ki)] <- cur.prob
    classes[,as.character(ki)] <- cur.cl
    all.probs[[as.character(ki)]] <- prob.train
  }
  invisible(list(cl = classes,prob = prob,all.probs = all.probs,rules.matrix=train.pairs))
}

.get.all.pairs.genes <- function(all.pairs){
  genes <- c()
  for (cp in strsplit(all.pairs,"<")){
    genes <- c(genes,cp)
  }
  unique(genes)
}
## Remove the duplicated Entrez by keeping the most higly expressed
## This is closer to the single sample selection
## D is a raw gene expression matrix rows == genes and columns patients
.removeDuplicatedEntrezPerPatients <- function(D,EntrezID,probes){
  ## Maybe we have nothing to do already
  if (all(!duplicated(EntrezID))){
    return(list(dataset=D,EntrezID=EntrezID))
  }
  else{
    uniqEntrez <- sort(unique(EntrezID))
    newD <- matrix(0,nrow=length(uniqEntrez),ncol=ncol(D))
    if (!missing(probes)){
      sel.probes <- matrix("",nrow=length(uniqEntrez),ncol=ncol(D))
    }
    for (i in 1:ncol(D)){
      curD <- D[,i]
      curEntrez <- EntrezID
      oi <- order(curD,decreasing=TRUE) ## order by raw expression
      curD <- curD[oi]
      curEntrez <- curEntrez[oi]
      cur.sel <- !duplicated(curEntrez) ## remove duplicated
      curD <- curD[cur.sel]
      curEntrez <- curEntrez[cur.sel]
      reorder <- match(uniqEntrez,curEntrez)
      newD[,i] <- curD[reorder]
      if (!missing(probes)){
        sel.probes[,i] <- probes[oi][cur.sel][reorder]
      }
    }
    colnames(newD) <- colnames(D)
    if (!missing(probes)){
      colnames(sel.probes) <- colnames(D)
      return(list(dataset=newD,EntrezID=uniqEntrez,probes=sel.probes))
    }
    else{
      return(list(dataset=newD,EntrezID=uniqEntrez))
    }
  }
}

.apply.nbc <- function(D,EntrezID,sel.nbc){
  ## Verify the number of rows of D and EntrezIDs have the same length 
  if (nrow(D) != length(EntrezID)){
    stop(sprintf("You need the same number of rows and EntrezID. Right now nrow(D) = %d and length(EntrezID) = %d",nrow(D),length(EntrezID)))
  }
  ## AIMS needs to be applied on expression values > 0. Require 95% of the values to be > 0
  if (!all(apply(D,2,function(x){(sum(x < 0,na.rm=TRUE)/length(x)) < 0.05}))){
    stop("AIMS needs to be applied on expressionn values > 0. Did you gene-centered your matrix D? You should not.")
  }
  
  ## Verify if D is a numeric matrix
  if (!all(apply(D,2,is.numeric))){
    stop("Verify D is a numeric matrix. apply(D,2,as.numeric) could do the job or verify the first column doesn't contain probe ids") 
  }
  
  D <- apply(D,2,as.numeric)
  EntrezID <- as.character(EntrezID)
  sel.ids.nb <- .get.all.pairs.genes(sel.nbc$all.pairs)
  sel.gn <- EntrezID %in% sel.ids.nb
  D <- D[sel.gn,,drop=FALSE]
  Entrez <- EntrezID[sel.gn]
  col.D.test <- .removeDuplicatedEntrezPerPatients(D,Entrez)
  pred.test <- .one.vs.all.tsp(D=col.D.test$dataset,
                              GeneName=col.D.test$EntrezID,
                              sel.nbc)

  ## Add more information to the output variable
  pred.test$data.used <- col.D.test$dataset
  pred.test$EntrezID.used <- col.D.test$EntrezID
  
  invisible(pred.test)
}

applyAIMS <- function(eset,EntrezID,AIMSmodel){
  D <- NA
  if (!is(eset,"ExpressionSet") & !is(eset,"matrix")){
    stop("eset argument should be either an ExpressionSet (Biobase) or a numerical matrix")
  }
  
  if (is(eset,"ExpressionSet")){
    D <- exprs(eset)
  }
  else if (is(eset,"matrix")){
    D <- eset
  }

  .apply.nbc(D,EntrezID,AIMSmodel)
}

