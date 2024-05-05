#' 
#' Functions adapted from original parker-based PAM50 subtyping
#' @name PCA50
#' @import gplots
#' @import RColorBrewer
#' @import lattice
#' @import genefilter
#' @import ctc
#' @import heatmap.plus
#' @import magrittr
#' @import dplyr
#' @import tidyverse
NULL

#' 

#### functions for subtypePrediction_distributed.R ####

#' function for boxplots of correlation by subtype
#' @param y description
#' @param short description
#' @param pro description
#' @noRd
myplot<-function(y,short,pro){
  par(mfrow=c(3,2),mar=c(5,3,2,2),las=3)
  y$prediction<-factor(y$prediction,levels=c("Basal","Her2","LumA","LumB","Normal"))
  boxplot(y$distances[,1]~y$prediction,border=8,ylab="Basal Correlation",main=paste(short,": Basal",sep=""))
  stripchart(y$distances[,1]~y$prediction,vertical=T,method="jitter",pch=3,add=T)
  boxplot(y$distances[,2]~y$prediction,border=8,ylab="Her2 Correlation",main=paste(short,": Her2",sep=""))
  stripchart(y$distances[,2]~y$prediction,vertical=T,method="jitter",pch=3,add=T)
  boxplot(y$distances[,3]~y$prediction,border=8,ylab="LumA Correlation",main=paste(short,": LumA",sep=""))
  stripchart(y$distances[,3]~y$prediction,vertical=T,method="jitter",pch=3,add=T)
  boxplot(y$distances[,4]~y$prediction,border=8,ylab="LumB Correlation",main=paste(short,": LumB",sep=""))
  stripchart(y$distances[,4]~y$prediction,vertical=T,method="jitter",pch=3,add=T)
  boxplot(y$distances[,5]~y$prediction,border=8,ylab="Normal Correlation",main=paste(short,": Normal",sep=""))
  stripchart(y$distances[,5]~y$prediction,vertical=T,method="jitter",pch=3,add=T)
  boxplot(pro~y$prediction,border=8,ylab="Proliferation Index",main=paste(short,": Proliferation Index",sep=""))
  stripchart(pro~y$prediction,vertical=T,method="jitter",pch=3,add=T)
}

#' function for central median
#' @param x gene expression matrix
#' @noRd 
medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

#' function for quantile central
#' this is adapted from genefu package
#' @param x gene expression matrix or vector
#' @noRd 
rescale <- function(x, na.rm=FALSE, q=0) {
    if(q == 0) {
      ma <- max(x, na.rm=na.rm)
      mi <- min(x, na.rm=na.rm)
    } else {
      ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
      mi <- quantile(x, probs=q/2, na.rm=na.rm)
    }
    xx <- (x - mi) / (ma - mi)
    attributes(xx) <- list("names"=names(x), "q1"=mi,"q2"=ma)
    return(xx)
  }

#' function for calibration methods
#' @noRd 
docalibration = function( y, df.al,calibration = "None", internal=internal, external=external){
  mq = 0.05 ## presetting in genefu robust model
  switch( calibration,
          "None" = {print("No calibration") }, # genefu none
          "Internal" = { ## internal
            if(internal == "medianCtr"){ y = medianCtr(y) } ## parker default method
            else if(internal == "meanCtr") {y = scale(y, center=TRUE, scale=TRUE) } ## "scale" in genefu method
            else if(internal == "qCtr") {y = apply(y, 2, function(x) { return((rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2) }) } ## "robust" in genefu method; mp = 0.05
            else if(internal == internal) { ## which column to use
              medians =  readarray(df.al)
              #print(paste("calibration to:",internal))
              tm<-overlapSets(medians$xd,y)
              y<-(tm$y-tm$x[,internal])
              }
            else { print( "Please choose internal calibration stategy: medianCtr, meanCtr, qCtr or ER+/- relevant ")}
          },
          "External" = { ## external
            medians =  readarray(df.al) 
            #print(paste("calibration to:",external)) ## pre-prepared medians and givenmedians
            tm<-overlapSets(medians$xd,y)
            y<-(tm$y-tm$x[,external]) }
          
  )
  
  return(as.matrix( y) )
  
}



#' function for mapping ID and supplementing missing data if necessary
#' @param x Gene expression matrix
#' @param y Feature data provided by user. The table should contain at least three column, which are probe(probeid or transcriptID), EntrezGene.ID and symbol. 
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select IQR for Affy and select mean for Agilent)
#' @param impute Logic. Please specify if there are NA data adn want keep them
#' @param verbose Logic. 
#' @noRd
domapping = function(x ,y, method = "mean", mapping = TRUE,impute = TRUE, verbose = TRUE ){
       
  # ## test data for microarray
  # library(genefu)
  # # load VDX dataset
  # data(vdxs)
  # x = t(data.vdxs)
  # y = annot.vdxs # "probe" "EntrezGene.ID"
  # method="iqr"
  # mapping=TRUE
  # impute = TRUE
  # verbose = TRUE
  # ## test data if NA cells exist to impute
  # x[c(1,2),c(1,2)] = NA
  # 
  # 
  # ## test data for RNAseq
  # x = data ## from IBC_ensembl TCGA cohort
  # y = anno_feature ## "SYMBOL"   "ENTREZID"
  # method="mean"
  # mapping = FALSE
  # impute = TRUE
  # verbose = TRUE
  # ## test data if NA cells exist to impute
  # #x[c(1,2),c(1,2)] = NA
  
  ## loading genes.signature
  data("genes.signature")
  
  ## first step 
  ## for empty cells. imput or not ?
  if(sum(apply(x,2,is.na))>0 & impute){
    
    library(impute)
    if(verbose){
      probeid_NA = rownames(x)[rowSums(is.na(x))]
      sample_NA = colnames(x)[colSums(is.na(x))]
      print(paste0("The imput objects: ", probeid_NA, " in ", sample_NA))
    }
    
    x = impute.knn(x)
    x = x$data
  }
  
  ## no provided y and !mapping(RNAseq)
  if(length(y) == 0 & !mapping ){
    y =AnnotationDbi::select(org.Hs.eg.db, keys =rownames(x), columns = c( "ENTREZID","SYMBOL"), keytype='SYMBOL' ) 
  } else if(length(y) == 0 & mapping) {
    print("Please provide feature annotation")
  }
  
  ## second step 
  ## if mapping (microarray or transcript )
  
  if( mapping){
    
    probeid = rownames(x)
    entrezid = y$EntrezGene.ID
    names(entrezid) = y$probe
    
    ## process probeid in input data
    entrezid = entrezid[probeid]
    ##remove NA
    entrezid = entrezid[!(is.na(entrezid))]
    x = x[names(entrezid),]
    entrezid = factor( entrezid, levels =  unique(entrezid) ) ## names are unique probeid and content are redundant entrezid 
  
    
    ## This is for probeID or transcriptID 
    ## split expression matrix
    split_mat <- split( as.data.frame(x), entrezid, drop = F)
    
    # function to calculate the desired statistic
    calculate_stat <- function(mat, method) {
      switch(method,
             "mean" = apply(mat, 2, mean, na.rm = TRUE),
             "median" = apply(mat, 2, median, na.rm = TRUE),
             "stdev" = {
               temp = mat
               stdevs = apply(mat, 1, sd, na.rm = TRUE)
               vals =  temp[match(max(stdevs),stdevs),]
               return( vals)
             },
             "iqr" = {
               temp = mat
               iqrs = apply(mat, 1, function(x) { quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE) })
               vals =  temp[match(max(iqrs),iqrs),]
               return( vals)
             }
      )
    }
    
    ## keep processed x
    x = mapply( calculate_stat, split_mat, MoreArgs = list(method = method), SIMPLIFY = T, USE.NAMES = T)
    x = as.data.frame(t(x))
    
    ##print necessary information
    ##Parker
    missing_ID_parker = setdiff( IBC.parker$genes.sig50$EntrezGene.ID, rownames(x) )
    if( length(missing_ID_parker) == 0 & verbose ){ 
      print("PAM50 signatures are covered")
    } else if(verbose) {
      print("These signatures are missing :")
      print(missing_ID_parker)
    }
    
    ##AIMS
    missing_ID_AIMS = setdiff( genes.signature[ genes.signature$AIMS_based == "Yes",]$EntrezGene.ID, rownames(x) )
    if( length(missing_ID_AIMS) == 0 & verbose ){ 
      print("AMIS-based signatures are covered")
    } else if(verbose) {
      print("These signatures are missing for AIMS-based methods :")
      print(missing_ID_AIMS)
    }
    
    } else { ## mapping FALSE 
    
    ## for RNAseq, gene expression at gene level
    ## check if y is empty

    x = merge(y, x, by.y = "row.names", by.x = "SYMBOL", all.x = TRUE)
    

    genes.signature_check = separate_rows(genes.signature, Alias, sep = ", ")
    
    ## find entrezID first
    x_temp_entrezid = x[ x$ENTREZID %in% genes.signature_check$EntrezGene.ID, ]
    
    missing_ID = setdiff( genes.signature_check$EntrezGene.ID, x_temp_entrezid$ENTREZID )
    
    if( length(x$ENTREZID[is.na(x$ENTREZID)])>0 ){
      
      symbol.input = x$SYMBOL[is.na(x$ENTREZID)] 
      
      res_ID_missing = sapply(symbol.input, function(symbol){
        
        ID = genes.signature_check[(symbol == genes.signature_check$Symbol ) | 
                              (symbol == genes.signature_check$Alias ), ]$EntrezGene.ID
        
        if(length(ID) ==0){ID=NULL}
        
        return(ID)
      },simplify = FALSE, USE.NAMES = TRUE)
      
      ## assign entrezID
      x[is.na(x$ENTREZID),]$ENTREZID = res_ID_missing
      
    }
    

    ## find entrez ID again
    x_temp_entrezid = x[ x$ENTREZID %in% genes.signature_check$EntrezGene.ID ,]
    rownames(x_temp_entrezid) = x_temp_entrezid$ENTREZID
    x = x_temp_entrezid[,-c(1,2)]

    ## create two parker matrix and AIMS matrix
    ## print necessary information
    ## process data
    ## print necessary info
    ## Parker
    missing_ID_parker = setdiff( IBC.parker$genes.sig50$EntrezGene.ID, rownames(x) )
    if( length(missing_ID_parker) == 0 & verbose ){ 
      print("PAM50 signatures are covered")
    } else if(verbose) {
      print("These signatures are missing :")
      print(missing_ID_parker)
    }
    
    ##AIMS
    missing_ID_AIMS = setdiff( genes.signature[ genes.signature$AIMS_based == "Yes",]$EntrezGene.ID, rownames(x) )
    if( length(missing_ID_AIMS) == 0 & verbose ){ 
      print("AMIS-based signatures are covered")
    } else if(verbose) {
      print("These signatures are missing for AIMS-based methods :")
      print(missing_ID_AIMS)
    }
    
    }
  
  
    ## get matrix for parker (symbol as colnames)
    x_parker = x[rownames(x) %in% as.character(IBC.parker$genes.sig50$EntrezGene.ID),]
    #x_parker_temp = x_parker
    rownames(x_parker) = IBC.parker$genes.sig50$Symbol[which( rownames(x_parker) %in% IBC.parker$genes.sig50$EntrezGene.ID ) ]

    ## get matrix for AIMS (entrezID as colnames)
    x_AMIS = x[rownames(x) %in% as.character( genes.signature[ genes.signature$AIMS_based == "Yes",]$EntrezGene.ID) ,]
  
    result = list(x_parker = x_parker,x_AMIS = x_AMIS )
  
    return(result)
  }

#' function for PCA plots
#' @noRd
pcaEA<-function(x,classes,size=1,showLegend=T,legendloc="topright",mainStr="",startPC=1,stopPC=2,showNames=T,showClasses=F,axisExpansion=0,groupColors=NA){
  
  features<- dim(x)[1]
  samples<- dim(x)[2]
  sampleNames<- dimnames(x)[[2]]
  featureNames<-dimnames(x)[[1]]
  x<-apply(x,2,as.numeric)
  
  #principal components plots
  data.pca<-prcomp(as.matrix(x))
  
  # Proportion of total variance distributed over 10 first components:
  tmp<-data.pca$sdev[1:10]^2/sum(data.pca$sdev^2)
  
  gr.labels<-as.vector(t(classes))
  gr.labels.fac<-factor(as.vector(t(classes)),exclude="")
  nlabels<-nlevels(gr.labels.fac)
  legendLabels<-vector()
  legendColors<-vector()
  for(k in 1:nlabels){
    group<-levels(gr.labels.fac)[k]
    legendLabels[k]<-group
    if(length(groupColors)>1){
      gr.labels[gr.labels.fac==group]<-groupColors[k]
      legendColors[k]<-groupColors[k]
    }else{
      gr.labels[gr.labels.fac==group]<-k
      legendColors[k]<-k
    }
  }
  
  if(length(groupColors)==1){
    gr.labels<-as.numeric(gr.labels)
  }
  
  #plot 2pcs by each other
  i<-startPC
  j<-stopPC
  
  #graphing parameters
  par(lab=c(3,4,3))
  par(mgp=c(.3,.5,.0))
  par(mai=c(.5,.5,.5,.5))
  par(xaxt="n",yaxt="n")
  
  strM<-mainStr
  strX<-paste("PC",i,paste("(",round(tmp[i],4)*100,"%)",sep=""),sep=" ")
  strY<-paste("PC",j,paste("(",round(tmp[j],4)*100,"%)",sep=""),sep=" ")
  xmin<-min(data.pca$rotation[,i])-abs(axisExpansion*min(data.pca$rotation[,i]))
  xmax<-max(data.pca$rotation[,i])+abs(axisExpansion*max(data.pca$rotation[,i]))
  ymin<-min(data.pca$rotation[,j])-abs(axisExpansion*min(data.pca$rotation[,j]))
  ymax<-max(data.pca$rotation[,j])+abs(axisExpansion*max(data.pca$rotation[,j]))
  plot(data.pca$rotation[,i],data.pca$rotation[,j], xlab=strX, ylab=strY, 
       main=strM, col=gr.labels,cex=size,pch="",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
  if(showNames){
    text(data.pca$rotation[,i],data.pca$rotation[,j],labels=names(data.pca$rotation[,i]),cex=size*.6)
  }else{
    if(showClasses){
      text(data.pca$rotation[,i],data.pca$rotation[,j],labels=gr.labels.fac,cex=size*.6)
    }else{
      points(data.pca$rotation[,i],data.pca$rotation[,j],col=gr.labels,cex=size*1.5,pch=19)
    }
    if(showLegend){
      legend(legendloc,legend=legendLabels,col=legendColors,pch=19,x.intersp=.3,yjust=.5,bty="n",cex=size)
    }
  }
}

#' assign color for heatmap
#' @noRd
cols <- function(lowi = "green", highi = "red", ncolors = 20) {
  low <- col2rgb(lowi)/255
  high <- col2rgb("black")/255
  col1 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                       high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  low <- col2rgb("black")/255
  high <- col2rgb(highi)/255
  col2 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                       high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  col<-c(col1[1:(ncolors-1)],col2)
  return(col)
}

#' function for heatmap plot
#' @noRd
myHeatmap<-function(x,t.colors=NA,fileName="cluster.cdt",linkage="average",distance="spearman",contrast=2,returnSampleClust=F,rowNames=NA,rightMar=7,bottomMar=1,colNames=NA){
  
  temp<-hclust2treeview(x,method=distance,file=fileName,link=linkage,keep.hclust=T)
  gTree<-temp[[1]]
  sTree<-temp[[2]]
  
  imageVals<-x
  imageVals[x > contrast] <- contrast
  imageVals[x < -1 * contrast] <- -1 * contrast
  
  if(sum(is.na(t.colors))>0){
    heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
            col=cols(),labCol=colNames, scale="none",
            margins=c(bottomMar,rightMar),labRow=rowNames)
  }else{
    if(length(t.colors)>dim(imageVals)[2]){
      heatmap.plus(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
                   col=cols(),labCol=colNames,labRow=rowNames,scale="none",
                   ColSideColors=t.colors, margins=c(bottomMar,rightMar))
    }else{
      heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
              col=cols(),labCol=colNames,labRow=rowNames,scale="none",
              ColSideColors=as.vector(t(t.colors)), margins=c(bottomMar,rightMar))
    }
  }
  if(returnSampleClust){
    return(sTree)
  }
}


#' function for standardize
#' @param x description
#' @noRd
standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}


#' Function for ordering gene in expression matrix as PAM50 gene signature 
#' @param x PAM50 centroid matrix
#' @param y enpression matrix
#' @noRd
overlapSets<-function(x,y){
  
  # subset the two lists to have a commonly ordered gene list
  x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
  y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]
  
  #and sort such that thing are in the correct order
  x<-x[sort.list(row.names(x)),]
  y<-y[sort.list(row.names(y)),]
  
  return(list(x=x,y=y))
}

#'
#' Function for de-duplicated genes in gene expression matrix
#' method : mean, median, stdev, iqr
#' @noRd 
collapseIDs<-function(x,method="mean"){
  
  allids<-as.vector(row.names(x))
  ids<- levels(as.factor(allids))
  x.col<- NULL
  
  if(length(ids)==dim(x)[1]){ 
    dimnames(x)[[1]]<-allids
    return(x) 
  }
  
  for(i in 1:length(ids)){
    if(sum(allids==ids[i])>1){
      indices <- allids==ids[i] 
      if(method=="mean"){
        vals<-apply(x[indices,],2,mean,na.rm=T)
      }
      if(method=="median"){
        vals<-apply(x[indices,],2,median,na.rm=T)
      }
      if(method=="stdev"){   
        temp<- x[indices,]
        stdevs<- apply(temp,1,sd,na.rm=T)
        vals<- temp[match(max(stdevs),stdevs),]
      }
      if(method=="iqr"){   
        temp<- x[indices,]
        iqrs<- apply(temp,1,function(x){quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T)})
        vals<- temp[match(max(iqrs),iqrs),]
      }
      x.col <- rbind(x.col,vals)
    }else{
      x.col <- rbind(x.col,x[allids==ids[i],])
    }
    
    # can be mean or iqr (probe with max iqr is selected)
    # typically, mean is preferred for long oligo (microarray) and
    # iqr is preferred for short oligo platforms (affy)
    # ID should be unique by RNAseq
    
  }
  
  dimnames(x.col)<- list(ids,dimnames(x)[[2]])
  return(x.col)
  
}

#' Function for new data structure
#' @param data gene expression matrix or median of train dataset
#' @noRd
readarray<-function(data,designFile=NA,impute=T,method="mean"){
  
  # do de-duplicated for gene
  ## But not working for RNAseq, only for microarray sequencing data
  
  xd = collapseIDs(data,method)
  
  features<- dim(xd)[1]
  samples<- dim(xd)[2]
  geneNames<-rownames(xd)
  sampleNames = colnames(xd)
  xd<-apply(xd,2,as.numeric)
  rownames(xd)<-geneNames
  colnames(xd)<-sampleNames
  classes = NULL
  
  if(!is.na(designFile)){
    x<-read.table(designFile,sep="\t",header=T,row.names=1,fill=T,stringsAsFactors=FALSE)
    xd<-xd[,sort.list(colnames(xd))]
    xd<-xd[,colnames(xd) %in% rownames(x)]
    x<-x[rownames(x) %in% colnames(xd),]
    x<-x[sort.list(rownames(x)),]
    classes<-as.data.frame(x)
  }

  if(sum(apply(xd,2,is.na))>0 & impute){
    library(impute)
    allAnn<-dimnames(xd)
    data.imputed<-impute.knn(as.matrix(xd))$data
    xd<-data.imputed[1:features,]
    dimnames(xd)<-allAnn
  }
  
  return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}

#' Function for suffix of medians for gene centering
#' @noRd 
getsurffix = function( calibration,internal=internal, external=external){
  

  if(calibration == "None" ){surffix =calibration } else{
    switch( calibration,
            "None" = {surffix = calibration},
            "Internal" = {surffix = internal },
            "External" = { surffix = external }
    )
    }
  return(surffix)
}

#' Function for consensus subtype
#' @noRd 
get_consensus_subtype <- function(patient_row) {
  patient_row = unlist(patient_row, use.names = FALSE)
  counts <- table(patient_row)
  max_subtype <- names(counts)[which.max(counts)]
  return(max_subtype)
}


#' Function to obtain entropy of each subtype per patient
#' @noRd 
get_entropy_subtype <- function(patient_row) {
  patient_row = unlist(patient_row, use.names = FALSE)
  probs = table(patient_row) / length(patient_row)
  res = -sum(probs * log2(probs), na.rm = TRUE)
  return(res)
  
  }

#' Function to get the average correlation and ROR
get_average_subtype <- function(res_ihc_iterative, consensus_subtypes, hasClinical =F) {

  ## correlation and ROR to be averaged
  ## if hasclini ?? need to be added later
  sum_colnames = c("Basal","Her2","LumA", "LumB", "Normal")
  
  all_patients = names(res_ihc_iterative[[1]]$predictions )
  
  sum_cols_list = mapply(function(res_ihc){
    
    #res_ihc = as.data.frame( res_ihc_iterative[[1]])
    
    res_ihc$distances = as.data.frame(res_ihc$distances )
    
    ## if FALSE, make the cell as NULL
    keep = res_ihc$predictions == consensus_subtypes
    res_ihc$distances[ !keep, ] = as.list(rep(NA, 5 ))
    
    res = res_ihc$distances %>%
      mutate_at(vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
    
    return(res )
    
  }, res_ihc_iterative, SIMPLIFY = FALSE, USE.NAMES = FALSE )
  
  
  ## count_na
  count_weight_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
    x[!is.na(x)] <- 1
    x[is.na(x)] <- 0
    return(x)
  }))
  
  ## sum all for each cell
  sum_cols_save = Reduce(`+`, lapply(sum_cols_list, function(x) {
    
    ## change NA cell to 0 cell
    x[is.na(x)] = 0
    
    return(x)
  }))
  
  ## get the mean for each cell
  ## only when subtype is supported by consensus_subtypes for each iteration and each patient
  mean_cols_save = sum_cols_save / count_weight_save
  
  ## get the mean testdata
  sum_cols_list.testdata <- mapply(function(res_ihc){
    res_ihc$testData
  }, res_ihc_iterative, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  ## sum all for each cell
  mean_cols_save.testdata = Reduce(`+`, sum_cols_list.testdata) / length(sum_cols_list.testdata)

 
  ## when hasClinical
  if(hasClinical) { 
  
    sum_cols_list.prosigna = mapply(function(res_ihc){
      
      #res_ihc = res_ihc_iterative[[1]]
      
      res_ihc$distances.prosigna = as.data.frame(res_ihc$distances.prosigna )
      
      ## if FALSE, make the cell as NULL
      keep = res_ihc$predictions == consensus_subtypes
      res_ihc$distances.prosigna[ !keep, ] = as.list(rep(NA, 5 ))
      
      res = res_ihc$distances.prosigna %>%
        mutate_at(vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
      
      return(res )
      
    }, res_ihc_iterative , SIMPLIFY = FALSE, USE.NAMES = FALSE )
    
    
    ## sum all for each cell
    sum_cols_save.prosigna = Reduce(`+`, lapply(sum_cols_list.prosigna, function(x) {
      
      ## change NA cell to 0 cell
      x[is.na(x)] = 0
      
      return(x)
    }))
    
    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes for each iteration and each patient
    sum_cols_save.prosigna = sum_cols_save.prosigna / count_weight_save
    
    
    res = list( mean_distance = mean_cols_save, mean_distance.prosigna = sum_cols_save.prosigna, testdata =  mean_cols_save.testdata) 
    
    
  } else{
    
    res = list(mean_distance = mean_cols_save,testdata =  mean_cols_save.testdata)
    
  }
  
  return(res)
}

#' Function for predict PAM50 subtyping
#' @param x median train file
#' @param y gene expression matrix
#' @param classes description
#' @param nGenes None
#' @param priors description
#' @param distm "euclidean" or "spearman"
#' @param std logical value
#' @param centrids logical value
#' @noRd
sspPredict<-function(x, classes="", y, nGenes="", priors="equal",std=F, distm="euclidean",centroids=F, hasClinical = F){
  
  # ## test data
  # x = IBC.parker$centroid
  # y = mat
  # classes=""
  # nGenes=""
  # priors="equal"
  # std=F
  # distm="spearman"
  # centroids=T

  dataMatrix<-x
  features<- dim(x)[1]
  samples<- dim(x)[2]
  sampleNames<- dimnames(x)[[2]]
  featureNames<- dimnames(x)[[1]]
  
  #parse the test file - same as train file but no rows of classes
  tdataMatrix<-y
  tfeatures<- dim(y)[1]
  tsamples<- dim(y)[2]
  tsampleNames<- dimnames(y)[[2]]
  tfeatureNames<- dimnames(y)[[1]]
  
  #dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
  temp <- overlapSets(dataMatrix,tdataMatrix)
  dataMatrix <- temp$x
  tdataMatrix <- temp$y
  sfeatureNames<-row.names(dataMatrix)
  
  # standardize both sets
  if(std){
    dataMatrix<-standardize(dataMatrix)
    tdataMatrix<-standardize(tdataMatrix)
  }
  
  if(!centroids){
    thisClass <- as.vector(classes[,1])
    nClasses<-nlevels(as.factor(thisClass))
    classLevels<-levels(as.factor(thisClass))
    for(j in 1:nClasses){
      thisClass[thisClass==classLevels[j]] <- j
    }
    thisClass<-as.numeric(thisClass)
    dataMatrix <- dataMatrix[,!(is.na(thisClass))]
    thisClass <- thisClass[!(is.na(thisClass))]
    
    scores<-apply(dataMatrix,1,bwss,thisClass)
    trainscores<-vector()	
    for(j in 1:dim(dataMatrix)[1]){			
      trainscores[j]<-scores[[row.names(dataMatrix)[j]]]$bss / scores[[row.names(dataMatrix)[j]]]$wss
    }
    
    dataMatrix<-dataMatrix[sort.list(trainscores,decreasing=T),]
    tdataMatrix<-tdataMatrix[sort.list(trainscores,decreasing=T),]	
    
    if(nGenes==""){
      nGenes<-dim(dataMatrix)[1]
    }
    #print(paste("Number of genes used:",nGenes))
    
    dataMatrix<-dataMatrix[1:nGenes,]
    tdataMatrix<-tdataMatrix[1:nGenes,]
    
    centroids<-matrix(nrow=nGenes,ncol=nClasses)
    for(j in 1:nClasses){
      centroids[,j]<-apply(dataMatrix[,thisClass==j],1,mean)
    }
    dimnames(centroids)<-list(row.names(dataMatrix),NULL)
    
  }else{
    nGenes<-dim(dataMatrix)[1]
    #print(paste("Number of genes used:",nGenes))
    centroids<-dataMatrix
    nClasses<-dim(centroids)[2] ## five subtypes; for prosigna, keep normal when calculating
    classLevels<-dimnames(centroids)[[2]]
  }
  
  distances<-matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
  for(j in 1:nClasses){
    if(distm=="euclidean"){
      distances[,j]<- dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
    }
    if(distm=="correlation" | distm=="pearson"){
      distances[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids[,j], x, method = "pearson", use = "pairwise.complete.obs"))
      }
    if(distm=="spearman"){
      distances[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids[,j], x, method = "spearman", use = "pairwise.complete.obs"))
      }
  }
  
  
  prediction = classLevels[apply(distances, 1, which.min,simplify = TRUE)]
  names(prediction)<-tsampleNames
  
  
  if(hasClinical){
    ## genes to be excluded
    genes.ex = c("BIRC5", "CCNB1", "GRB7","MYBL2")
    
    #parse the test file - same as train file but no rows of classes
    tdataMatrix<-y[ which( !(rownames(y) %in%  genes.ex ) ),]

    #dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
    temp <- overlapSets(dataMatrix,tdataMatrix)
    dataMatrix <- temp$x
    tdataMatrix <- temp$y

    nGenes<-dim(dataMatrix)[1]
    #print(paste("Number of genes used for prosigna-like:",nGenes))
    centroids_prosigna<-dataMatrix

    distances.prosigna <-matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
    for(j in 1:nClasses){
      if(distm=="euclidean"){
        distances.prosigna[,j]<-dist(t(cbind(dataMatrix[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
      }
      if(distm=="correlation" | distm=="pearson"){
        distances.prosigna[,j] = apply(tdataMatrix, 2, function(x) -cor(dataMatrix[,j], x, method = "pearson", use = "pairwise.complete.obs"))
        }
      if(distm=="spearman"){
        distances.prosigna[,j] = apply(tdataMatrix, 2, function(x) -cor(dataMatrix[,j], x, method = "spearman", use = "pairwise.complete.obs"))
        }
    }
    
  }
  

  ## return 
  if ( !hasClinical) {
  res = list(predictions=prediction,testData= as.matrix( y),distances=distances,centroids=centroids)
  }else {
    res = list(predictions=prediction,testData= as.matrix(y),distances=distances,distances.prosigna = distances.prosigna,centroids=centroids)
  }
  
  return(res)
}


#' Function for risk 
#' 
#' @param hasClinical provide clinical information for grouping 
#' @param out it is the result of sspPredict() function. 
#' @return ROR, ROR risk group and other indications
#' @noRd
RORgroup = function(out, df.cln , hasClinical = FALSE ){
  
  # ## test data
  # out$distances = -1 * out$distances
  # out$distances.prosigna = -1 * out$distances.prosigna
  # 
  # df.cln$T = rep(c(1,2,3), length.out = 141 )
  # df.cln$NODE = rep(c( 1,3,4,NA,5), length.out = 141 )
  # Clinical = df.cln
  # rownames(Clinical) = Clinical$PatientID
  # hasClinical = TRUE

  
  # ## test data
  # out
  # df.cln = df.pam
  # hasClinical = hasClinical
  # 
  
  Clinical = df.cln[which( df.cln$PatientID %in% names(out$predictions) ),]
  rownames(Clinical) = Clinical$PatientID
  ## prepare outtable
  # 
  # out$distances = -1 * out$distances
  # out$distances.prosigna = -1 * out$distances.prosigna
  # 
  
  distance = data.frame(out$distances, row.names = names(out$predictions) )
  names(distance) = c("Basal","Her2","LumA","LumB","Normal")
  
  Call = data.frame( "Call" = out$predictions, row.names = names(out$predictions))
  
  #providing proliferation signatures
  proliferationGenes<-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
  proliferationGenes.prosigna = c("ANLN", "CEP55", "ORC6L", "CCNE1", "EXO1", "PTTG1", "CDC20", "KIF2C", "RRM2", "CDC6", "KNTC2", "TYMS", "CDCA1", 
                                  "MELK", "UBE2C", "CENPF", "MKI67", "UBE2T")

  ###
  # some constants for ROR groups
  ###
  
  # for subtype only model
  glthreshold<- -0.15
  ghthreshold<-  0.1
  
  # for subtype + proliferation model
  gplthreshold<- -0.25
  gphthreshold<-  0.1
  
  # for combined model
  clthreshold<- -0.1
  chthreshold<-  0.2
  
  # for combined + proliferation model
  cplthreshold<- -0.2
  cphthreshold<-  0.2
  
  ## ??? 
  # for combined + proliferation model + prosigna with negative node status 
  cplthreshold.prosigna.NODE.0 = 40
  cphthreshold.prosigna.NODE.0 = 60
  
  # for combined + proliferation model + prosigna with node status
  ## first 4 or more nodes would be assigned as high-risk
  cplthreshold.prosigna.NODE = 15
  cphthreshold.prosigna.NODE = 40
  
  
  ## ER and HER2 score
  #out$testData = mat
  erScore = out$testData["ESR1",]
  her2Score = out$testData["ERBB2",]
  
  er_her2 = data.frame( "ER" = erScore, "HER2" = her2Score, row.names = colnames(out$testData ))

  # calculate the proliferation score
  prolifScore = apply(out$testData[ which( rownames(out$testData) %in% proliferationGenes), ], 2, mean, na.rm =T )
  
  prolifScore.prosigna = apply(out$testData[ which( rownames(out$testData) %in% proliferationGenes.prosigna) , ], 2, mean, na.rm =T )
  
  ## confidence
  call.conf<-c()
  for(j in 1:length(out$predictions)){
    call.conf[j]= 1-cor.test(out$testData[,j],out$centroids[,which(colnames(out$centroids)==out$predictions[j])],method="spearman", exact=FALSE)$p.value
    
    ## debug
    # tryCatch(
    #   {
    #     call.conf[j]= 1-cor.test(out$testData[,j],out$centroids[,which(colnames(out$centroids)==out$predictions[j])],method="spearman",exact=FALSE)$p.value
    #     
    #   },
    #   warning = function(w) {
    #     message("Warning: Tied data encountered.")
    #     print( colnames( out$testData)[j] )
    #   }
    # )
  }
  
  call.conf<-round(call.conf,2)
  call.conf = data.frame( "Confidence" = call.conf, row.names = names(out$predictions))
  

  ## calculate the risk scores
  ## genomic
  ## basal, her2, lumA, lumB in order
  genomic <- 0.04210193*out$distances[,1] + 0.12466938*out$distances[,2] + -0.35235561*out$distances[,3] + 0.14213283*out$distances[,4]
  
  ## ?? where weights ?
  genomicWprolif <- -0.0009299747*out$distances[,1] + 0.0692289192*out$distances[,2] + -0.0951505484*out$distances[,3] +  0.0493487685*out$distances[,4] + 0.3385116381*prolifScore
  
  
  # threshold the risk score
  griskgroups<-genomic
  griskgroups[genomic>ghthreshold]<-"high"
  griskgroups[genomic>glthreshold & genomic<ghthreshold]<-"med"
  griskgroups[genomic<glthreshold]<-"low"
  gpriskgroups<-genomicWprolif
  gpriskgroups[genomicWprolif>gphthreshold]<-"high"
  gpriskgroups[genomicWprolif>gplthreshold & genomicWprolif<gphthreshold]<-"med"
  gpriskgroups[genomicWprolif<gplthreshold]<-"low"
  
  genomic<- 100* (genomic + 0.35 ) / 0.85
  genomicWprolif<- 100* (genomicWprolif + 0.35 ) / 0.85
  

  ROR.genomic = data.frame("ROR-S (Subtype Only)" = genomic, "ROR-S Group (Subtype Only)" = griskgroups,
                           "ROR-P (Subtype + Proliferation)" = genomicWprolif, "ROR-P Group (Subtype + Proliferation)" = gpriskgroups,
                           check.names = FALSE)
  
  
  ## how to transfer this variable ?
  
  if(hasClinical){
    
    # ROR
    ## combined with T
    ## check order by patinets/rownames
    if ( "T" %in% colnames(Clinical) ) {
      
      xT= as.numeric(as.vector(Clinical$T))
      combined = 0.0442770*out$distances[,1] + 0.1170297*out$distances[,2] + -0.2608388*out$distances[,3] + 0.1055908*out$distances[,4] + 
        0.1813751*xT
      combinedWprolif = -0.009383416*out$distances[,1] +  0.073725503*out$distances[,2] + -0.090436516*out$distances[,3] + 0.053013865*out$distances[,4] + 
        0.131605960*xT + 0.327259375*prolifScore
      
      criskgroups<-combined
      criskgroups[combined>chthreshold]<-"high"
      criskgroups[combined>clthreshold & combined<chthreshold]<-"med"
      criskgroups[combined<clthreshold]<-"low"
      cpriskgroups<-combinedWprolif
      cpriskgroups[combinedWprolif>cphthreshold]<-"high"
      cpriskgroups[combinedWprolif>cplthreshold & combinedWprolif<cphthreshold]<-"med"
      cpriskgroups[combinedWprolif<cplthreshold]<-"low"
      
      combined<- 100* (combined + 0.35 ) / 0.85
      combinedWprolif<- 100* (combinedWprolif + 0.35 ) / 0.85
 
      ROR.combined = data.frame("ROR-C (Subtype + Clinical)" = combined, "ROR-C Group (Subtype + Clinical)" = cpriskgroups,
                                "ROR-PC (Subtype + Clinical + Proliferation)" = combinedWprolif, "ROR-PC Group (Subtype + Clinical + Proliferation)" = cpriskgroups,
                                check.names = FALSE)
      
      
      # ROR prosigna
      ## log2 of nCounter expression, FFPE
      ## log2 of FPKM, FF , illumina
      ## 46 genes, 
      ## 46 genes to calculate correlation or no
      combinedWprolif.prosigna = 54.7690 * (-0.0067 * out$distances.prosigna[,1] + 0.4317*out$distances.prosigna[,2] - 0.3172*out$distances.prosigna[,3] + 
                                              0.4894*out$distances.prosigna[,4] + 0.1981*prolifScore.prosigna + 0.1133*xT + 0.8826)
      
     
      ##grouping for ROR prosigna 
      if( "NODE" %in% colnames(Clinical) ) {
        
        
        cpriskgroups.prosigna = combinedWprolif.prosigna
        
        ## check if NODE 
        ## grouping by NODE
        patients.NODE = rownames(Clinical)[ !is.na(Clinical$NODE )]

        if( length( patients.NODE ) > 0  )  {
        
          ## NODE >= 4
          patients.NODE.4 = rownames(Clinical)[ Clinical$NODE >= 4 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.4 ) > 0  ) { 
            
            cpriskgroups.prosigna[patients.NODE.4] = "high"
          }
          
          ## NODE < 4 & >0
          patients.NODE.3 = rownames(Clinical)[ Clinical$NODE < 4 & Clinical$NODE > 0 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.3 ) > 0  ) { 
            
            temp = combinedWprolif.prosigna[patients.NODE.3]
            cpriskgroups.prosigna[patients.NODE.3][ temp> cphthreshold.prosigna.NODE ] <-"high"
            cpriskgroups.prosigna[patients.NODE.3][temp>cplthreshold.prosigna.NODE & temp<cphthreshold.prosigna.NODE ]<-"med"
            cpriskgroups.prosigna[patients.NODE.3][temp<cplthreshold.prosigna.NODE]<-"low"
            
          }
          
          ## negative node
          patients.NODE.0 = rownames(Clinical)[ Clinical$NODE == 0 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.0 ) > 0  ) { 
            
            temp = combinedWprolif.prosigna[patients.NODE.0]
            cpriskgroups.prosigna[patients.NODE.0][ temp> cphthreshold.prosigna.NODE.0 ] <-"high"
            cpriskgroups.prosigna[patients.NODE.0][temp>cplthreshold.prosigna.NODE.0 & temp<cphthreshold.prosigna.NODE.0 ]<-"med"
            cpriskgroups.prosigna[patients.NODE.0][temp<cplthreshold.prosigna.NODE.0]<-"low"
            
          }
          
          
          ## check if NA NODE 
          ## grouping by NA NODE
          patients.NA = rownames(Clinical)[ is.na(Clinical$NODE ) | is.na(Clinical$T) ]
          
          if ( length( patients.NA ) > 0  )
            
          temp = combinedWprolif.prosigna[patients.NA]
          cpriskgroups.prosigna[patients.NA] = NA

        }

        ROR.combined.prosigna = data.frame("ROR-PC (Subtype + Clinical + Proliferation.prosigna)" = combinedWprolif.prosigna, "ROR-PC Group (Subtype + Clinical + Proliferation.prosigna)" = cpriskgroups.prosigna,
                                           check.names = FALSE)
        
      } else {
        
        cpriskgroups.prosigna[combinedWprolif.prosigna> cphthreshold.prosigna]<-"high"
        cpriskgroups.prosigna[combinedWprolif.prosigna>cplthreshold.prosigna & combinedWprolif.prosigna<cphthreshold.prosigna]<-"med"
        cpriskgroups.prosigna[combinedWprolif.prosigna<cplthreshold.prosigna]<-"low"
        
        
        ROR.combined.prosigna = data.frame("ROR-PC (Subtype + Clinical + Proliferation.prosigna)" = combinedWprolif.prosigna, "ROR-PC Group (Subtype + Clinical + Proliferation.prosigna)" = cpriskgroups.prosigna,
                                  check.names = FALSE)
        
      }
      
      
      outtable = cbind( distance, Call, call.conf, ROR.genomic,ROR.combined, ROR.combined.prosigna,  er_her2)
      
      
    }
    
    
  } else {
    
    outtable = cbind( distance, Call, call.conf, ROR.genomic, er_her2)
    
  }
  
  return(outtable)
}


#' Function for visualization

visualition = function(y,out,short ){
  
  ## for visualization
  x = PCA_PAM50$x #mediansPerDataset_v2.txt
  
    
  inputDir = paste( getwd(), "Call_wth_IHCMd", sep = "/") 
  
  pdfname1<-paste(inputDir,paste("predictionScores_pam50RankCorrelation_1_",short,".pdf",sep=""),sep="/")
  pdfname2<-paste(inputDir,paste("predictionScores_pam50RankCorrelation_2_",short,".pdf",sep=""),sep="/")
  clustername<-paste(inputDir,paste(short,"_PAM50_normalized_heatmap",sep=""),sep="/")
  outFile<- paste(inputDir,paste(short,"_pam50scores.txt",sep=""),sep="/")
  
  # make some plots for evaluation
  print(paste("ER range:",quantile(erScore,.9,na.rm=T)-quantile(erScore,.1,na.rm=T)))
  
  subtypeColors<-out$predictions
  subtypeColors[subtypeColors=="Basal"]<-"red"
  subtypeColors[subtypeColors=="Her2"]<-"hotpink"
  subtypeColors[subtypeColors=="LumA"]<-"darkblue"
  subtypeColors[subtypeColors=="LumB"]<-"skyblue"
  subtypeColors[subtypeColors=="Normal"]<-"green"
  conf.colors<-call.conf
  conf.colors[call.conf>=0.95]<-"black"
  conf.colors[call.conf<0.95]<-"red"
  
  ## if make plots
  pdf(paste(clustername,".pdf",sep=""))
  myHeatmap(out$testData,cbind(subtypeColors,conf.colors),file=paste(clustername,".cdt",sep=""),rowNames=rownames(out$testData))
  dev.off()
  
  pdf(pdfname1,height=10,width=12)
  pars<-par(no.readonly=T)
  myplot(out,short,prolifScore)
  dev.off()
  
  tm<-overlapSets(x$xd,y$xd)
  tm$x<-tm$x[,!is.na(x$classes$subtype)]
  tm<-cbind(tm$x,	impute.knn(as.matrix(tm$y))$data)
  classes<-matrix(nrow=4,ncol=dim(tm)[2])
  nTrainSamples<-length(x$classes$subtype[!is.na(x$classes$subtype)])
  
  classes[1,]<-c(rep("train",nTrainSamples),rep(short,dim(y$xd)[2]))
  classes[2,]<-c(x$classes$subtype[!is.na(x$classes$subtype)],rep(NA,dim(tm)[2]-dim(x$xd[,!is.na(x$classes$subtype)])[2]))
  classes[3,]<-c(rep(NA,dim(tm)[2]-length(out$predictions)),out$predictions)
  
  pdf(pdfname2,height=6,width=12)
  tm<-scale(tm,center=F)
  par(mfrow=c(1,3))
  pcaEA(tm,classes[1,],mainStr="Traing and Test sets",showNames=F,showClasses=F)
  pcaEA(tm[,!is.na(as.vector(t(classes[2,])))],classes[2,!is.na(as.vector(t(classes[2,])))],mainStr="Training cases",showNames=F,showClasses=F,groupColors=c("red","hotpink","darkblue","skyblue","green"))
  pcaEA(tm[,!is.na(as.vector(t(classes[3,])))],classes[3,!is.na(as.vector(t(classes[3,])))],mainStr="Test cases",showNames=F,showClasses=F,groupColors=c("red","hotpink","darkblue","skyblue","green"))
  par(pars)
  dev.off()
    
  
}


#' 
#' Function for predicting PAM50 intrinsic subtypes and calculation of proliferation 
#' 


####function to make calls using near-centroid strategy #### 
#' Function for calling PAM50 subtypes by parker based methods
#' Here, we integrated parker-based methods and genefu PAM50 model
#' @param mat gene expression matrix, log of normalized
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd
#' 

makeCalls.parker = function(mat,df.cln, calibration = "None", internal = NA,external=NA, medians = NA, hasClinical =FALSE  ){
  
  ####
  # run the assignment algorithm
  # be careful to choose calibration strategy
  ####
  
  # mat = log2(Test.matrix + 0.001 )
  # df.cln = df.cln
  # calibration = "None"
  # internal = NA
  # external= NA
  # medians = NA
  # hasClinical = FALSE
  # 
  # mat = data_input$x_parker
  # df.cln = clinic
  # calibration = "Internal"
  # internal = "meanCtr"
  # external= NA
  # medians = NA
  # hasClinical = T
  
  
  fl.mdn = IBC.parker$medians
  ## prepare df.al
  
  if(calibration == "External" & external == "Given.mdns" ) {
    
    if( length(medians) == 1 || is.na( medians) ){
      stop("Please input prepared medians, genes in the first column and median values in the second column")
    } else {
      
      colnames(medians) = c("X","Given.mdns")
      
      df.al = merge(fl.mdn, medians, by = "X")
      rownames(df.al) = df.al$X
      df.al = df.al[,-1]
      
      surffix = external 
    }
    
  } else {
    
    df.al = fl.mdn
    
    rownames(df.al) = df.al$X
    df.al = df.al[,-1]
    
  }
  
  
  ## run 
  
  # normalization
  mat = docalibration( mat, df.al, calibration, internal=internal, external=external)
  
  out = sspPredict(IBC.parker$centroid, classes="", mat, std=F, distm="spearman", centroids=T, hasClinical = hasClinical)
  
  if(hasClinical) {
    out$distances.prosigna =  -1 * out$distances.prosigna
  }
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical )
  
  Int.sbs = data.frame(PatientID = names(out$predictions), Int.SBS.ihc = out$predictions, IHC = df.cln$IHC , row.names = NULL )
  
  return(list(Int.sbs=Int.sbs, score.fl=out.ROR, mdns.fl= df.al, outList=out))
  
  
}



#### function to form a ER-balance subet and derive its median---write it to PAM50 dir

#' Function to form a ER-balance subet and derive its median
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd
makeCalls.ihc = function(mat, df.cln, seed=118,calibration = "Internal", internal = "IHC.mdns", external=NA, medians = NA , hasClinical = FALSE){
  # message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  # message("##IHC subtype column should be named ---IHC---")
  # 
  # ##test data for TCGA Cell2015
  # mat = data_input$x_parker
  # df.cln = clinic
  # seed=118
  # calibration = "Internal"
  # internal = "IHC.mdns"
  # hasClinical = TRUE
  # 
  ERN.ihc = df.cln[which(df.cln$ER == "ER-"),] ### get ER- samples data.frame
  dim(ERN.ihc)	#[1] 153   9
  
  ERP.ihc = df.cln[which(df.cln$ER == "ER+"),]
  dim(ERP.ihc) #[1] 559   9
  
  #seed = 118
  set.seed(seed);i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1]) # take equal number of ER+ and ER- samples
  length(ERP.ihc$PatientID[i]) # ER positive samples
  length(ERN.ihc$PatientID)    # ER negative samples
  
  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  #mat = Test.matrix
  mbal.ihc = mat[,c(ERP.ihc$PatientID[i],ERN.ihc$PatientID)]
  
  dim(mbal.ihc) #[1]  50 306
  
  # Find median
  surffix = getsurffix(calibration = calibration, internal,external)
  
  mdns      = apply(mbal.ihc,1,median,na.rm=T) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
  
  colnames(df.mdns) = c("X",surffix)
  
  ## merge mdns
  fl.mdn = IBC.parker$medians
  
  df.al = merge(fl.mdn, df.mdns , by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  ## centroids
  centroids = IBC.parker$centroid #pam50_centroids.txt
  

  ## normalization
  mat = docalibration( mat, df.al, calibration,internal)
  
  out = sspPredict( centroids, classes="", mat, std=F, distm="spearman", centroids=T, hasClinical = hasClinical)
  
  if(hasClinical) {
    out$distances.prosigna =  -1 * out$distances.prosigna
  }
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical )
  
  Int.sbs = data.frame(PatientID = names(out$predictions), Int.SBS.ihc = out$predictions, IHC = df.cln$IHC , row.names = NULL )
  
  return(list(Int.sbs=Int.sbs, score.fl=out.ROR, mdns.fl= df.al,outList=out))
  
}



#' Function for iterative ER subset gene centering 
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param iterative times to predict subtypes
#' @param ratio The options are either 1:1 or 54(ER+):64(ER-). The latter is ER ratio used for UNC230 train cohort
#' @param calibration surffix name for output
#' @param hasClinical provide clinical information, default is NA. 
#' @noRd

makeCalls.ihc.iterative = function( mat, df.cln, iterative = 100, ratio = 54/64, calibration = "Internal", internal = "ER.mdns", external=NA, medians = NA , hasClinical = FALSE){

  # ## test data
  # mat = Test.matrix
  # ratio=1
  # ratio= 54/64
  # ##  46% ER-positive (54/118) and 54% ER-negative (64/118) from Zhao paper ,
  # dependent = "ER.mdns"
  # calibration = "Internal"
  # iterative = 5
  # hasClinical = T
  
  seed = 118
  set.seed(seed)
  
  
  # load the published centroids for classifcation
  centroids = IBC.parker$centroid #pam50_centroids.txt
  
  ## preprocess the input matrix
  ### get ER- samples
  ERN.ihc = df.cln[which(df.cln$ER =="ER-"),] ### get ER- samples data.frame
  dim(ERN.ihc)	#[1] 30 2
  
  ### get ER+ samples
  ERP.ihc = df.cln[which(df.cln$ER =="ER+"),]
  dim(ERP.ihc) #[1] 111 2
  
  ## check the ER composition
  if(dim(ERP.ihc)[1] < dim(ERN.ihc)[1]  ){
    temp = ERP.ihc
    ERP.ihc = ERN.ihc
    ERN.ihc = temp
  }
  
  ## check ratio
  ## make sure a reasonable ratio
  if(ratio > dim(ERP.ihc)[1]/dim(ERN.ihc)[1] ) {stop( "please specify a ratio less than max(n_ER-, n_ER+ )/min(n_ER-, n_ER+)" ) }
  
  
  res_ihc_iterative = mapply(function (itr){
    
    ## get samples via ratio
    
    i = sample( dim(ERP.ihc)[1], ceiling( ratio * dim(ERN.ihc)[1]) )
    length(ERP.ihc$PatientID[i]) # ER positive samples
    length(ERN.ihc$PatientID)    # ER negative samples
    
    ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
    #mat = Test.matrix
    mbal.ihc = mat[,c(ERP.ihc$PatientID[i],ERN.ihc$PatientID)]
    
    dim(mbal.ihc) #[1]  50 60
    
    
    surffix = getsurffix(calibration = calibration, internal)
    
    # Calculate median
    mdns      = apply(mbal.ihc,1,median,na.rm=T) # compute median of each row i.e gene
    mdns.df   = as.data.frame(mdns)
    df.mdns   = data.frame(X=rownames(mdns.df),mdns.ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
    colnames(df.mdns) = c("X",surffix)
    
    ## integrate ihc.mdns
    fl.mdn = IBC.parker$medians
    
    df.al = merge(fl.mdn, df.mdns , by = "X")
    rownames(df.al) = df.al$X
    df.al = df.al[,-1]


    ## normalization
    mat = docalibration( mat, df.al, calibration,internal)
    
    out = sspPredict(centroids, classes="", mat, std=F, distm="spearman", centroids=T, hasClinical = hasClinical)
    
    return( out )
    
  },paste0("itr.", seq(iterative) ) ,SIMPLIFY = FALSE,USE.NAMES = TRUE)
  
  
  ## get consensus intrinsic subtype
  Call_subtypes = mapply(function(res_ihs){res_ihs$predictions }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
  consensus_subtypes = apply(Call_subtypes, 1, get_consensus_subtype)
  
  Int.sbs = data.frame(PatientID = df.cln$PatientID, Int.SBS.ihc.itr = consensus_subtypes, IHC = df.cln$IHC , row.names = NULL )


  ## get the entropy
  entropy_subtype  = apply(Call_subtypes, 1, get_entropy_subtype)

  ## get correlation and ROR score for each patient
  ## do average for each cell, matched with subtype
  ## how about correlation? average
  ## how about testdata? average
  mean_eve = get_average_subtype(res_ihc_iterative, consensus_subtypes, hasClinical = hasClinical)
  
  ## calculate ROR 
  ## group by ROR
  ## save all scores
  if(hasClinical ) {
    out = list(predictions = consensus_subtypes, testData = mean_eve$testdata,  distances = -1 * mean_eve$mean_distance, distances.prosigna = -1 * mean_eve$mean_distance.prosigna, centroids = centroids )
  }else{
    out = list(predictions = consensus_subtypes, testData = mean_eve$testdata,  distances = -1 * mean_eve$mean_distance, centroids = centroids )
  }
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical )
  
  out.ROR$Entropy = entropy_subtype
  out.ROR %<>% select( "Basal","Her2","LumA","LumB","Normal","Call","Confidence", "Entropy",everything() )
  
  ## each time, changing medians, hence,
  ## no saving medians in result (if want to save, do average then)
  ## in out, distances are average, subtypes are consensus subtype
  return( list(Int.sbs = Int.sbs, score.fl= out.ROR, outList = out, sbs.itr = Call_subtypes ))
  
}
  

#### form Secondary ER-balanced set (refer to paper) leveraing PCA and subsequent intermediate intrinsic subtype

#' Function for the first step of PCA-PAM50 approach
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd

makeCalls.PC1ihc = function(mat, df.cln, seed=118, calibration = "Internal", internal ="PC1ihc.mdns", external=NA, medians = NA , hasClinical = FALSE){
  # message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  # message("##IHC subtype column should be named ---IHC---")
  # 
  
  # ## test data
  # mat = Test.matrix
  # df.cln = df.cln
  # calibration = "Internal"
  # internal="PC1ihc.mdns"
  # seed=118
  # hasClinical = FALSE

  # # ## test data
  # mat = data_input$x_parker
  # df.cln = clinic
  # calibration = "Internal"
  # internal="PC1ihc.mdns"
  # seed=118
  # hasClinical = T


  #rv      = rowVars(mat)
  #select  = order(rv, decreasing = TRUE)[seq_len(dim(mat)[1])] # the input is PAM50 matrix --50 genes -- get from dimension
  pca     = prcomp(t(mat))#[select,]
  pc12    = pca$x[,1:2] #get two principal 
  df.pc1  = data.frame(PatientID=rownames(pc12),PC1 = pc12[,1],stringsAsFactors=F)
  df.pca1 = merge(df.cln,df.pc1,by="PatientID")
  
  # ## visualization
  # ## http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  # factoextra::fviz_eig(pca)
  # factoextra::fviz_pca_var(pca, col.var = "black")
  # factoextra::fviz_cos2(pca, choice = "var", axes = 1:2)
  # # Contributions of variables to PC1
  # factoextra::fviz_contrib(pca, choice = "var", axes = 1, top = 10)
  # # Contributions of variables to PC2
  # factoextra::fviz_contrib(pca, choice = "var", axes = 2, top = 10)
  # 
  # factoextra::fviz_pca_var(pca, col.var = "contrib",
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07") )
  # 
  # factoextra::fviz_pca_ind(pca, col.ind = "cos2", 
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #              repel = FALSE # Avoid text overlapping (slow if many points)
  # )
  # factoextra::fviz_pca_ind(pca,
  #              geom.ind = "point", # show points only (nbut not "text")
  #              col.ind = df.cln$IHC, # color by groups
  #              palette = c("#00AFBB", "#E7B800", "#FC4E07", "black","red"),
  #              addEllipses = FALSE, # Concentration ellipses
  #              legend.title = "Groups"
  # )
  # 
  # factoextra::fviz_cos2(pca, choice = "ind")
  # 
  # factoextra::fviz_contrib(pca, choice = "ind", axes = 1:2)
  # 
  ### function to count the number of mis-classified cases on a given PC1 point ---find the cutoff
  
  # getno = function(x){
  #   #x = 20
  #   p.rgt = length(which(df.pca1$IHC %in% c("LB1","LB2","LA") & df.pca1$PC1 >x))/length(which(df.pca1$IHC %in% c("LB1","LB2","LA") ))
  #   n.lft = length(which(df.pca1$IHC %in% c("Her2+","TN") & df.pca1$PC1 <x))/ length(which(df.pca1$IHC %in% c("Her2+","TN")))
  #   tot   = (p.rgt + n.lft) * 100
  #   return(list(PC1=x,Mis=tot))
  # }
  # 
  
  getno = function(x){
    p.rgt = length(which(df.pca1$ER == "ER+" & df.pca1$PC1 >x))/length(which(df.pca1$ER == "ER+"))
    n.lft = length(which(df.pca1$ER == "ER-" & df.pca1$PC1 <x))/ length(which(df.pca1$ER == "ER-"))
    tot   = (p.rgt + n.lft) * 100
    return(list(PC1=x,Mis=tot))
  }
  
  
  ## calculate the 
  df.mis  = do.call(rbind.data.frame,lapply(seq(-20,20,by=0.1),getno))
  
  ## ??? improve
  if(min(df.mis$Mis) < 100 ) {
    num.min = df.mis$PC1[which(df.mis$Mis == min(df.mis$Mis))]
  } else {
  num.min = df.mis$PC1[which(df.mis$Mis == max(df.mis$Mis))]
  }
  # 
  # ## visualizaton
  # ## Here, in train/example dataset, Her2, TN are on right side of PCA; LA LB1 LB2 are on left side
  # ## but in TCGA cell2015, it is on another way around. 
  # plot(x = df.mis$PC1, y =  df.mis$Mis)
  # 
  ERP.pc1ihc = df.pca1[which(df.pca1$ER == "ER+" & df.pca1$PC1 <= mean(num.min)),] # used mean to overcome situation where there are two minimum
  ERN.pc1ihc = df.pca1[which(df.pca1$ER == "ER-" & df.pca1$PC1 > mean(num.min)),]
  
  dim(ERP.pc1ihc)#  543   3
  dim(ERN.pc1ihc)#  118   3
  
  if(dim(ERP.pc1ihc)[1] < dim(ERN.pc1ihc)[1] ){
    temp = ERN.pc1ihc
    ERN.pc1ihc =  ERP.pc1ihc
    ERP.pc1ihc = temp
    rm(temp)
    }
  
  set.seed(seed);i = sample(dim(ERP.pc1ihc)[1],dim(ERN.pc1ihc)[1]) # take equal number of ER+ and ER- samples
  length(ERP.pc1ihc$PatientID[i]) # ER positive samples
  length(ERN.pc1ihc$PatientID)    # ER negative samples
  
  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pc1ihc = mat[,c(ERP.pc1ihc$PatientID[i],ERN.pc1ihc$PatientID)]
  
  dim(mbal.pc1ihc) #[1]  50 236
  
  # Find median
  
  surffix = getsurffix(calibration = calibration, internal)
  
  mdns      = apply(mbal.pc1ihc,1,median,na.rm=T) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.pc1ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
  colnames(df.mdns) = c("X",surffix)
  
  
  
  # # merge this median to PAM50-medians file "mediansPerDataset_v2.txt"
  # fl.nm     = paste(PAM50dir,"mediansPerDataset_v2.txt",sep="/")
  # fl.mdn    = read.table(fl.nm,sep="\t",header=T,stringsAsFactors=F)
  # 
  fl.mdn = IBC.parker$medians 
  
  df.al = merge(fl.mdn, df.mdns, by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  
  # normalization
  mat = docalibration( mat, df.al, calibration, internal)
  
  ## slow ??? yes
  out = sspPredict(IBC.parker$centroid, classes="", mat, std=F, distm="spearman", centroids=T, hasClinical = hasClinical)

  if(hasClinical) {
    out$distances.prosigna =  -1 * out$distances.prosigna
  }
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical )

  Int.sbs = data.frame(PatientID = names(out$predictions), Int.SBS.PC1ihc = out$predictions, IHC = df.cln$IHC , row.names = NULL )
  
  return(list(Int.sbs=Int.sbs, score.fl=out.ROR, mdns.fl= df.al,outList=out))
  
  
}

#' Function for the second step of PCA-PAM50 approach
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd

makeCalls.v1PAM = function(mat, df.pam, seed=118, calibration = "Internal", internal ="v1PAM.mdns", external=NA, medians = NA , hasClinical = FALSE ){
  
  # message("###df.pam data.frame should have a column --PatientID-- with which mat cols are also named")
  # message("##v1PAM subtype column should be named ---PAM50---")
  # 
  # print( paste0( "in VIPAM : ",internal ) )
  
  # mat = data_input$x_parker
  # df.pam = df.pc1pam
  # hasClinical = hasClinical
  # seed = 118
  # calibration = "Internal"
  # internal ="v1PAM.mdns"
  # 
  
  ERN.pam = df.pam[which(df.pam$PAM50 %in% c("Basal")),] ### get ER- samples data.frame
  dim(ERN.pam)
  
  ERP.pam = df.pam[which(df.pam$PAM50 %in% c("LumA")),]
  dim(ERP.pam)
  
  
  
  set.seed(seed);i = sample(dim(ERP.pam)[1],dim(ERN.pam)[1]) # take equal number of ER+ and ER- samples
  
  length(ERP.pam$PatientID[i]) # ER positive samples
  length(ERN.pam$PatientID)    # ER negative samples
  
  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pam = mat[,c(ERP.pam$PatientID[i],ERN.pam$PatientID)]
  
  dim(mbal.pam) #[1]  50 306
  
  # Find median
  surffix = getsurffix(calibration = calibration, internal)
  
  mdns      = apply(mbal.pam,1,median,na.rm=T) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.pam=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
  colnames(df.mdns) = c("X",surffix)
  
  
  # # merge this median to PAM50-medians file "mediansPerDataset_v2.txt"
  # fl.nm     = paste(PAM50dir,"mediansPerDataset_v2.txt",sep="/")
  # fl.mdn    = read.table(fl.nm,sep="\t",header=T,stringsAsFactors=F)
  # 
  # df.al     = merge(fl.mdn,df.mdns,by="X")
  # 
  
  fl.mdn = IBC.parker$medians

  df.al = merge(fl.mdn, df.mdns, by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  

  ## normalization
  mat= docalibration( mat, df.al, calibration,internal)
  
  out = sspPredict(IBC.parker$centroid, classes="", mat, std=F, distm="spearman", centroids=T, hasClinical = hasClinical)
  
  if(hasClinical) {
    out$distances.prosigna =  -1 * out$distances.prosigna
  }
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln = df.pam, hasClinical = hasClinical )
  
  Int.sbs = data.frame(PatientID = names(out$predictions), Int.SBS.ihc = out$predictions, IHC = df.pam$IHC , row.names = NULL )
  
  return(list(Int.sbs=Int.sbs, score.fl=out.ROR, mdns.fl= df.al,outList=out))
  
}



#' Function for calling PAM50 subtypes by ssBC methods
#' This function is adapted from ssBC TNBC-BreastCancerRes2015 and subgrou specific TNBC-JAMAOncol2024 
#' @param mat gene expression matrix, log of normalized
#' @param df.cln clicnical information table with PatientID and IHC column, it should include ER column with "ER-" or "ER+". Or it includes TN column with "TN".
#' @param s Specify "ER" or "TN" or "ER_JAMA" or "HER2+" or "TNBC". The original quantile is "ER" and "TN" of TNBC-BreastCancerRes2015.  If you choose "ER_JAMA" or "HER2+" or "TNBC", it means you choose quantile from TNBC-JAMAOncol2024. 
#' @param hasClinical Logic. 
#' @noRd
#' 

## add one more module sub-group specific BC 
## https://github.com/afernan4/PAM50_ER_HER2_normalization/blob/main/SIGMA.txt
## https://ascopubs.org/doi/suppl/10.1200/JCO.20.01276/suppl_file/ds_jco.20.01276-2.pdf

makeCalls.ssBC = function(mat, df.cln, s , hasClinical =FALSE  ){
  
  # ## test data
  # data("IBC.parker")
  # gene.sigma = IBC.parker$ssBC.subgroupQuantile
  # mat = exprs(UNC232)
  # df.cln= pData(UNC232)
  # df.cln$IHC = ""
  # 
  # s = "ER"
  
  # ## test data
  # mat = data_input$x_parker
  # df.cln = clinic.scanb
  # s = "ER"
  # hasClinical =T

  gene.sigma = IBC.parker$ssBC.subgroupQuantile
  

  if( s == "ER") { ## use ER selected strategy
    
    ## if there is no sample in either of both, wont influence the code
    ERN_samples = rownames(df.cln)[which(df.cln$ER == "ER-" )]
    ERP_samples = rownames(df.cln)[which(df.cln$ER == "ER+" )]
    samples_selected = list( ER_neg = ERN_samples, ER_pos = ERP_samples)
    
  } else if(s == "TN"){
    
    ## if there is no sample in either of both, wont influence the code
    TN_samples = rownames(df.cln)[which(df.cln$TN == "TN" )]
    samples_selected = list( TN = TN_samples)
    
  }  else if( s == "ER_JAMA" ){
    
    ## if there is no sample in either of both, wont influence the code
    ERN_HER2N_samples = rownames(df.cln)[which(df.cln$ER == "ER-" & df.cln$HER2 =="HER2-"  )]
    ERP_HER2N_samples = rownames(df.cln)[which(df.cln$ER == "ER+" & df.cln$HER2 =="HER2-" )]
    
    ERN_HER2P_samples = rownames(df.cln)[which(df.cln$ER == "ER-" & df.cln$HER2 =="HER2+"  )]
    ERP_HER2P_samples = rownames(df.cln)[which(df.cln$ER == "ER+" & df.cln$HER2 =="HER2+"  )]
    
    
    samples_selected = list( ERneg_HER2neg = ERN_HER2N_samples, ERpos_HER2neg = ERP_HER2N_samples,
                             HER2pos_ERneg = ERN_HER2P_samples, HER2pos_ERpos = ERP_HER2P_samples)
    
    
  } else if ( s == "TNBC") { ## selected cohort; TNBC-JAMAOncol2024
    
    ## if there is no sample in either of both, wont influence the code
    TN_samples = rownames(df.cln)[which(df.cln$TN == "TN" )]
    #nTN_samples = rownames(df.cln)[which(df.cln$TN == "nTN" )] 
    samples_selected = list( TNBC = TN_samples)
    
  } else {
    stop("Please enter valid varaible for s ")
  }
  
  res = mapply(function(element){
    
    #print(element)
    
    x.m = mat[, samples_selected[[element]] ]
    
    ## ensure gene signatures are used in centroid table
    ## to be improve later

    ## do calibration
    ## the same as the python 
    ## https://github.com/afernan4/PAM50_ER_HER2_normalization/"PAM50 gene normalization by clinically defined breast cancer subgroups.ipynb"
    gene.sigma.o = gene.sigma[rownames(x.m), element]
    x.sigma = unlist(lapply(1:nrow(x.m), function(i) quantile(x.m[i,], probs = gene.sigma.o[i], na.rm = T)))
    x.m = sweep(x.m, 1, x.sigma) ## subtract
    ## it has been calibrated by selected quantile 
    
    
    out = sspPredict(IBC.parker$centroid, classes="", x.m , std=F, distm="spearman", centroids=T, hasClinical = hasClinical)
    
    
  
  }, names(samples_selected), SIMPLIFY = F, USE.NAMES = T )
  
  ## prepare data for ROR

  if( s == "ER") { ## use ER selected strategy
    
    prediction = c(res$ER_neg$predictions, res$ER_pos$predictions)
    distances = rbind(res$ER_neg$distances, res$ER_pos$distances )
    testData = cbind(res$ER_neg$testData, res$ER_pos$testData )
    
    if(hasClinical ){
      distances.prosigna = rbind(res$ER_neg$distances.prosigna, res$ER_pos$distances.prosigna )
    }
    
  } else if(s == "ER_JAMA" ){
    
    prediction = c(res$ERneg_HER2neg$predictions, res$ERpos_HER2neg$predictions, 
                   res$HER2pos_ERneg$predictions, res$HER2pos_ERpos$predictions)
    distances = rbind(res$ERneg_HER2neg$distances, res$ERpos_HER2neg$distances,
                      res$HER2pos_ERneg$distances, res$HER2pos_ERpos$distances)
    testData = cbind(res$ERneg_HER2neg$testData, res$ERpos_HER2neg$testData,
                     res$HER2pos_ERneg$testData, res$HER2pos_ERpos$testData)
    
    if(hasClinical ){
      distances.prosigna = rbind(res$ERneg_HER2neg$distances.prosigna, res$ERpos_HER2neg$distances.prosigna,
                                 res$HER2pos_ERneg$distances.prosigna, res$HER2pos_ERpos$distances.prosigna)
    }
    
    
  } else if(s == "TN"){
    
    prediction = res$TN$predictions
    distances = res$TN$distances
    testData = res$TN$testData

    if(hasClinical ){
      distances.prosigna = res$TN$distances.prosigna
    }
    
  } else if (s == "TNBC")  {
    
    prediction = res$TNBC$predictions
    distances = res$TNBC$distances
    testData = res$TNBC$testData
    
    if(hasClinical ){
      distances.prosigna = res$TNBC$distances.prosigna
    }
  } 
  
  
  if ( !hasClinical) {
    out = list(predictions=prediction,testData=testData,distances=distances,centroids= IBC.parker$centroid)
  }else {
    out= list(predictions=prediction,testData=testData,distances=distances,distances.prosigna = distances.prosigna,centroids=IBC.parker$centroid)
  }
  
  
  if(hasClinical) {
    out$distances.prosigna =  -1 * out$distances.prosigna
  }
  out$distances = -1 * out$distances
  
  
  ##calculate and grouping
  out.ROR = RORgroup(out,df.cln, hasClinical = hasClinical )
  
  Int.sbs = data.frame(PatientID = names(out$predictions), Int.SBS.ihc = out$predictions, IHC = df.cln[ which( rownames(df.cln) %in% names(out$predictions) ) ,"IHC"], row.names = NULL )
  
  return(list(Int.sbs=Int.sbs, score.fl=out.ROR, mdns.fl= gene.sigma,outList=out))
  
}
  

