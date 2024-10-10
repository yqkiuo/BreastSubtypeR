#' 
#' Functions adapted from original parker-based PAM50 subtyping
#' @name PCA50
#' @import ggplot2
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import circlize
#' @import magrittr
#' @import factoextra
#' @import impute
#' @importFrom dplyr select
#' @importFrom dplyr mutate_at
NULL

#' 
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
#' @param x Gene expression matrix or vector
#' @noRd 
rescale <- function(x, na.rm=FALSE, q=0) {
  # x = mat[1,]
  # q = 0.05
  # na.rm =TRUE
  
  if(q == 0) {
      ma <- max(x, na.rm=na.rm)
      mi <- min(x, na.rm=na.rm)
    } else {
      ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
      mi <- quantile(x, probs=q/2, na.rm=na.rm)
    }
    xx <- (x - mi) / (ma - mi)
    #attributes(xx) <- list("names"=names(x), "q1"=mi,"q2"=ma)
    return(xx)
}


#' function for calibration methods
#' @param y Gene expression matrix 
#' @param df.al Medians for calibration
#' @param calibration How to do calibration, "None"(default) means no calibration for gene expression matrix. When setting calibration =None, you dont need to set internal and external parameters.  "Internal" means calibration for gene expression matrix by itself. "External" means calibration by external cohort. 
#' @param internal Specify the strategy for internal calibration, medianCtr(default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians calculated by train cohorts. When users want to use Medians prepared by user selves, this parameter should be "Given.mdns", not platform name. 
#' @noRd 
docalibration = function( y, df.al,calibration = "None", internal=internal, external=external){
  # y = mat
  # df.al = df.al
  # internal = "meanCtr"

  mq = 0.05 ## presetting in genefu robust model
  switch( calibration,
          "None" = {print("No calibration") }, # genefu none
          "Internal" = { ## internal
            if(internal == "medianCtr"){ y = medianCtr(y) } ## parker default method
            else if(internal == "meanCtr") { y =  t(scale( t(y), center=TRUE, scale=TRUE))} ## "scale" in genefu method
            else if(internal == "qCtr") { y = t(apply(y, 1, function(x) {  return( (rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2 ) }) )  } ## "robust" in genefu method; mp = 0.05
            else if(internal == internal) { ## which column to use
              medians =  readarray(df.al)
              #print(paste("calibration to:",internal))
              tm = overlapSets(medians$xd,y)
              y = (tm$y-tm$x[,internal])
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
#' @param x Gene expression matrix. Probe or gene symbol in row and sample in column.
#' @param y Feature data provided by user. The table should contain at least three column, which are probe(probeid or transcriptID), EntrezGene.ID and symbol. 
#' @param method Method to deduplicated probes for microarray or RNAseq. Please select "IQR" for Affy, "mean" for Agilent and "max" for RNAseq. 
#' @param impute Logic. Please specify whether to perform impute.knn on NA data
#' @param verbose Logic. 
#' @noRd

#### wrap the code 
#### use probeid as input in feature table

domapping = function(x, y = NA, method = "max", impute = TRUE, verbose = TRUE ){
  
  # x = expr
  # y = anno_feature #probe ENTREZID
  # method="max"
  # mapping= TRUE
  # impute = TRUE
  # verbose = TRUE
  # 
  ## loading genes.signature
  genes.signature = BreastSubtypeR$genes.signature
  
  ## first step 
  ## for empty cells. imput or not ?
  if(sum(apply(x,2,is.na))>0 & impute){
    
    if(verbose){
      probeid_NA = rownames(x)[rowSums(is.na(x)) > 0]
      sample_NA = colnames(x)[colSums(is.na(x)) > 0]
      print(paste0("The imput objects: ", probeid_NA, " in ", sample_NA))
    }
    
    x = impute.knn(x)
    x = x$data
  }
  
  ## no provided y and !mapping(RNAseq)
  if(length(y) == 0 & !mapping ){
    y =AnnotationDbi::select(org.Hs.eg.db, keys =rownames(x), columns = c( "ENTREZID","SYMBOL"), keytype='SYMBOL' ) 
  } else if(length(y) == 0 & mapping) {
    print("Please provide feature annotation to map probeID or transcripID ")
  }
  
  ## deduplicated
  probeid = rownames(x) ## probes, ensembl, symbol
  entrezid = y$ENTREZID
  names(entrezid) = y$probe
  
  ## process probeid in input data
  entrezid = entrezid[probeid]
  ##remove NA
  entrezid = entrezid[!(is.na(entrezid))]
  x = x[names(entrezid),]
  entrezid = factor( entrezid, levels =  unique(entrezid) ) 
  ## names are unique probeid and content are redundant entrezid 
  
  
  ## This is for probeID or transcriptID or symbol
  ## split expression matrix
  split_mat = split( data.frame(x), entrezid, drop = F)
  
  # function to calculate the desired statistic
  calculate_stat = function(mat, method) {
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
           },
           "max" = mat[ which.max( rowSums(mat) ),]
    )
  }
  
  ## keep processed x
  x = mapply( calculate_stat, split_mat, MoreArgs = list(method = method), SIMPLIFY = T, USE.NAMES = T)
  x = apply(x, 1, unlist) 
  
  ##print necessary information
  ##Parker
  missing_ID_parker = setdiff( BreastSubtypeR$genes.sig50$EntrezGene.ID, rownames(x) )
  if( length(missing_ID_parker) == 0 & verbose ){ 
    print("Genes used in NC-based methods are covered")
  } else {
    print("These genes are missing for NC-based methods:")
    print(missing_ID_parker)
  }
  
  ##AIMS
  missing_ID_AIMS = setdiff( genes.signature[ genes.signature$SSP_based == "Yes",]$EntrezGene.ID, rownames(x) )
  if( length(missing_ID_AIMS) == 0 & verbose ){ 
    print("Genes used in SSP-based methods are covered")
  } else  {
    print("These genes are missing for SSP-based methods:")
    print(missing_ID_AIMS)
  }
  
  
  ## get matrix for NC (symbol as rows, sample as col)
  x_NC = x[ match( BreastSubtypeR$genes.sig50$EntrezGene.ID, rownames(x) ) ,]
  rownames(x_NC) = BreastSubtypeR$genes.sig50$Symbol[ match( rownames(x_NC), BreastSubtypeR$genes.sig50$EntrezGene.ID ) ]

  ## get matrix for AIMS (entrezID as colnames)
  x_SSP = x[ match(as.character( genes.signature[ genes.signature$SSP_based == "Yes",]$EntrezGene.ID) , rownames(x)  ),]
  
  result = list(x_NC = data.frame(x_NC), x_NC.log = log2( data.frame(x_NC) +1) , x_SSP =data.frame( x_SSP) )
  
  return(result)
  
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
#' @param y Gene expression matrix
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
#' Function for de-duplicated genes in medians of train cohort
#' @param x Median of train cohort
#' @param method mean, median, stdev, iqr
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

#' Function for data structure of medians of train cohort
#' @param data Median of train cohort
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



#' Function to get the average correlation and ROR
get_average_subtype = function(res_ihc_iterative, consensus_subtypes) {

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
    
    res = mutate_at(res_ihc$distances, vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
    
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

 
  ## distances.prosigna for ROR 
  sum_cols_list.prosigna = mapply(function(res_ihc){
    
    #res_ihc = res_ihc_iterative[[1]]
    
    res_ihc$distances.prosigna = as.data.frame(res_ihc$distances.prosigna )
    
    ## if FALSE, make the cell as NULL
    keep = res_ihc$predictions == consensus_subtypes
    res_ihc$distances.prosigna[ !keep, ] = as.list(rep(NA, 4 ))
    
    res = mutate_at(res_ihc$distances.prosigna, vars( everything() ), ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA))
    
    return(res )
    
  }, res_ihc_iterative , SIMPLIFY = FALSE, USE.NAMES = FALSE )
  
  ## count_na
  count_weight_save.prosigna <- Reduce(`+`, lapply(sum_cols_list.prosigna, function(x) {
    x[!is.na(x)] <- 1
    x[is.na(x)] <- 0
    return(x)
  }))
  
  
  ## sum all for each cell
  sum_cols_save.prosigna = Reduce(`+`, lapply(sum_cols_list.prosigna, function(x) {
    
    ## change NA cell to 0 cell
    x[is.na(x)] = 0
    
    return(x)
  }))
  
  ## get the mean for each cell
  ## only when subtype is supported by consensus_subtypes for each iteration and each patient
  sum_cols_save.prosigna = sum_cols_save.prosigna / count_weight_save.prosigna
  
  
  res = list( mean_distance = mean_cols_save, mean_distance.prosigna = sum_cols_save.prosigna, testdata =  mean_cols_save.testdata) 
  

  return(res)
}

#' Function for predict PAM50 subtyping
#' @param x median train file
#' @param y gene expression matrix
#' @param classes description
#' @param nGenes None
#' @param distm "euclidean" or "spearman"
#' @param Prosigna Logic. Please specify if it predicts prosigna-like subtype
#' @param std Logical value. 
#' @param centrids Logical value
#' @noRd
sspPredict<-function(x, y, std=FALSE, distm="euclidean",centroids=FALSE, Prosigna = TRUE){
  
  # ## test data
  # x = BreastSubtypeR$centroid
  # y = data_input$x_parker
  # classes=""
  # nGenes=""
  # std=F
  # distm="spearman"
  # centroids=T
  # Prosigna=T

  dataMatrix = x
  tdataMatrix = y

  #dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
  temp = overlapSets(dataMatrix,tdataMatrix)
  dataMatrix = temp$x
  tdataMatrix = temp$y
  sfeatureNames = row.names(dataMatrix)

  # standardize both sets
  if(std){
    dataMatrix<-standardize(dataMatrix)
    tdataMatrix<-standardize(tdataMatrix)
  }

  nGenes<-dim(dataMatrix)[1]
  #print(paste("Number of genes used:",nGenes))
  centroids<-dataMatrix
  nClasses<-dim(centroids)[2] ## five subtypes; for prosigna, keep normal when calculating
  classLevels<-dimnames(centroids)[[2]]

  distances = matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
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
  names(prediction) = colnames(tdataMatrix)
  
  ## run Prosigna like subtype
  if( Prosigna){
    
    nClasses = nClasses-1 ## omitting normal
    classLevels = classLevels[1:4]
    
    distances.prosigna.subtype = matrix(ncol= nClasses, nrow=dim(tdataMatrix)[2])
    for(j in 1:nClasses){
      if(distm=="euclidean"){
        distances.prosigna.subtype[,j]<- dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
      }
      if(distm=="correlation" | distm=="pearson"){
        distances.prosigna.subtype[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids[,j], x, method = "pearson", use = "pairwise.complete.obs"))
      }
      if(distm=="spearman"){
        distances.prosigna.subtype[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids[,j], x, method = "spearman", use = "pairwise.complete.obs"))
      }
    }
    
    prediction.prosigna = classLevels[apply(distances.prosigna.subtype, 1, which.min,simplify = TRUE)]
    names(prediction.prosigna) =  colnames(tdataMatrix)

  }
  
  ## run Prosigna like ROR
  ## prepare prosigna.distances

  ## genes to be excluded
  genes.ex = c("BIRC5", "CCNB1", "GRB7","MYBL2")
  
  #parse the test file - same as train file but no rows of classes
  tdataMatrix = y[ which( !(rownames(y) %in%  genes.ex ) ),]

  #dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
  temp = overlapSets(dataMatrix,tdataMatrix)
  dataMatrix = temp$x
  tdataMatrix = temp$y


  ## omitting normal
  dataMatrix = dataMatrix[,1:4] ## omitting normal or not for ROR???
  nGenes = dim(dataMatrix)[1]
  centroids.prosigna = dataMatrix
  nClasses =dim(centroids.prosigna)[2] ## four subtypes
  classLevels = dimnames(centroids.prosigna)[[2]]

  distances.prosigna  =  matrix(ncol= nClasses,nrow=dim(tdataMatrix)[2])
  for(j in 1:nClasses){
    if(distm=="euclidean"){
      distances.prosigna[,j] = dist(t(cbind(centroids.prosigna[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
    }
    if(distm=="correlation" | distm=="pearson"){
      distances.prosigna[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids.prosigna[,j], x, method = "pearson", use = "pairwise.complete.obs"))
      }
    if(distm=="spearman"){
      distances.prosigna[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids.prosigna[,j], x, method = "spearman", use = "pairwise.complete.obs"))
      }
  }

  

  ## return 
  if ( Prosigna) {
    res = list(predictions=prediction, predictions.prosigna = prediction.prosigna,testData= as.matrix(y),distances=distances,distances.prosigna = distances.prosigna, distances.prosigna.subtype = distances.prosigna.subtype, centroids=centroids)
  
  }else {
    res = list(predictions=prediction,testData= as.matrix(y), distances=distances, distances.prosigna = distances.prosigna,centroids=centroids)
  }
  
  return(res)
}


#' Function for risk 
#' 
#' @param hasClinical provide clinical information for grouping 
#' @param out it is the result of sspPredict() function. 
#' @return ROR, ROR risk group and other indications
#' @noRd

RORgroup = function(out, df.cln , hasClinical = FALSE, Prosigna = FALSE ){
  
  # # ## test data
  # out
  # df.cln = df.cln
  # hasClinical = hasClinical


  Clinical = df.cln[match(names(out$predictions),df.cln$PatientID ),]
  rownames(Clinical) = Clinical$PatientID
  ## prepare outtable
  # 
  # 
  
  distance = data.frame(out$distances, row.names = names(out$predictions) )
  colnames(distance) = c("Basal","Her2","LumA","LumB","Normal")

  Call = data.frame( "Call" = out$predictions, row.names = names(out$predictions))
  
  #providing proliferation signatures
  proliferationGenes = c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
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
  call.conf = c()
  for(j in 1:length(out$predictions)){
    call.conf[j]= 1-cor.test(out$testData[,j],out$centroids[,which(colnames(out$centroids)==out$predictions[j])],method="spearman", exact=FALSE)$p.value
  }
  
  call.conf = round(call.conf,2)
  call.conf = data.frame( "Confidence" = call.conf, row.names = names(out$predictions))
  

  ## calculate the risk scores
  ## genomic
  ## basal, her2, lumA, lumB in order
  genomic = 0.04210193*out$distances[,1] + 0.12466938*out$distances[,2] + -0.35235561*out$distances[,3] + 0.14213283*out$distances[,4]
  
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
  
  if (hasClinical){
    
    if ( "T" %in% colnames(Clinical)) {
      
      xT= as.numeric(as.vector(Clinical$T))
      
      combined = 0.0442770*out$distances[,1] + 0.1170297*out$distances[,2] + -0.2608388*out$distances[,3] + 0.1055908*out$distances[,4] + 
        0.1813751*xT
      combinedWprolif = -0.009383416*out$distances[,1] +  0.073725503*out$distances[,2] + -0.090436516*out$distances[,3] + 0.053013865*out$distances[,4] + 
        0.131605960*xT + 0.327259375*prolifScore
      
      criskgroups = combined
      criskgroups[combined>chthreshold] = "high"
      criskgroups[combined>clthreshold & combined<chthreshold] = "med"
      criskgroups[combined<clthreshold] = "low"
      cpriskgroups = combinedWprolif
      cpriskgroups[combinedWprolif>cphthreshold] = "high"
      cpriskgroups[combinedWprolif>cplthreshold & combinedWprolif<cphthreshold] = "med"
      cpriskgroups[combinedWprolif<cplthreshold] = "low"
      
      combined =  100* (combined + 0.35 ) / 0.85
      combinedWprolif =  100* (combinedWprolif + 0.35 ) / 0.85
      
      ROR.combined = data.frame("ROR-C (Subtype + Clinical)" = combined, "ROR-C Group (Subtype + Clinical)" = cpriskgroups,
                                "ROR-PC (Subtype + Clinical + Proliferation)" = combinedWprolif, "ROR-PC Group (Subtype + Clinical + Proliferation)" = cpriskgroups,
                                check.names = FALSE)
      

      #View(ROR.combined)
      
      # ROR prosigna
      ## log2 of nCounter expression, FFPE
      ## log2 of FPKM, FF , illumina
      ## 46 genes, 
      ## 46 genes to calculate correlation or no
 
      combinedWprolif.prosigna = 54.7690 * (-0.0067 * out$distances.prosigna[,1] + 0.4317*out$distances.prosigna[,2] - 0.3172*out$distances.prosigna[,3] + 
                                              0.4894*out$distances.prosigna[,4] + 0.1981*prolifScore.prosigna + 0.1133*xT + 0.8826)
    

      if ("NODE" %in% colnames(Clinical)  ) {
        
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
            cpriskgroups.prosigna[patients.NODE.3][ temp> cphthreshold.prosigna.NODE ]  = "high"
            cpriskgroups.prosigna[patients.NODE.3][temp>cplthreshold.prosigna.NODE & temp<cphthreshold.prosigna.NODE ] = "med"
            cpriskgroups.prosigna[patients.NODE.3][temp<cplthreshold.prosigna.NODE] = "low"
            
          }
          
          ## negative node
          patients.NODE.0 = rownames(Clinical)[ Clinical$NODE == 0 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.0 ) > 0  ) { 
            
            temp = combinedWprolif.prosigna[patients.NODE.0]
            cpriskgroups.prosigna[patients.NODE.0][ temp> cphthreshold.prosigna.NODE.0 ]  = "high"
            cpriskgroups.prosigna[patients.NODE.0][temp>cplthreshold.prosigna.NODE.0 & temp<cphthreshold.prosigna.NODE.0 ] = "med"
            cpriskgroups.prosigna[patients.NODE.0][temp<cplthreshold.prosigna.NODE.0] = "low"
            
          }
          
          
          ## check if NA NODE 
          ## grouping by NA NODE
          patients.NA = rownames(Clinical)[ is.na(Clinical$NODE ) | is.na(Clinical$T) ]
          
          if ( length( patients.NA ) > 0  ){
            
            cpriskgroups.prosigna[patients.NA] = NA
            
          }
          
          ROR.combined.prosigna = data.frame("ROR-PC (Subtype + Clinical + Proliferation.prosigna)" = combinedWprolif.prosigna, 
                                             "ROR-PC Group (Subtype + Clinical + Proliferation.prosigna)" = cpriskgroups.prosigna,
                                             check.names = FALSE)
        
   
      } else {
        
        cpriskgroups.prosigna[combinedWprolif.prosigna > cphthreshold.prosigna] = "high"
        cpriskgroups.prosigna[combinedWprolif.prosigna > cplthreshold.prosigna & combinedWprolif.prosigna < cphthreshold.prosigna] = "med"
        cpriskgroups.prosigna[combinedWprolif.prosigna < cplthreshold.prosigna] = "low"
        

        ROR.combined.prosigna = data.frame("ROR-PC (Subtype + Clinical + Proliferation.prosigna)" = combinedWprolif.prosigna, 
                                           "ROR-PC Group (Subtype + Clinical + Proliferation.prosigna)" = cpriskgroups.prosigna,
                                           check.names = FALSE)
        
      }
      
      } 
      
      }
    

    outtable = cbind( distance, Call, call.conf, ROR.genomic,ROR.combined, ROR.combined.prosigna, er_her2)
    
  } else {
    outtable = cbind( distance, Call, call.conf, ROR.genomic, er_her2)
  }
  
  
  if(Prosigna){
    
    ## pass distances (just omitting normal)
    Call.prosigna = data.frame( "Call.prosigna" =out$predictions.prosigna, row.names = names(out$predictions.prosigna))
    outtable = cbind(outtable, Call.prosigna)
    outtable =outtable %>% dplyr::select(colnames(distance), Call, Call.prosigna, everything() )
    
  }

  return(outtable)
}



#' Functions for visualization
#' Function for boxplot of correlation per subtype
#' @param out a data frame includes "patientID" and "Subtype"
#' @param correlations  correlations table from NC-based methods
#' @export
#' 

Vis_boxpot = function(out, correlations ){
  
  df = data.frame( predictions = out$Subtype, cor = apply(correlations, 1, max))
  
  plot =  ggplot( df, aes( x = predictions, y = cor) ) +
    geom_boxplot()
  
 return(plot)
  
}


#' Function for heatmap visualizayion
#' @param x gene expression matrix, log2 transformed
#' @param out a data frame includes "patientID" and "Subtype"
#' @export
#' 

Vis_heatmap = function(x, out){

  # ## test data
  # x = data_input$x_NC
  # out= data.frame(PatientID = res$results$parker.median$BS.all$PatientID,
  #                 Subtype = res$results$parker.median$BS.all$BS )
  # 
  
  scaled_mat = t(scale(t(x)))
  
  col_fun = colorRamp2(c(min(scaled_mat), 0 , max(scaled_mat)), c("green", "black", "red"))

  ## column annotation
  col_anno = data.frame( row.names = out$PatientID, Subtype = out$Subtype )

  anno_col = HeatmapAnnotation(df= col_anno, show_legend = FALSE,col = list(Subtype = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" ) ))
  
  heatmap = Heatmap(scaled_mat, name = "Subtype",
          col = col_fun,
          ## annotation
          top_annotation = anno_col,
          
          ## clustering
          ## as original heatmap plot
          clustering_distance_rows = "pearson",
          clustering_method_rows = "complete",
          
          cluster_column_slices = TRUE,
          column_split = col_anno$Subtype,
          clustering_distance_columns  = "pearson",
          clustering_method_columns = "complete",

          ## general
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          row_names_gp = gpar(fontsize = 8))
  
  return(heatmap)
}
  



#' Function for PCAplot 
#' @param x gene expression matrix, log2 transformed
#' @param out a data table includes "patientID" and "Subtype"
#' @param screeplot Logic. Please specify if show screeplot
#' @export
#' 

## reduce dependency ???
Vis_PCA = function(x, out, Eigen = FALSE){
  
  
  Subtype.color = c( "Basal" = "red", "Her2" = "hotpink","LumA" = "darkblue", "LumB" = "skyblue" , "Normal" = "green" )
  
  #x = data_input$x_parker
  
  x_pca = prcomp(t(x),scale. = TRUE )
  
  screeplot = fviz_eig(x_pca, addlabels = TRUE, ylim = c(0, 50))
  
  pcaplot = fviz_pca_ind(x_pca, label="none",mean.point = FALSE, pointshape = 16 ,
               col.ind = as.factor(out$Subtype )) + 
    scale_color_manual( name = "Subtype", values = Subtype.color )
 
  if(Eigen){
    return(screeplot)
  } else {
    return(pcaplot)
  }

  
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

makeCalls.parker = function(mat, df.cln, calibration = "None", internal = NA,external=NA, medians = NA,Prosigna = FALSE, hasClinical =FALSE  ){
  
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

  # mat = data_input$x_NC.log
  # df.cln = clinic.oslo
  # calibration = "Internal"
  # internal = "qCtr"
  # external= NA
  # medians = NA
  # Prosigna = T
  # hasClinical = T



  
  fl.mdn = BreastSubtypeR$medians

  
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
  
  out = sspPredict(BreastSubtypeR$centroid, mat, std=F, distm="spearman", centroids=T, Prosigna = Prosigna)
  

  if (Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna , IHC = df.cln$IHC , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS= out$predictions, IHC = df.cln$IHC , row.names = NULL )
  }
  out$distances.prosigna =  -1 * out$distances.prosigna
  out$distances = -1 * out$distances
  

  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical )

  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns = df.al, outList=out))
  
  
}



#### function to form a ER-balance subet and derive its median---write it to PAM50 dir

#' Function to form a ER-balance subet and derive its median
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd
makeCalls.ihc = function(mat, df.cln, seed=118,calibration = "Internal", internal = "IHC.mdns", external=NA, medians = NA , Prosigna = FALSE , hasClinical = FALSE){
  # message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  # message("##IHC subtype column should be named ---IHC---")
  # 
  #test data for TCGA Cell2015
  # mat = data_input$x_NC.log
  # df.cln = clinic
  # seed=118
  # calibration = "Internal"
  # internal = "IHC.mdns"
  # hasClinical = TRUE
  # Prosigna = TRUE
  
  ERN.ihc = df.cln[which(df.cln$ER == "ER-"),] ### get ER- samples data.frame
  dim(ERN.ihc)	#[1] 153   9
  
  ERP.ihc = df.cln[which(df.cln$ER == "ER+"),]
  dim(ERP.ihc) #[1] 559   9
  
  #seed = 118
  if( dim(ERN.ihc)[1] > dim(ERP.ihc)[1] ) {  temp = ERP.ihc; ERP.ihc = ERN.ihc; ERN.ihc = temp }
  set.seed(seed); i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1]) # take equal number of ER+ and ER- samples
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
  fl.mdn = BreastSubtypeR$medians
  
  df.al = merge(fl.mdn, df.mdns , by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  ## centroids
  centroids = BreastSubtypeR$centroid #pam50_centroids.txt
  

  ## normalization
  mat = docalibration( mat, df.al, calibration, internal)
  
  out = sspPredict( centroids, mat, std=F, distm="spearman", centroids=T, Prosigna = Prosigna)
  
  if (Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna ,IHC = df.cln$IHC , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, IHC = df.cln$IHC , row.names = NULL )
  }
  out$distances.prosigna =  -1 * out$distances.prosigna
  out$distances = -1 * out$distances

  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical, Prosigna = Prosigna)

  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns= df.al, outList=out))
  
}



#' Function for iterative ER subset gene centering 
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param iteration times to predict subtypes
#' @param ratio The options are either 1:1 or 54(ER+):64(ER-). The latter is ER ratio used for UNC230 train cohort
#' @param calibration surffix name for output
#' @param hasClinical provide clinical information, default is NA. 
#' @noRd

makeCalls.ihc.iterative = function( mat, df.cln, iteration = 100, ratio = 54/64, calibration = "Internal", internal = "ER.mdns", external=NA, medians = NA , Prosigna = FALSE, hasClinical = FALSE, seed=118){

  # # ## test data
  # mat = data_input$x_NC.log
  # df.cln = clinic.oslo
  # ratio=1
  # ratio= 54/64
  # ##  46% ER-positive (54/118) and 54% ER-negative (64/118) from Zhao paper ,
  # calibration = "Internal"
  # internal = "ER.mdns"
  # iterative = 5
  # hasClinical = T
  # Prosigna=F
  # seed = 118
  # iteration = 100
  # 
  # set.seed(seed)
  
  
  # load the published centroids for classifcation
  centroids = BreastSubtypeR$centroid #pam50_centroids.txt
  
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
    fl.mdn = BreastSubtypeR$medians
    
    df.al = merge(fl.mdn, df.mdns , by = "X")
    rownames(df.al) = df.al$X
    df.al = df.al[,-1]


    ## normalization
    mat = docalibration( mat, df.al, calibration,internal)
    
    out = sspPredict(centroids, mat, std=F, distm="spearman", centroids=T, Prosigna = Prosigna)
    
    return( out )
    
  },paste0("itr.", seq(iteration) ) ,SIMPLIFY = FALSE,USE.NAMES = TRUE)
  
  
  ## get consensus intrinsic subtype
  Call_subtypes = mapply(function(res_ihs){res_ihs$predictions }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
  consensus_subtypes = apply(Call_subtypes, 1, get_consensus_subtype)
 
  if (Prosigna) {
    Call_subtypes.prosigna = mapply(function(res_ihs){res_ihs$predictions.prosigna }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
    consensus_subtypes.prosigna = apply(Call_subtypes.prosigna, 1, get_consensus_subtype)
    
    Int.sbs = data.frame(PatientID = names(consensus_subtypes), BS = consensus_subtypes, BS.prosigna = consensus_subtypes.prosigna , IHC = df.cln$IHC , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID =names(consensus_subtypes), BS = consensus_subtypes, IHC = df.cln$IHC , row.names = NULL )
  }
  
  ## get correlation and ROR score for each patient
  ## do average for each cell, matched with subtype
  ## how about correlation? average
  ## how about testdata? average
  mean_eve = get_average_subtype(res_ihc_iterative, consensus_subtypes)
  
  ## calculate ROR 
  ## group by ROR
  ## save all scores
  
  if (Prosigna) {
    out = list(predictions = consensus_subtypes, predictions.prosigna = consensus_subtypes.prosigna, testData = mean_eve$testdata, distances = -1 * mean_eve$mean_distance, distances.prosigna = -1 * mean_eve$mean_distance.prosigna, centroids = centroids )
  } else {
    out = list(predictions = consensus_subtypes, testData = mean_eve$testdata,  distances = -1 * mean_eve$mean_distance,  distances.prosigna = -1 * mean_eve$mean_distance.prosigna, centroids = centroids )
  }

  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical,Prosigna = Prosigna )
  
  ## each time, changing medians, hence,
  ## no saving medians in result (if want to save, do average then)
  ## in out, distances are average, subtypes are consensus subtype
  if(Prosigna ) {
    res =  list(BS.all = Int.sbs, score.ROR= out.ROR, outList = out, BS.itr.keep = Call_subtypes, BS.itr.keep = Call_subtypes.prosigna )
  } else {
    res = list(BS.all = Int.sbs, score.ROR= out.ROR, outList = out, BS.itr.keep = Call_subtypes )
  }
  
  return(res)
}
  

#### form Secondary ER-balanced set (refer to paper) leveraing PCA and subsequent intermediate intrinsic subtype

#' Function for the first step of PCA-PAM50 approach
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd

makeCalls.PC1ihc = function(mat, df.cln, seed=118, calibration = "Internal", internal ="PC1ihc.mdns", external=NA, medians = NA ,Prosigna =FALSE, hasClinical = FALSE){
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
  
  dim(ERP.pc1ihc)
  dim(ERN.pc1ihc)
  
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
  
  dim(mbal.pc1ihc) 
  
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
  fl.mdn = BreastSubtypeR$medians 
  
  df.al = merge(fl.mdn, df.mdns, by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  
  # normalization
  mat = docalibration( mat, df.al, calibration, internal)

  out = sspPredict(BreastSubtypeR$centroid, mat, std=F, distm="spearman", centroids=T, Prosigna = Prosigna)

  if(Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna ,IHC = df.cln$IHC , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, IHC = df.cln$IHC , row.names = NULL )
  }
  
  out$distances.prosigna =  -1 * out$distances.prosigna
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical,Prosigna = Prosigna )

  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns= df.al,outList=out))
  
  
}

#' Function for the second step of PCA-PAM50 approach
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param NM.mdns surffix name for medains
#' @param short surffix name for output
#' @noRd

makeCalls.v1PAM = function(mat, df.pam, calibration = "Internal", internal ="v1PAM.mdns", external=NA, medians = NA ,Prosigna =FALSE, hasClinical = FALSE,seed=118 ){
  
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
  
  fl.mdn = BreastSubtypeR$medians

  df.al = merge(fl.mdn, df.mdns, by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  

  ## normalization
  mat= docalibration( mat, df.al, calibration, internal)
  
  out = sspPredict(BreastSubtypeR$centroid, mat, std=F, distm="spearman", centroids=T, Prosigna = Prosigna)
  
  if(Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna, IHC = df.pam$IHC , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, IHC = df.pam$IHC , row.names = NULL )
  }
  out$distances.prosigna =  -1 * out$distances.prosigna
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln = df.pam, hasClinical = hasClinical,Prosigna = Prosigna )
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns.fl= df.al,outList=out))
  
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

makeCalls.ssBC = function(mat, df.cln, s , Prosigna = FALSE , hasClinical =FALSE  ){
  
  # ## test data
  # data("BreastSubtypeR")
  # gene.sigma = BreastSubtypeR$ssBC.subgroupQuantile
  # mat = exprs(UNC232)
  # df.cln= pData(UNC232)
  # df.cln$IHC = ""
  # 
  # s = "ER"
  
  # ## test data
  # mat = data_input$x_parker
  # df.cln = clinic.oslo
  # s = "ER_JAMA"
  # hasClinical =T

  gene.sigma = BreastSubtypeR$ssBC.subgroupQuantile
  

  if( s == "ER") { ## use ER selected strategy
    
    ## if there is no sample in either of both, wont influence the code
    ERN_samples = rownames(df.cln)[which(df.cln$ER == "ER-" )]
    ERP_samples = rownames(df.cln)[which(df.cln$ER == "ER+" )]
    samples_selected = list( ER_neg = ERN_samples, ER_pos = ERP_samples)
    
  } else if(s == "TN"){
    
    ## if there is no sample in either of both, wont influence the code
    TN_samples = rownames(df.cln)[which(df.cln$TN == "TN" )]
    samples_selected = list( TN = TN_samples)
    
  }  else if( s == "ER_JAMA" ){ #TNBC-JAMAOncol2024
    
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
    
    
    out = sspPredict(BreastSubtypeR$centroid, x.m , std=F, distm="spearman", centroids=T, Prosigna = Prosigna)

  
  }, names(samples_selected), SIMPLIFY = F, USE.NAMES = T )
  
  ## prepare data for ROR

  if( s == "ER") { ## use ER selected strategy
    
    predictions = c(res$ER_neg$predictions, res$ER_pos$predictions)
    distances = rbind(res$ER_neg$distances, res$ER_pos$distances )
    testData = cbind(res$ER_neg$testData, res$ER_pos$testData )
    distances.prosigna = rbind(res$ER_neg$distances.prosigna, res$ER_pos$distances.prosigna )
    
    if(Prosigna ){
      predictions.prosigna = c(res$ER_neg$predictions.prosigna, res$ER_pos$predictions.prosigna)
    }
    
  } else if(s == "ER_JAMA" ){
    
    predictions = c(res$ERneg_HER2neg$predictions, res$ERpos_HER2neg$predictions, 
                   res$HER2pos_ERneg$predictions, res$HER2pos_ERpos$predictions)
    distances = rbind(res$ERneg_HER2neg$distances, res$ERpos_HER2neg$distances,
                      res$HER2pos_ERneg$distances, res$HER2pos_ERpos$distances)
    testData = cbind(res$ERneg_HER2neg$testData, res$ERpos_HER2neg$testData,
                     res$HER2pos_ERneg$testData, res$HER2pos_ERpos$testData)
    distances.prosigna = rbind(res$ERneg_HER2neg$distances.prosigna, res$ERpos_HER2neg$distances.prosigna,
                               res$HER2pos_ERneg$distances.prosigna, res$HER2pos_ERpos$distances.prosigna)
    
    if(Prosigna ){
      predictions.prosigna = c(res$ERneg_HER2neg$predictions.prosigna, res$ERpos_HER2neg$predictions.prosigna, 
                      res$HER2pos_ERneg$predictions.prosigna, res$HER2pos_ERpos$predictions.prosigna)
    }
    
    
  } else if(s == "TN"){
    
    predictions = res$TN$predictions
    distances = res$TN$distances
    testData = res$TN$testData
    distances.prosigna = res$TN$distances.prosigna

    if(Prosigna ){
      predictions.prosigna = res$TN$predictions.prosigna
    }
    
  } else if (s == "TNBC")  {
    
    predictions = res$TNBC$predictions
    distances = res$TNBC$distances
    testData = res$TNBC$testData
    distances.prosigna = res$TNBC$distances.prosigna
    
    if(Prosigna ){
      predictions.prosigna = res$TNBC$predictions.prosigna
    }
  } 
  
  
  
  if(Prosigna){
    out = list(predictions=predictions,predictions.prosigna = predictions.prosigna ,testData=testData,distances=distances, distances.prosigna = distances.prosigna, centroids= BreastSubtypeR$centroid)
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna , IHC = df.cln[ which( rownames(df.cln) %in% names(out$predictions) ) ,"IHC"], row.names = NULL )
  } else {
    out= list(predictions=predictions, testData=testData,distances=distances,  distances.prosigna = distances.prosigna, centroids=BreastSubtypeR$centroid)
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, IHC = df.cln[ which( rownames(df.cln) %in% names(out$predictions) ) ,"IHC"], row.names = NULL )
    }
  out$distances.prosigna =  -1 * out$distances.prosigna
  out$distances = -1 * out$distances
  
  
  ##calculate and grouping
  out.ROR = RORgroup(out,df.cln, hasClinical = hasClinical, Prosigna = Prosigna )
  
  ## reorder
  orde.No = match(colnames(mat), Int.sbs$PatientID )
  Int.sbs = Int.sbs[orde.No,]
  out.ROR = out.ROR[colnames(mat),]
  out$distances.prosigna = out$distances.prosigna[orde.No,]
  out$distances = out$distances[orde.No,]
  out$testData =  out$testData[,colnames(mat) ]
  out$predictions =  out$predictions[colnames(mat)]
  
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns = gene.sigma, outList=out))
  
}
  

