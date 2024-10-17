#' 
#' Functions adapted from NC-based subtyping methods
#' @name NC-based
#' @import ggplot2
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import magrittr
#' @import impute
#' @importFrom dplyr select
#' @importFrom dplyr mutate_at
#' @noRd 
NULL

#' 
#' Function for central median
#' @param x gene expression matrix
#' @noRd 
medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

#' Function for quantile central
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


#' Function for calibration methods
#' @param y Gene expression matrix 
#' @param df.al Medians for calibration
#' @param calibration How to do calibration, "None"(default) means no calibration for gene expression matrix. When setting calibration =None, you dont need to set internal and external parameters.  "Internal" means calibration for gene expression matrix by itself. "External" means calibration by external cohort. 
#' @param internal Specify the strategy for internal calibration, medianCtr(default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians calculated by train cohorts. When users want to use Medians prepared by user selves, this parameter should be "Given.mdns", not platform name. 
#' @noRd 
docalibration = function( y, df.al, calibration = "None", internal=internal, external=external){
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
              medians = df.al
              #print(paste("calibration to:",internal))
              tm = overlapSets(medians,y)
              y = (tm$y-tm$x[,internal])
            }
            else { print( "Please choose internal calibration stategy: medianCtr, meanCtr, qCtr or ER+/- relevant ")}
          },
          "External" = { ## external
            medians = df.al 
            #print(paste("calibration to:",external)) ## pre-prepared medians and givenmedians
            tm = overlapSets(medians,y)
            y = (tm$y-tm$x[,external]) }
          
  )
  
  return(as.matrix( y) )
  
}


#' Function for standardize
#' @param x gene expression matrix
#' @noRd
standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}


#' Function for ordering gene in expression matrix as PAM50 genes 
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


#' Function for predict PAM50 subtyping
#' @param x median train file
#' @param y gene expression matrix
#' @param classes description
#' @param distm "euclidean" or "spearman" (default)
#' @param centrids Logic.
#' @param Prosigna Logic. Please specify if it predicts prosigna-like subtype
#' @noRd
sspPredict<-function(x, y, std=FALSE, distm="spearman", Prosigna = TRUE){
  
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
  
  temp = overlapSets(dataMatrix,tdataMatrix)
  dataMatrix = temp$x
  tdataMatrix = temp$y
  sfeatureNames = row.names(dataMatrix)
  
  # standardize both sets
  if(std){
    dataMatrix = standardize(dataMatrix)
    tdataMatrix = standardize(tdataMatrix)
  }

  
  nGenes = dim(dataMatrix)[1]
  #print(paste("Number of genes used:",nGenes))
  centroids = dataMatrix
  nClasses = dim(centroids)[2]
  classLevels = dimnames(centroids)[[2]]
  
  
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
        distances.prosigna.subtype[,j] = dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
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


#' Function for risk calculation
#' 
#' @param out The result of sspPredict() function. 
#' @param Prosigna Logic.
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @return ROR, ROR risk group and other indications
#' @noRd

RORgroup = function(out, df.cln , hasClinical = FALSE, Prosigna = FALSE ){
  
  # # ## test data
  # out
  # df.cln = df.cln
  # hasClinical = hasClinical

  ## prepare outtable

  sample = data.frame( patientID = names(out$predictions) )
  
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
    
    Clinical = df.cln[match(names(out$predictions),df.cln$PatientID ),]
    rownames(Clinical) = Clinical$PatientID

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
        
        Clinical$NODE = as.numeric(Clinical$NODE)
        
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
    
    
    outtable = cbind(sample, distance, Call, call.conf, ROR.genomic,ROR.combined, ROR.combined.prosigna, er_her2)
    
  } else {
    outtable = cbind(sample, distance, Call, call.conf, ROR.genomic, er_her2)
  }
  
  
  if(Prosigna){
    
    ## pass distances (just omitting normal)
    Call.prosigna = data.frame( "Call.prosigna" =out$predictions.prosigna, row.names = names(out$predictions.prosigna))
    outtable = cbind(outtable, Call.prosigna)
    outtable =outtable %>% dplyr::select(patientID, colnames(distance), Call, Call.prosigna, everything() )
    
  }
  
  return(outtable)
}



#' 
#' Function for predicting PAM50 intrinsic subtypes and calculation of proliferation 
#' 


####function to make calls using near-centroid strategy #### 
#' Function for calling PAM50 subtypes by parker based methods
#' Here, we integrated parker-based methods and genefu PAM50 model
#' @param mat gene expression matrix, log of normalized
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param calibration How to do calibration, "None"(default) means no calibration for gene expression matrix. When setting calibration =None, you dont need to set internal and external parameters.  "Internal" means calibration for gene expression matrix by itself. "External" means calibration by external cohort. 
#' @param internal Specify the strategy for internal calibration, medianCtr(default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians calculated by train cohorts. When users want to use Medians prepared by user selves, this parameter should be "Given.mdns", not platform name. 
#' @param medians If you specify "external" parameter as "Given.mdns", you should input matrix/table, 50 signatures in the first column and "Given.mdns" values in the second column.
#' @param Prosigna Logic. 
#' @param hasClinical Logic. Please specify if you prepared clinical information, like Tumore size as T column, lymphatic node status as NODE column. 
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
  
  out = sspPredict(BreastSubtypeR$centroid, mat, std=F, distm="spearman", Prosigna = Prosigna)
  
  
  if (Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS= out$predictions, row.names = NULL )
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
#' @param df.cln clicnical information table with PatientID 
#' @param calibration The calibration method to use. Options are "None", "Internal", or "External". If "Internal" is selected, see the "internal" parameter for further details. If "External" is selected, see the "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are median-centered ("medianCtr", default), mean-centered ("meanCtr"), or quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for external medians, which are calculated by the training cohort. If you want to use user-provided medians, set this parameter to "Given.mdns" and provide the medians via the "medians" parameter. 
#' @param medians If "Given.mdns" is specified for the "external" parameter, input a matrix/table where the first column contains 50 genes and the second column contains the corresponding "Given.mdns" values.
#' @param Prosigna Logic.
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @noRd
makeCalls.ihc = function(mat, df.cln, seed=118, calibration = "Internal", internal = "IHC.mdns", external=NA, medians = NA , Prosigna = FALSE , hasClinical = FALSE){
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
  
  out = sspPredict( centroids, mat, std=F, distm="spearman",  Prosigna = Prosigna)
  
  if (Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna, row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
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
#' @param calibration The calibration method to use. Options are "None", "Internal", or "External". If "Internal" is selected, see the "internal" parameter for further details. If "External" is selected, see the "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are median-centered ("medianCtr", default), mean-centered ("meanCtr"), or quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for external medians, which are calculated by the training cohort. If you want to use user-provided medians, set this parameter to "Given.mdns" and provide the medians via the "medians" parameter. 
#' @param medians If "Given.mdns" is specified for the "external" parameter, input a matrix/table where the first column contains 50 genes and the second column contains the corresponding "Given.mdns" values.
#' @param Prosigna Logic. 
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
  if(ratio > dim(ERP.ihc)[1]/dim(ERN.ihc)[1] ) {
    n_ERP = dim(ERP.ihc)[1]
    n_ERN = dim(ERN.ihc)[1]
    stop( paste0( "please specify a ratio less than ", max(n_ERP, n_ERN )/min(n_ERP, n_ERN )   )  ) 
    }
  
  
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
    
    out = sspPredict(centroids, mat, std=F, distm="spearman",Prosigna = Prosigna)
    
    return( out )
    
  },paste0("itr.", seq(iteration) ) ,SIMPLIFY = FALSE,USE.NAMES = TRUE)
  
  
  ## get consensus intrinsic subtype
  Call_subtypes = mapply(function(res_ihs){res_ihs$predictions }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
  consensus_subtypes = apply(Call_subtypes, 1, get_consensus_subtype)
  
  if (Prosigna) {
    Call_subtypes.prosigna = mapply(function(res_ihs){res_ihs$predictions.prosigna }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
    consensus_subtypes.prosigna = apply(Call_subtypes.prosigna, 1, get_consensus_subtype)
    
    Int.sbs = data.frame(PatientID = names(consensus_subtypes), BS = consensus_subtypes, BS.prosigna = consensus_subtypes.prosigna , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID =names(consensus_subtypes), BS = consensus_subtypes, row.names = NULL )
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
#' @param df.cln clinical information table
#' @param calibration The calibration method to use. Options are "None", "Internal", or "External". If "Internal" is selected, see the "internal" parameter for further details. If "External" is selected, see the "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are median-centered ("medianCtr", default), mean-centered ("meanCtr"), or quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for external medians, which are calculated by the training cohort. If you want to use user-provided medians, set this parameter to "Given.mdns" and provide the medians via the "medians" parameter. 
#' @param medians If "Given.mdns" is specified for the "external" parameter, input a matrix/table where the first column contains 50 genes and the second column contains the corresponding "Given.mdns" values.
#' @param Prosigna Logic. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
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
  
  # ## test data
  # mat = data_input$x_NC.log
  # df.cln = clinic.oslo
  # calibration = "Internal"
  # internal="PC1ihc.mdns"
  # seed=118
  # hasClinical = T
  
  
  # Pull the PCA components
  #rv      = rowVars(mat)
  #select  = order(rv, decreasing = TRUE)[seq_len(dim(mat)[1])] # the input is PAM50 matrix --50 genes -- get from dimension
  pca     = prcomp(t(mat))#[select,]
  pc12    = pca$x[,1:2] #get two principal 
  df.pc1  = data.frame(PatientID=rownames(pc12),PC1 = pc12[,1],stringsAsFactors=F)
  df.pca1 = left_join(df.cln,df.pc1,by="PatientID")
  
  ### function to count the number of mis-classified cases on a given PC1 point ---find the cutoff
  getno = function(x){
    p.rgt = length(which(df.pca1$ER %in% c("ER+") & df.pca1$PC1 >x))/length(which(df.pca1$ER %in% c("ER+") ))
    n.lft = length(which(df.pca1$ER %in% c("ER-") & df.pca1$PC1 <x))/ length(which(df.pca1$ER %in% c("ER-")))
    tot   = (p.rgt + n.lft) * 100
    return(list(PC1=x,Mis=tot))
  }
  
  
  df.mis  = do.call(rbind.data.frame,lapply(seq(-20,20,by=0.1),getno))
  plot(df.mis)
  
  num.min = df.mis$PC1[which(df.mis$Mis == min(df.mis$Mis))]
  
  ERP.pc1ihc = df.pca1[which(df.pca1$ER %in% c("ER+") & df.pca1$PC1 <= mean(num.min)),] # used mean to overcome situation where there are two minimum
  ERN.pc1ihc = df.pca1[which(df.pca1$ER %in% c("ER-") & df.pca1$PC1 > mean(num.min)),]
  
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
  
  out = sspPredict(BreastSubtypeR$centroid, mat, std=F, distm="spearman", Prosigna = Prosigna)
  
  if(Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
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
#' @param mat gene expression matrix 
#' @param df.pam clinical information table created using makeCalls.PC1ihc().  
#' @param calibration The calibration method to use, "Internal". 
#' @param internal Specify the strategy for internal calibration, "v1PAM.mdns".
#' @param external NA
#' @param medians NA
#' @param Prosigna Logic. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
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
  
  out = sspPredict(BreastSubtypeR$centroid, mat, std=F, distm="spearman", Prosigna = Prosigna)
  
  if(Prosigna) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna, row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
  }
  out$distances.prosigna =  -1 * out$distances.prosigna
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln = df.pam, hasClinical = hasClinical,Prosigna = Prosigna )
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns.fl= df.al,outList=out))
  
}



#' Function for calling PAM50 subtypes by ssBC methods
#' This function is adapted from ssBC TNBC-BreastCancerRes2015 and subgrou specific TNBC-JAMAOncol2024 
#' @param mat gene expression matrix
#' @param df.cln clinical information table. The first column must be named "PatientID".
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @param s Options are "ER" or "TN" or "ER_JAMA" or "HER2+" or "TNBC". Specify the medians you want. The original quantile is "ER" and "TN" of TNBC-BreastCancerRes2015.  If you choose "ER_JAMA" or "HER2+" or "TNBC", it means you choose quantile from TNBC-JAMAOncol2024. 
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
    
    
    out = sspPredict(BreastSubtypeR$centroid, x.m , std=F, distm="spearman", Prosigna = Prosigna)
    
    
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
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.prosigna = out$predictions.prosigna , row.names = NULL )
  } else {
    out= list(predictions=predictions, testData=testData,distances=distances,  distances.prosigna = distances.prosigna, centroids=BreastSubtypeR$centroid)
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
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


