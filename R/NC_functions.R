#' 
#' Functions adapted from NC-based subtyping methods
#' @import magrittr
#' @import impute
#' @importFrom dplyr select
#' @noRd 
NULL

#' 
#' Function for central median
#' @param x gene expression matrix
#' @noRd 
medianCtr=function(x){
  annAll = dimnames(x)
  medians = apply(x,1,median,na.rm=TRUE)
  x = t(scale(t(x),center=medians,scale=FALSE))
  dimnames(x) = annAll
  return(x)
}

#' Function to rescale values based on quantiles
#' This is adapted from genefu package
#' @param x Gene expression matrix or vector
#' @noRd 
rescale = function(x, na.rm=FALSE, q=0) {

  if(q == 0) {
    ma = max(x, na.rm=na.rm)
    mi = min(x, na.rm=na.rm)
  } else {
    ma = quantile(x, probs=1-(q/2), na.rm=na.rm)
    mi = quantile(x, probs=q/2, na.rm=na.rm)
  }
  xx = (x - mi) / (ma - mi)
  
  return(xx)
}

#' Function for ordering genes in expression matrix as PAM50 genes 
#' @param x centroid matrix
#' @param y Gene expression matrix
#' @noRd
overlapSets=function(x,y){
  
  # subset the two lists to have a commonly ordered gene list
  x=x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
  y=y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]
  
  #and sort such that thing are in the correct order
  x=x[sort.list(row.names(x)),]
  y=y[sort.list(row.names(y)),]
  
  return(list(x=x,y=y))
}

#' Function for calibration methods
#' @param y Gene expression matrix 
#' @param df.al Medians for calibration
#' @param calibration How to do calibration, "None"(default) means no calibration for gene expression matrix. When setting calibration =None, you dont need to set internal and external parameters.  "Internal" means calibration for gene expression matrix by itself. "External" means calibration by external cohort. 
#' @param internal Specify the strategy for internal calibration, medianCtr(default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians calculated by train cohorts. When users want to use Medians prepared by user selves, this parameter should be "Given.mdns", not platform name. 
#' @noRd 
docalibration = function( y, df.al, calibration = "None", internal=internal, external=external){

  mq = 0.05 ## presetting in genefu robust model
  switch( calibration,
          "None" = {print("No calibration") }, 
          "Internal" = { ## internal
            if(internal == "medianCtr"){ y = medianCtr(y) } ## parker default method
            else if(internal == "meanCtr") { y =  t(scale( t(y), center=TRUE, scale=TRUE))} ## "scale" in genefu method
            else if(internal == "qCtr") { y = t(apply(y, 1, function(x) {  return( (rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2 ) }) )  } ## "robust" in genefu method; mp = 0.05
            else if(internal == internal) { ## which column to use
              medians = df.al
              
              tm = overlapSets(medians,y)
              y = (tm$y-tm$x[,internal])
            }
            else { message( "Please choose internal calibration stategy: medianCtr, meanCtr, qCtr.")}
          },
          "External" = { ## external
            medians = df.al 
            #print(paste("calibration to:",external)) ## pre-prepared medians and givenmedians
            tm = overlapSets(medians,y)
            y = (tm$y-tm$x[,external]) }
          
  )
  
  return(as.matrix( y) )
  
}


#' Function for standardization
#' @param x gene expression matrix
#' @noRd
standardize=function(x){
  annAll=dimnames(x)
  x=scale(x)
  dimnames(x)=annAll
  return(x)
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


#' Function for predicting PAM50 subtypes
#' @param x median train file
#' @param y gene expression matrix
#' @param std Logic
#' @param distm "euclidean" or "spearman" (default)
#' @param centrids Logic
#' @param Subtype Logic. Please specify if it predicts Subtype-like subtype
#' @noRd
sspPredict=function(x, y, std=FALSE, distm="spearman", Subtype = TRUE){
  
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
      distances[,j]= dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
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
  
  ## run four Subtypes
  if( Subtype){
    
    nClasses = nClasses-1 ## omitting normal
    classLevels = classLevels[1:4]
    
    distances.Subtype.subtype = matrix(ncol= nClasses, nrow=dim(tdataMatrix)[2])
    for(j in 1:nClasses){
      if(distm=="euclidean"){
        distances.Subtype.subtype[,j] = dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
      }
      if(distm=="correlation" | distm=="pearson"){
        distances.Subtype.subtype[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids[,j], x, method = "pearson", use = "pairwise.complete.obs"))
      }
      if(distm=="spearman"){
        distances.Subtype.subtype[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids[,j], x, method = "spearman", use = "pairwise.complete.obs"))
      }
    }
    
    prediction.Subtype = classLevels[apply(distances.Subtype.subtype, 1, which.min,simplify = TRUE)]
    names(prediction.Subtype) =  colnames(tdataMatrix)
    
  }
  
  
  ## prepare Subtype.distances
  ## genes to be excluded
  genes.ex = c("BIRC5", "CCNB1", "GRB7","MYBL2")
  
  #parse the test file - same as train file but no rows of classes
  tdataMatrix = y[ which( !(rownames(y) %in%  genes.ex ) ),]
  
  #dimnames(tdataMatrix)[[2]]=paste("x",seq(1,471))
  temp = overlapSets(dataMatrix,tdataMatrix)
  dataMatrix = temp$x
  tdataMatrix = temp$y

  ## omitting normal
  dataMatrix = dataMatrix[,1:4] ## omitting normal
  nGenes = dim(dataMatrix)[1]
  centroids.Subtype = dataMatrix
  nClasses =dim(centroids.Subtype)[2] ## four subtypes
  classLevels = dimnames(centroids.Subtype)[[2]]
  
  distances.Subtype  =  matrix(ncol= nClasses,nrow=dim(tdataMatrix)[2])
  for(j in 1:nClasses){
    if(distm=="euclidean"){
      distances.Subtype[,j] = dist(t(cbind(centroids.Subtype[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
    }
    if(distm=="correlation" | distm=="pearson"){
      distances.Subtype[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids.Subtype[,j], x, method = "pearson", use = "pairwise.complete.obs"))
    }
    if(distm=="spearman"){
      distances.Subtype[,j] = apply(tdataMatrix, 2, function(x) -cor(centroids.Subtype[,j], x, method = "spearman", use = "pairwise.complete.obs"))
    }
  }

  ## return 
  if ( Subtype) {
    res = list(predictions=prediction, predictions.Subtype = prediction.Subtype,testData= as.matrix(y),distances=distances,distances.Subtype = distances.Subtype, distances.Subtype.subtype = distances.Subtype.subtype, centroids=centroids)
    
  }else {
    res = list(predictions=prediction,testData= as.matrix(y), distances=distances, distances.Subtype = distances.Subtype,centroids=centroids)
  }
  
  return(res)
}


#' Function for risk calculation
#' 
#' @param out The result of sspPredict() function. 
#' @param out clinical table
#' @param Subtype Logic.
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @return ROR, ROR risk group and other indications
#' @noRd

RORgroup = function(out, df.cln, Subtype = FALSE, hasClinical = FALSE ){

  sample = data.frame( patientID = names(out$predictions) )
  
  distance = data.frame(out$distances, row.names = names(out$predictions) )
  colnames(distance) = c("Basal","Her2","LumA","LumB","Normal")
  
  Call = data.frame( "Call" = out$predictions, row.names = names(out$predictions))
  
  #providing proliferation signatures
  proliferationGenes = c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
  proliferationGenes.Subtype = c("ANLN", "CEP55", "ORC6L", "CCNE1", "EXO1", "PTTG1", "CDC20", "KIF2C", "RRM2", "CDC6", "KNTC2", "TYMS", "CDCA1", 
                                  "MELK", "UBE2C", "CENPF", "MKI67", "UBE2T")
  
  ###
  # some constants for ROR groups
  ###
  
  # for subtype only model
  glthreshold= -0.15
  ghthreshold=  0.1
  
  # for subtype + proliferation model
  gplthreshold= -0.25
  gphthreshold=  0.1
  
  # for combined model
  clthreshold= -0.1
  chthreshold=  0.2
  
  # for combined + proliferation model
  cplthreshold= -0.2
  cphthreshold=  0.2
  
  ## ??? 
  # for combined + proliferation model + Subtype with negative node status 
  cplthreshold.Subtype.NODE.0 = 40
  cphthreshold.Subtype.NODE.0 = 60
  
  # for combined + proliferation model + Subtype with node status
  ## first 4 or more nodes would be assigned as high-risk
  cplthreshold.Subtype.NODE = 15
  cphthreshold.Subtype.NODE = 40
  
  
  ## ER and HER2 score
  #out$testData = mat
  erScore = out$testData["ESR1",]
  her2Score = out$testData["ERBB2",]
  
  er_her2 = data.frame( "ER" = erScore, "HER2" = her2Score, row.names = colnames(out$testData ))
  
  # calculate the proliferation score
  prolifScore = apply(out$testData[ which( rownames(out$testData) %in% proliferationGenes), ], 2, mean, na.rm =TRUE )
  
  prolifScore.Subtype = apply(out$testData[ which( rownames(out$testData) %in% proliferationGenes.Subtype) , ], 2, mean, na.rm =TRUE )
  
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

  genomicWprolif = -0.0009299747*out$distances[,1] + 0.0692289192*out$distances[,2] + -0.0951505484*out$distances[,3] +  0.0493487685*out$distances[,4] + 0.3385116381*prolifScore
  
  
  # threshold the risk score
  griskgroups=genomic
  griskgroups[genomic>ghthreshold]="high"
  griskgroups[genomic>glthreshold & genomic<ghthreshold]="med"
  griskgroups[genomic<glthreshold]="low"
  gpriskgroups=genomicWprolif
  gpriskgroups[genomicWprolif>gphthreshold]="high"
  gpriskgroups[genomicWprolif>gplthreshold & genomicWprolif<gphthreshold]="med"
  gpriskgroups[genomicWprolif<gplthreshold]="low"
  
  genomic= 100* (genomic + 0.35 ) / 0.85
  genomicWprolif= 100* (genomicWprolif + 0.35 ) / 0.85
  
  
  ROR.genomic = data.frame("ROR-S (Subtype Only)" = genomic, "ROR-S Group (Subtype Only)" = griskgroups,
                           "ROR-P (Subtype + Proliferation)" = genomicWprolif, "ROR-P Group (Subtype + Proliferation)" = gpriskgroups,
                           check.names = FALSE)
  
  if (hasClinical){
    
    Clinical = df.cln[match(names(out$predictions),df.cln$PatientID ),]
    rownames(Clinical) = Clinical$PatientID

    if ( "TSIZE" %in% colnames(Clinical)) {
      
      xT= as.numeric(as.vector(Clinical$TSIZE))
      
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

      combinedWprolif.Subtype = 54.7690 * (-0.0067 * out$distances.Subtype[,1] + 0.4317*out$distances.Subtype[,2] - 0.3172*out$distances.Subtype[,3] + 
                                              0.4894*out$distances.Subtype[,4] + 0.1981*prolifScore.Subtype + 0.1133*xT + 0.8826)
      
      
      ## check NODE 
      ## grouping by NODE
      if ("NODE" %in% colnames(Clinical)  ) {
        
        Clinical$NODE = as.numeric(Clinical$NODE)
        
        cpriskgroups.Subtype = combinedWprolif.Subtype
        
        patients.NODE = rownames(Clinical)[ !is.na(Clinical$NODE )]
        
        if( length( patients.NODE ) > 0  )  {
          
          ## NODE >= 4
          patients.NODE.4 = rownames(Clinical)[ Clinical$NODE >= 4 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.4 ) > 0  ) { 
            
            cpriskgroups.Subtype[patients.NODE.4] = "high"
          }
          
          ## NODE < 4 & >0
          patients.NODE.3 = rownames(Clinical)[ Clinical$NODE < 4 & Clinical$NODE > 0 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.3 ) > 0  ) { 
            
            temp = combinedWprolif.Subtype[patients.NODE.3]
            cpriskgroups.Subtype[patients.NODE.3][ temp> cphthreshold.Subtype.NODE ]  = "high"
            cpriskgroups.Subtype[patients.NODE.3][temp>cplthreshold.Subtype.NODE & temp<cphthreshold.Subtype.NODE ] = "med"
            cpriskgroups.Subtype[patients.NODE.3][temp<cplthreshold.Subtype.NODE] = "low"
            
          }
          
          ## negative node
          patients.NODE.0 = rownames(Clinical)[ Clinical$NODE == 0 & !is.na(Clinical$NODE )]
          
          if (length( patients.NODE.0 ) > 0  ) { 
            
            temp = combinedWprolif.Subtype[patients.NODE.0]
            cpriskgroups.Subtype[patients.NODE.0][ temp> cphthreshold.Subtype.NODE.0 ]  = "high"
            cpriskgroups.Subtype[patients.NODE.0][temp>cplthreshold.Subtype.NODE.0 & temp<cphthreshold.Subtype.NODE.0 ] = "med"
            cpriskgroups.Subtype[patients.NODE.0][temp<cplthreshold.Subtype.NODE.0] = "low"
            
          }
          
          
          ## check if NA NODE 
          ## grouping by NA NODE
          patients.NA = rownames(Clinical)[ is.na(Clinical$NODE ) | is.na(Clinical$T) ]
          
          if ( length( patients.NA ) > 0  ){
            
            cpriskgroups.Subtype[patients.NA] = NA
            
          }
          
          ROR.combined.Subtype = data.frame("ROR-PC (Subtype + Clinical + Proliferation.Subtype)" = combinedWprolif.Subtype, 
                                             "ROR-PC Group (Subtype + Clinical + Proliferation.Subtype)" = cpriskgroups.Subtype,
                                             check.names = FALSE)
          
          
        } else {
          
          cpriskgroups.Subtype[combinedWprolif.Subtype > cphthreshold.Subtype] = "high"
          cpriskgroups.Subtype[combinedWprolif.Subtype > cplthreshold.Subtype & combinedWprolif.Subtype < cphthreshold.Subtype] = "med"
          cpriskgroups.Subtype[combinedWprolif.Subtype < cplthreshold.Subtype] = "low"
          
          
          ROR.combined.Subtype = data.frame("ROR-PC (Subtype + Clinical + Proliferation.Subtype)" = combinedWprolif.Subtype, 
                                             "ROR-PC Group (Subtype + Clinical + Proliferation.Subtype)" = cpriskgroups.Subtype,
                                             check.names = FALSE)
          
        }
        
      } 
      
    }
    
    
    outtable = cbind(sample, distance, Call, call.conf, ROR.genomic,ROR.combined, ROR.combined.Subtype, er_her2)
    
  } else {
    outtable = cbind(sample, distance, Call, call.conf, ROR.genomic, er_her2)
  }
  
  
  if(Subtype){
    
    ## pass distances (just omitting normal)
    Call.Subtype = data.frame( "Call.Subtype" =out$predictions.Subtype, row.names = names(out$predictions.Subtype))
    outtable = cbind(outtable, Call.Subtype)
    outtable =outtable %>% dplyr::select(patientID, colnames(distance), Call, Call.Subtype, everything() )
    
  }
  
  return(outtable)
}



#' 
#' Functions for predicting PAM50 intrinsic subtypes and calculation of proliferation 
#' 


#### functions to make calls using nearest-centroid (NC-based) strategies #### 

#' Function for calling PAM50 subtypes by Parker et al.-based methods
#' Here, we integrated Parker et al.-based methods and genefu PAM50 implementation
#' @param mat gene expression matrix, log of normalized
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param calibration How to do calibration, "None"(default) means no calibration for gene expression matrix. When setting calibration =None, you dont need to set internal and external parameters.  "Internal" means calibration for gene expression matrix by itself. "External" means calibration by external cohort. 
#' @param internal Specify the strategy for internal calibration, medianCtr(default), meanCtr and qCtr
#' @param external Specify the platform name(which column) of external medians calculated by train cohorts. When users want to use Medians prepared by user selves, this parameter should be "Given.mdns", not platform name. 
#' @param medians If you specify "external" parameter as "Given.mdns", you should input matrix/table, 50 signatures in the first column and "Given.mdns" values in the second column.
#' @param Subtype Logic. 
#' @param hasClinical Logic. Please specify if you prepared clinical information, like Tumore size as T column, lymphatic node status as NODE column. 
#' @noRd
#' 

makeCalls.parker = function(mat, df.cln, calibration = "None", internal = NA,external=NA, medians = NA,Subtype = FALSE, hasClinical =FALSE  ){

  ## loading dataset
  data("BreastSubtypeRobj")
  
  fl.mdn = BreastSubtypeRobj$medians
  
  
  if(calibration == "External" & external == "Given.mdns" ) {
    
    if( length(medians) == 1 || is.na( medians) ){
      stop("Please input prepared medians, genes in the first column and median values in the second column. ")
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
  
  out = sspPredict(BreastSubtypeRobj$centroid, mat, std=FALSE, distm="spearman", Subtype = Subtype)
  
  
  if (Subtype) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.Subtype = out$predictions.Subtype , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS= out$predictions, row.names = NULL )
  }
  out$distances.Subtype =  -1 * out$distances.Subtype
  out$distances = -1 * out$distances
  
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical )
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns = df.al, outList=out))
  
  
}


#### function to form an ER-balanced subset and derive its median

#' Function to form an ER-balanced subset and derive its median
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID 
#' @param calibration The calibration method to use. Options are "None", "Internal", or "External". If "Internal" is selected, see the "internal" parameter for further details. If "External" is selected, see the "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are median-centered ("medianCtr", default), mean-centered ("meanCtr"), or quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for external medians, which are calculated by the training cohort. If you want to use user-provided medians, set this parameter to "Given.mdns" and provide the medians via the "medians" parameter. 
#' @param medians If "Given.mdns" is specified for the "external" parameter, input a matrix/table where the first column contains 50 genes and the second column contains the corresponding "Given.mdns" values.
#' @param Subtype Logic.
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @param seed An integer value is used to set the random seed.
#' @noRd

makeCalls.ihc = function(mat, df.cln, calibration = "Internal", internal = "IHC.mdns", external=NA, medians = NA , Subtype = FALSE , hasClinical = FALSE, seed=118){
  ## loading dataset
  data("BreastSubtypeRobj")
  
  ERN.ihc = df.cln[which(df.cln$ER == "ER-"),] ### get ER- samples data.frame
  dim(ERN.ihc)	#[1] 153   9
  
  ERP.ihc = df.cln[which(df.cln$ER == "ER+"),]
  dim(ERP.ihc) #[1] 559   9
  
  #seed = 118
  if( dim(ERN.ihc)[1] > dim(ERP.ihc)[1] ) {  temp = ERP.ihc; ERP.ihc = ERN.ihc; ERN.ihc = temp }
  # set.seed(seed); i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1]) # take equal number of ER+ and ER- samples
  withr::with_seed(seed, {
    i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1]) # take equal number of ER+ and ER- samples
  })
  
  length(ERP.ihc$PatientID[i]) # ER positive samples
  length(ERN.ihc$PatientID)    # ER negative samples
  

  mbal.ihc = mat[,c(ERP.ihc$PatientID[i],ERN.ihc$PatientID)]
  
  dim(mbal.ihc) 
  
  # Find median
  surffix = getsurffix(calibration = calibration, internal,external)
  
  mdns      = apply(mbal.ihc,1,median,na.rm=TRUE) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
  
  colnames(df.mdns) = c("X",surffix)
  
  ## merge mdns
  fl.mdn = BreastSubtypeRobj$medians
  
  df.al = merge(fl.mdn, df.mdns , by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  ## centroids
  centroids = BreastSubtypeRobj$centroid #pam50_centroids.txt
  
  ## normalization
  mat = docalibration( mat, df.al, calibration, internal)
  
  out = sspPredict( centroids, mat, std=FALSE, distm="spearman",  Subtype = Subtype)
  
  if (Subtype) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.Subtype = out$predictions.Subtype, row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
  }
  
  out$distances.Subtype =  -1 * out$distances.Subtype
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical, Subtype = Subtype)
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns= df.al, outList=out))
  
}



#' Function for iterative ER-balanced subset gene centering 
#' @param mat gene expression matrix 
#' @param df.cln clicnical information table with PatientID and IHC column
#' @param iterative Times to do iterative ER balanced procedure with certain ratio. 
#' @param ratio The options are either 1:1 or 54 (ER+) : 64 (ER-) (default). The latter was ER ratio used for UNC230 train cohort.
#' @param calibration The calibration method to use. Options are "None", "Internal", or "External". If "Internal" is selected, see the "internal" parameter for further details. If "External" is selected, see the "external" parameter.
#' @param internal Specify the strategy for internal calibration. Options are median-centered ("medianCtr", default), mean-centered ("meanCtr"), or quantile-centered ("qCtr").
#' @param external Specify the platform name (i.e., the column name) for external medians, which are calculated by the training cohort. If you want to use user-provided medians, set this parameter to "Given.mdns" and provide the medians via the "medians" parameter. 
#' @param medians If "Given.mdns" is specified for the "external" parameter, input a matrix/table where the first column contains 50 genes and the second column contains the corresponding "Given.mdns" values.
#' @param Subtype Logic. If `TRUE`, the function predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from the `phenodata` table. Required columns include:
#' @param seed An integer value is used to set the random seed.
#' @noRd

makeCalls.ihc.iterative = function( mat, df.cln, iteration = 100, ratio = 54/64, calibration = "Internal", internal = "ER.mdns", external=NA, medians = NA , Subtype = FALSE, hasClinical = FALSE, seed=118){
  
  ## loading dataset
  data("BreastSubtypeRobj")
  
  # load the published centroids for classifcation
  centroids = BreastSubtypeRobj$centroid #pam50_centroids.txt
  
  ## preprocess the input matrix
  ### get ER- samples
  ERN.ihc = df.cln[which(df.cln$ER =="ER-"),] 
  dim(ERN.ihc)
  
  ### get ER+ samples
  ERP.ihc = df.cln[which(df.cln$ER =="ER+"),]
  dim(ERP.ihc) 
  
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

    
    mbal.ihc = mat[,c(ERP.ihc$PatientID[i],ERN.ihc$PatientID)]

    surffix = getsurffix(calibration = calibration, internal)
    
    # Calculate median
    mdns      = apply(mbal.ihc,1,median,na.rm=TRUE) # compute median of each row i.e gene
    mdns.df   = as.data.frame(mdns)
    df.mdns   = data.frame(X=rownames(mdns.df),mdns.ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
    colnames(df.mdns) = c("X",surffix)
    
    ## integrate ihc.mdns
    fl.mdn = BreastSubtypeRobj$medians
    
    df.al = merge(fl.mdn, df.mdns , by = "X")
    rownames(df.al) = df.al$X
    df.al = df.al[,-1]
    
    
    ## normalization
    mat = docalibration( mat, df.al, calibration,internal)
    
    out = sspPredict(centroids, mat, std=FALSE, distm="spearman",Subtype = Subtype)
    
    return( out )
    
  },paste0("itr.", seq(iteration) ) ,SIMPLIFY = FALSE,USE.NAMES = TRUE)
  
  
  ## get consensus intrinsic subtype
  Call_subtypes = mapply(function(res_ihs){res_ihs$predictions }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
  consensus_subtypes = apply(Call_subtypes, 1, get_consensus_subtype)
  
  if (Subtype) {
    Call_subtypes.Subtype = mapply(function(res_ihs){res_ihs$predictions.Subtype }, res_ihc_iterative, SIMPLIFY = TRUE,USE.NAMES = TRUE )
    consensus_subtypes.Subtype = apply(Call_subtypes.Subtype, 1, get_consensus_subtype)
    
    Int.sbs = data.frame(PatientID = names(consensus_subtypes), BS = consensus_subtypes, BS.Subtype = consensus_subtypes.Subtype , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID =names(consensus_subtypes), BS = consensus_subtypes, row.names = NULL )
  }

  mean_eve = get_average_subtype(res_ihc_iterative, consensus_subtypes)

  if (Subtype) {
    out = list(predictions = consensus_subtypes, predictions.Subtype = consensus_subtypes.Subtype, testData = mean_eve$testdata, distances = -1 * mean_eve$mean_distance, distances.Subtype = -1 * mean_eve$mean_distance.Subtype, centroids = centroids )
  } else {
    out = list(predictions = consensus_subtypes, testData = mean_eve$testdata,  distances = -1 * mean_eve$mean_distance,  distances.Subtype = -1 * mean_eve$mean_distance.Subtype, centroids = centroids )
  }
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical,Subtype = Subtype )

  if(Subtype ) {
    res =  list(BS.all = Int.sbs, score.ROR= out.ROR, outList = out, BS.itr.keep = Call_subtypes, BS.itr.keep = Call_subtypes.Subtype )
  } else {
    res = list(BS.all = Int.sbs, score.ROR= out.ROR, outList = out, BS.itr.keep = Call_subtypes )
  }
  
  return(res)
}


#### form secondary ER-balanced set (refer to paper) leveraging PCA and subsequent intermediate intrinsic subtypes

#' Function for the first step of PCA-PAM50 approach
#' @param mat gene expression matrix 
#' @param df.cln clinical information table
#' @param calibration The calibration method to use, "Internal". 
#' @param internal Specify the strategy for internal calibration, "PC1ihc.mdns".
#' @param external NA
#' @param medians NA
#' @param Subtype Logic
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @param seed 
#' @noRd

makeCalls.PC1ihc = function(mat, df.cln, calibration = "Internal", internal ="PC1ihc.mdns", external=NA, medians = NA ,Subtype =FALSE, hasClinical = FALSE, seed=118){
  ## loading dataset
  data("BreastSubtypeRobj")
  
  # Initial checks for 'df.cln' and 'mat'
  if (is.null(df.cln) || !'PatientID' %in% colnames(df.cln) || !'IHC' %in% colnames(df.cln)) {
    stop("Clinical data 'df.cln' is missing or does not contain required columns 'PatientID' and 'IHC'. Refer to the vignette to create 'test.clinical' data.")
  }
  
  if (is.null(mat) || !is.matrix(mat)) {
    stop("Gene expression matrix 'mat' is missing or not correctly formatted. Refer to the vignette to create 'test.matrix' data.")
  }
  
  # Pull the PCA components
  #rv      = rowVars(mat)
  #select  = order(rv, decreasing = TRUE)[seq_len(dim(mat)[1])] # the input is PAM50 matrix --50 genes -- get from dimension
  pca     = prcomp(t(mat))#[select,]
  pc12    = pca$x[,1:2] #get two principal
  df.pc1  = data.frame(PatientID=rownames(pc12),PC1 = pc12[,1],stringsAsFactors=FALSE)
  df.pca1 = merge(df.cln,df.pc1,by="PatientID")
  
  
  #--our function works best if majority of ER- cases fall in the positive PC1 axis--check
  # Identify ER-negative cases
  er_negative <- !grepl("^L", df.pca1$IHC)
  
  # Determine if the majority of ER-negative cases fall in the negative axis of PC1
  if (sum(df.pca1$PC1[er_negative] < 0) > sum(df.pca1$PC1[er_negative] >= 0)) {
    #print("yes")
    df.pca1$PC1 <- -df.pca1$PC1
  }
  
  # Ensure that IHC is not a factor or has all necessary levels defined
  if (is.factor(df.pca1$IHC)) {
    df.pca1$IHC = as.character(df.pca1$IHC)
  }
  
  # Convert IHC column to uppercase to handle case insensitivity
  df.pca1$IHC <- toupper(df.pca1$IHC)

  # Function to count the number of misclassified cases on a given PC1 point ---find the cutoff
  getno = function(x) {
    p.rgt = length(which(grepl("^L", df.pca1$IHC) & df.pca1$PC1 > x)) / length(which(grepl("^L", df.pca1$IHC)))
    n.lft = length(which(!grepl("^L", df.pca1$IHC) & df.pca1$PC1 < x)) / length(which(!grepl("^L", df.pca1$IHC)))
    tot = (p.rgt + n.lft) * 100
    return(list(PC1 = x, Mis = tot))
  }
  
  df.mis  = do.call(rbind.data.frame,lapply(seq(-20,20,by=0.1),getno))

  num.min = df.mis$PC1[which(df.mis$Mis == min(df.mis$Mis))]
  
  ERP.pc1ihc = df.pca1[which(grepl("^L", df.pca1$IHC) & df.pca1$PC1 <= mean(num.min)), ] # used mean to overcome situation where there are two minimum
  ERN.pc1ihc = df.pca1[which(!grepl("^L", df.pca1$IHC) & df.pca1$PC1 > mean(num.min)), ]
  
  dim(ERP.pc1ihc)
  dim(ERN.pc1ihc)
  
  if(dim(ERP.pc1ihc)[1] < dim(ERN.pc1ihc)[1] ){
    temp = ERN.pc1ihc
    ERN.pc1ihc =  ERP.pc1ihc
    ERP.pc1ihc = temp
    rm(temp)
  }
  
  # set.seed(seed);i = sample(dim(ERP.pc1ihc)[1],dim(ERN.pc1ihc)[1]) # take equal number of ER+ and ER- samples
  withr::with_seed(seed, {
    i = sample(dim(ERP.pc1ihc)[1],dim(ERN.pc1ihc)[1]) # take equal number of ER+ and ER- samples
  })
  length(ERP.pc1ihc$PatientID[i]) # ER positive samples
  length(ERN.pc1ihc$PatientID)    # ER negative samples
  
  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pc1ihc = mat[,c(ERP.pc1ihc$PatientID[i],ERN.pc1ihc$PatientID)]
  
  dim(mbal.pc1ihc) 
  
  # Find median
  
  surffix = getsurffix(calibration = calibration, internal)
  
  mdns      = apply(mbal.pc1ihc,1,median,na.rm=TRUE) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.pc1ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
  colnames(df.mdns) = c("X",surffix)
  
  ## medians
  fl.mdn = BreastSubtypeRobj$medians 
  
  df.al = merge(fl.mdn, df.mdns, by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  
  # normalization
  mat = docalibration( mat, df.al, calibration, internal)
  
  out = sspPredict(BreastSubtypeRobj$centroid, mat, std=FALSE, distm="spearman", Subtype = Subtype)
  
  if(Subtype) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.Subtype = out$predictions.Subtype , row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
  }
  
  out$distances.Subtype =  -1 * out$distances.Subtype
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln, hasClinical = hasClinical,Subtype = Subtype )
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns= df.al,outList=out))
  

}

#' Function for the second step of PCA-PAM50 approach
#' @param mat gene expression matrix 
#' @param df.pam clinical information table created using makeCalls.PC1ihc().  
#' @param calibration The calibration method to use, "Internal". 
#' @param internal Specify the strategy for internal calibration, "v1PAM.mdns".
#' @param external NA
#' @param medians NA
#' @param Subtype Logic. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @param seed An integer value is used to set the random seed.
#' @noRd

makeCalls.v1PAM = function(mat, df.pam, calibration = "Internal", internal ="v1PAM.mdns", external=NA, medians = NA ,Subtype =FALSE, hasClinical = FALSE,seed=118 ){
  
  ## loading dataset
  data("BreastSubtypeRobj")
  
  ERN.pam = df.pam[which(df.pam$PAM50 %in% c("Basal")),] ### get ER- samples data.frame
  dim(ERN.pam)
  
  ERP.pam = df.pam[which(df.pam$PAM50 %in% c("LumA")),]
  dim(ERP.pam)
  
  # set.seed(seed);i = sample(dim(ERP.pam)[1],dim(ERN.pam)[1]) # take equal number of ER+ and ER- samples
  withr::with_seed(seed, {
    i = sample(dim(ERP.pam)[1],dim(ERN.pam)[1]) # take equal number of ER+ and ER- samples
  })
  
  length(ERP.pam$PatientID[i]) # ER positive samples
  length(ERN.pam$PatientID)    # ER negative samples
  
 
  mbal.pam = mat[,c(ERP.pam$PatientID[i],ERN.pam$PatientID)]
  
  dim(mbal.pam) 
  
  # Find median
  surffix = getsurffix(calibration = calibration, internal)
  
  mdns      = apply(mbal.pam,1,median,na.rm=TRUE) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.pam=mdns.df$mdns) # ER-blanced set based on IHC status alone--- 
  colnames(df.mdns) = c("X",surffix)

  ## median
  fl.mdn = BreastSubtypeRobj$medians
  
  df.al = merge(fl.mdn, df.mdns, by = "X")
  rownames(df.al) = df.al$X
  df.al = df.al[,-1]
  
  
  ## normalization
  mat= docalibration( mat, df.al, calibration, internal)
  
  out = sspPredict(BreastSubtypeRobj$centroid, mat, std=FALSE, distm="spearman", Subtype = Subtype)
  
  if(Subtype) {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.Subtype = out$predictions.Subtype, row.names = NULL )
  } else {
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
  }
  out$distances.Subtype =  -1 * out$distances.Subtype
  out$distances = -1 * out$distances
  
  ##calculate and grouping
  out.ROR = RORgroup(out, df.cln = df.pam, hasClinical = hasClinical,Subtype = Subtype )
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns.fl= df.al,outList=out))
  
}


#' Function for calling PAM50 subtypes by subgroup-specifc (ssBC) methods
#' This function is adapted from ssBC TNBC-BreastCancerRes2015 and subgroup specific TNBC-JAMAOncol2024 
#' @param mat gene expression matrix
#' @param df.cln clinical information table. The first column must be named "PatientID".
#' @param s Options are "ER" or "TN" or "ER.v2" or "TN.v2". Specify the medians you want. The original quantile is "ER" and "TN" of TNBC-BreastCancerRes2015.  If you choose "ER.v2" or "TN.v2", it means you choose quantile from TNBC-JAMAOncol2024. 
#' @param Subtype Logic. Specify whether to predict Subtype-like subtyping. 
#' @param hasClinical Logic. Specify whether clinical information is included. For example, tumor size should be in the "T" column, and lymph node status should be in the "NODE" column.
#' @noRd
#' 

makeCalls.ssBC = function(mat, df.cln, s, Subtype = FALSE , hasClinical =FALSE  ){
  
  ## loading dataset
  data("BreastSubtypeRobj")
  
  if( !(dim(mat)[2] == dim(df.cln)[1]) ){
    stop("Please input equal number of patient clinical information to the number of patient in gene expression matrix.")
  }
    
  gene.sigma = BreastSubtypeRobj$ssBC.subgroupQuantile

  if( s == "ER") { ## use ER selected strategy
    
    ## if there is no sample in either of both, wont influence the code
    ERN_samples = rownames(df.cln)[which(df.cln$ER == "ER-" )]
    ERP_samples = rownames(df.cln)[which(df.cln$ER == "ER+" )]
    samples_selected = list( ER_neg = ERN_samples, ER_pos = ERP_samples)
    
  } else if(s == "TN"){
    
    ## if there is no sample in either of both, wont influence the code
    TN_samples = rownames(df.cln)[which(df.cln$TN == "TN" )]
    samples_selected = list( TN = TN_samples)
    
  }  else if( s == "ER.v2" ){ #TNBC-JAMAOncol2024
    
    ## if there is no sample in either of both, wont influence the code
    ERN_HER2N_samples = rownames(df.cln)[which(df.cln$ER == "ER-" & df.cln$HER2 =="HER2-"  )]
    ERP_HER2N_samples = rownames(df.cln)[which(df.cln$ER == "ER+" & df.cln$HER2 =="HER2-" )]
    
    ERN_HER2P_samples = rownames(df.cln)[which(df.cln$ER == "ER-" & df.cln$HER2 =="HER2+"  )]
    ERP_HER2P_samples = rownames(df.cln)[which(df.cln$ER == "ER+" & df.cln$HER2 =="HER2+"  )]
    
    
    samples_selected = list( ERneg_HER2neg = ERN_HER2N_samples, ERpos_HER2neg = ERP_HER2N_samples,
                             HER2pos_ERneg = ERN_HER2P_samples, HER2pos_ERpos = ERP_HER2P_samples)
    
    
  } else if ( s == "TN.v2") { ## selected cohort; TNBC-JAMAOncol2024
    
    ## if there is no sample in either of both, wont influence the code
    TN_samples = rownames(df.cln)[which(df.cln$TN == "TN" )]
    samples_selected = list( TNBC = TN_samples)
    
  } else {
    stop("Please enter valid varaible for s ")
  }
  
  res = mapply(function(element){
    

    x.m = mat[, samples_selected[[element]] ]

    gene.sigma.o = gene.sigma[rownames(x.m), element]
    x.sigma = unlist(lapply(1:nrow(x.m), function(i) quantile(x.m[i,], probs = gene.sigma.o[i], na.rm = TRUE)))
    x.m = sweep(x.m, 1, x.sigma) 
    ## it has been calibrated by selected quantile 

    out = sspPredict(BreastSubtypeRobj$centroid, x.m , std=FALSE, distm="spearman", Subtype = Subtype)
    
    
  }, names(samples_selected), SIMPLIFY = FALSE, USE.NAMES = TRUE )
  
  ## prepare data for ROR
  
  if( s == "ER") { ## use ER selected strategy
    
    predictions = c(res$ER_neg$predictions, res$ER_pos$predictions)
    distances = rbind(res$ER_neg$distances, res$ER_pos$distances )
    testData = cbind(res$ER_neg$testData, res$ER_pos$testData )
    distances.Subtype = rbind(res$ER_neg$distances.Subtype, res$ER_pos$distances.Subtype )
    
    if(Subtype ){
      predictions.Subtype = c(res$ER_neg$predictions.Subtype, res$ER_pos$predictions.Subtype)
    }
    
  } else if(s == "ER.v2" ){
    
    predictions = c(res$ERneg_HER2neg$predictions, res$ERpos_HER2neg$predictions, 
                    res$HER2pos_ERneg$predictions, res$HER2pos_ERpos$predictions)
    distances = rbind(res$ERneg_HER2neg$distances, res$ERpos_HER2neg$distances,
                      res$HER2pos_ERneg$distances, res$HER2pos_ERpos$distances)
    testData = cbind(res$ERneg_HER2neg$testData, res$ERpos_HER2neg$testData,
                     res$HER2pos_ERneg$testData, res$HER2pos_ERpos$testData)
    distances.Subtype = rbind(res$ERneg_HER2neg$distances.Subtype, res$ERpos_HER2neg$distances.Subtype,
                               res$HER2pos_ERneg$distances.Subtype, res$HER2pos_ERpos$distances.Subtype)
    
    if(Subtype ){
      predictions.Subtype = c(res$ERneg_HER2neg$predictions.Subtype, res$ERpos_HER2neg$predictions.Subtype, 
                               res$HER2pos_ERneg$predictions.Subtype, res$HER2pos_ERpos$predictions.Subtype)
    }
    
    
  } else if(s == "TN"){
    
    predictions = res$TN$predictions
    distances = res$TN$distances
    testData = res$TN$testData
    distances.Subtype = res$TN$distances.Subtype
    
    if(Subtype ){
      predictions.Subtype = res$TN$predictions.Subtype
    }
    
  } else if (s == "TN.v2")  {
    
    predictions = res$TNBC$predictions
    distances = res$TNBC$distances
    testData = res$TNBC$testData
    distances.Subtype = res$TNBC$distances.Subtype
    
    if(Subtype ){
      predictions.Subtype = res$TNBC$predictions.Subtype
    }
  } 
  
  
  if(Subtype){
    out = list(predictions=predictions,predictions.Subtype = predictions.Subtype ,testData=testData,distances=distances, distances.Subtype = distances.Subtype, centroids= BreastSubtypeRobj$centroid)
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, BS.Subtype = out$predictions.Subtype , row.names = NULL )
  } else {
    out= list(predictions=predictions, testData=testData,distances=distances,  distances.Subtype = distances.Subtype, centroids=BreastSubtypeRobj$centroid)
    Int.sbs = data.frame(PatientID = names(out$predictions), BS = out$predictions, row.names = NULL )
  }
  out$distances.Subtype =  -1 * out$distances.Subtype
  out$distances = -1 * out$distances
  
  
  ##calculate and grouping
  out.ROR = RORgroup(out,df.cln, hasClinical = hasClinical, Subtype = Subtype )

  ## reorder
  orde.No = match(colnames(mat), Int.sbs$PatientID )
  Int.sbs = Int.sbs[orde.No,]
  out.ROR = out.ROR[colnames(mat),]
  out$distances.Subtype = out$distances.Subtype[orde.No,]
  out$distances = out$distances[orde.No,]

  out$testData =  out$testData[,colnames(mat) ]
  out$predictions =  out$predictions[colnames(mat)]
  
  return(list(BS.all=Int.sbs, score.ROR=out.ROR, mdns = gene.sigma, outList=out))
  
}


