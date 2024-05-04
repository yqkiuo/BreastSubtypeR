
## create ssBC Rdata

# toolpath = "/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/2.projts/0.KI/2.Rtools/"
# 
# load(paste(toolpath,"R_code_ssBC/","ssBCsubtyping.Rdata", sep = ""))
# 
# rm(pam50.symbol2symbol.v2)
# ## edit llaply to lapply pam50.symbol2symbol.v2.R
# path = paste(toolpath,"R_code_ssBC/","pam50.symbol2symbol.v2.R", sep = "")
# 
# rm(toolpath)
# save.image(file = paste("/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/ProjectsAtKI/6.PAM50/data/" ,
#                   "ssBCsubtyping.RData",sep = "" ))


## create PCA-PAM50 Rdata
# toolpath = "/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/2.projts/0.KI/2.Rtools/"
# 
# PAM50dir  = paste( toolpath,"R_code_PCA-PAM50","PCA-PAM50", "PAM50/bioclassifier_R", sep = "/" ) 
# fl.nm     = paste(PAM50dir,"mediansPerDataset_v2.txt",sep="/")
# fl.mdn    = read.table(fl.nm,sep="\t",header=T,stringsAsFactors=F)
# 
# 
# trainCentroids<- paste(PAM50dir,"pam50_centroids.txt",sep="/")#_PNmdfd
# # load the published centroids for classifcation
# pamout.centroids<-read.table(trainCentroids,sep="\t",header=T,row.names=1)
# 
# trainFile<- paste(PAM50dir,"220arrays_nonUBCcommon+12normal_50g.txt",sep="/")
# 
# x<-readarray(trainFile,hr=2)
# 
# PCA_PAM50 = list(fl.mdn = fl.mdn, pamout.centroids = pamout.centroids,x = x )
# save(PCA_PAM50,
#      file = paste("/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/ProjectsAtKI/6.PAM50/data/" , "PCA-PAM50.RData",sep = "" ) )




#### Collect geneID ####
# 
# library(stringr)
# library(dplyr)
# 
# library(PAM50subtyping)
# library(org.Hs.eg.db)
# 
# library(AIMS)
# genes = str_extract_all(AIMSmodel$all.pairs,"\\d+" )
# 
# Genes_Used_AIMS = data.frame(genes_AIMS = unique( unlist( genes)))

# 
# data("PCA-PAM50")
# 
# pathParker50genes = "/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/2.projts/0.KI/2.Rtools/R_code_parkerPAM50/PAM50_R"
# 
# genes_anno = read.table( paste(pathParker50genes,"pam50_annotation.txt",sep = "/" ), sep = "\t", header = T)
# 
# library(sspbc)
# library(stringr)
# 
# genes_pam50 = str_extract_all(sspbc.models$ssp.pam50$all.pairs, "\\d+" )
# genes_subtype = str_extract_all(sspbc.models$ssp.subtype$all.pairs, "\\d+" )
# Genes_Used_sspbc = data.frame(genes_sspbc = unique( c (unlist( genes_pam50) ,unlist(genes_subtype ) ) ))

# #
# ## combine all together and get all possible annotations
# 
# genes_all_entrezID =  data.frame( IDs = unique( c( genes_anno$EntrezGene, Genes_Used_AIMS$genes_AIMS, Genes_Used_sspbc$genes_sspbc) ) )
# 
# 
# res = select(org.Hs.eg.db, keys =genes_all_entrezID$IDs, columns = c( "SYMBOL","ALIAS" ), keytype='ENTREZID' , multiVals = "CharacterList")
# 
# 
# res_collapsed <- res %>%
#   group_by(ENTREZID, SYMBOL) %>%
#   summarise(ALIAS = paste(ALIAS, collapse = ", "))
# 
# ## save as a text for manual check
# pathres="/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/ProjectsAtKI/6.PAM50/test"
# write.table(res_collapsed, file = paste(pathres, "GeneID_mapping.text",sep = "/" ), sep = "\t",row.names = F )
# 



#### create new data set for our R package ####
## centroid 
## medians
## gene.signature for parker-based
## all gene.signature 

# 
# 
# path = "/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/6.PAM50/test"
# genes = read.table( paste( path, "GeneID_mapping.txt", sep = "/"), sep = "\t", header = T)
# genes$Alias =ifelse( genes$Supplements == "",genes$Alias, paste0(genes$Alias, ", ", genes$Supplements )  )
# 
# 
# centroid = PCA_PAM50$pamout.centroids
# medians = PCA_PAM50$fl.mdn ## if have the latest median table??
# 
# ## gene.50
# genes.sig50 = intersect( rownames(centroid), medians$X)
# zid = genes$EntrezGene.ID
# names(zid) = genes$Symbol
# 
# genes.sig50.tl = data.frame( Symbol = genes.sig50, EntrezGene.ID = '' )
# genes.sig50.tl$EntrezGene.ID = zid[match( genes.sig50.tl$Symbol , names(zid) )  ]
# 
# genes.sig50.tl$EntrezGene.ID[is.na(genes.sig50.tl$EntrezGene.ID)] = lapply(genes.sig50.tl$Symbol[is.na(genes.sig50.tl$EntrezGene.ID)], function(x){
#   y = genes$EntrezGene.ID[  grepl(x, genes$Alias ) ]
#   return(y)
# })
# 
# ## proliferation genes
# proliferationGenes<-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
# proliferationGenes.prosigna = c("ANLN", "CEP55", "ORC6L", "CCNE1", "EXO1", "PTTG1", "CDC20", "KIF2C", "RRM2", "CDC6", "KNTC2", "TYMS", "CDCA1",
#                                 "MELK", "UBE2C", "CENPF", "MKI67", "UBE2T")
# genes.sig50.tl$Inproliferation = ifelse( genes.sig50.tl$Symbol %in% proliferationGenes, "yes", "no"  )
# genes.sig50.tl$Inproliferation.prosigna = ifelse( genes.sig50.tl$Symbol %in% proliferationGenes.prosigna, "yes", "no"  )
# genes.sig50 = genes.sig50.tl
# genes.sig50$EntrezGene.ID = unlist( genes.sig50$EntrezGene.ID)
# 
# ## ssBC
# path="/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/2.projts/0.KI/2.Rtools/R_code_ssBC"
# load(file = paste(path, "/","ssBCsubtyping.Rdata",sep = "" ))
# 
# ## reading latest ssBC ER/HER2 balanced quantile dataset
# path = "/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/6.PAM50/test"
# SIGMA = read.table( file = paste( path, "SIGMA.txt", sep = "/" ) , header = TRUE)
# 
# N_inter = intersect( rownames(SIGMA ), rownames(Pam50$subgroupQuantile))
# length(N_inter)
# 
# ##ordering
# SIGMA = SIGMA[rownames(Pam50$subgroupQuantile) ,]
# ## cbind
# SIGMA_ssBC = cbind( Pam50$subgroupQuantile, SIGMA)
# 
# ## create IBC.parker and save
# IBC.parker = list(medians = medians, centroid = centroid, genes.sig50 = genes.sig50, platfrom = Pam50$platform, subtypeCode = Pam50$subtypeCode, ssBC.subgroupQuantile = SIGMA_ssBC,  UNC232 = UNC232 )
# save(IBC.parker,
#      file = paste("/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/ProjectsAtKI/6.PAM50/data/" , "IBC.parker.RData",sep = "" ) )
# 
# 
# ## all genes.signature
# path = "/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/6.PAM50/test"
# genes = read.table( paste( path, "GeneID_mapping.txt", sep = "/"), sep = "\t", header = T)
# genes$Alias =ifelse( genes$Supplements == "",genes$Alias, paste0(genes$Alias, ", ", genes$Supplements )  )
# 
# genes.signature = genes[,c("EntrezGene.ID","Symbol" ,"Alias" )]
# 
# library(AIMS)
# genes = str_extract_all(AIMSmodel$all.pairs,"\\d+" )
# 
# Genes_Used_AIMS = data.frame(genes_AIMS = unique( unlist( genes)))
# 
# library(sspbc)
# library(stringr)
# 
# genes_pam50 = str_extract_all(sspbc.models$ssp.pam50$all.pairs, "\\d+" )
# genes_subtype = str_extract_all(sspbc.models$ssp.subtype$all.pairs, "\\d+" )
# Genes_Used_sspbc = data.frame(genes_sspbc = unique( c (unlist( genes_pam50) ,unlist(genes_subtype ) ) ))
# 
# 
# genes.signature$AIMS_based = ifelse( genes.signature$EntrezGene.ID %in% unique( c(Genes_Used_AIMS$genes_AIMS,Genes_Used_sspbc$genes_sspbc ) ),"Yes", "No" )
# genes.signature$AIMS = ifelse( genes.signature$EntrezGene.ID %in% unique( c(Genes_Used_AIMS$genes_AIMS) ),"Yes", "No" )
# 
# genes.signature_split = separate_rows(genes.signature, Alias, sep = ", ")
# 
# save(genes.signature,
#      file = paste("/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/ProjectsAtKI/6.PAM50/data/" , "genes.signature.RData",sep = "" ) )


