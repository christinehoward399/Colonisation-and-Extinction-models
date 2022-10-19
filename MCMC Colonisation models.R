##~~~~~~~~~~~~~~~~~~~ ##
##  COLONISATION MODELS
## MCMCglmms - want to account for both spatial and phylogenetic autocorrelation
##~~~~~~~~~~~~~~~~~~~ ##

rm(list=ls())
setwd("X:\\EBBA2")

library(stringr)
library(MCMCglmm)
library(reshape2)
library(phytools)
library(parallel)

results.dir.path.main<-"X:\\EBBA2\\Model Outputs\\MCMC Colonisation models"
dir.create(results.dir.path.main)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DATA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## required data frames 
cell # cell specific data including land cover, sampling effort etc. 
species ## species trait data

### ~~~~~~~~~~~~~~~~~~ MODEL FITTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Need to fit to 10 different absence samples

datFile<-"X:\\EBBA2\\Data\\Colonisation Random Samples"
its<-list.files(datFile)
its

for (DatI in its){
  
  iNum<-gsub(".csv","",strsplit(DatI," ")[[1]][3])
  
  ## SAMPLED CELLS
  Colonisation<-read.csv(paste0(datFile,"\\",DatI))
  
  ## Need to add in the climate suitability for those cells
  Colon<-do.call(rbind, lapply(unique(Colonisation$scientific_name), function(x){
    print(x)
    ## reduce down to the species
    subDat<-Colonisation[Colonisation$scientific_name==x,]
    ## Change in climate suitability
    climDat<-read.table(paste0("X:\\EBBA2\\Climate Suitability Modelling\\SDM results\\EBBA\\",x,"\\",x," climate prediction summary.txt", sep=""))
    try(subDat<-merge(subDat, climDat[,c("CGRSNAME", "pred1970","pred1995","ClimChange")], 
                      by.x="cell50x50", by.y="CGRSNAME", all.x=T), silent=T)
    return(subDat)
  }))
  
  # Format remaining data
  MyDat<-merge(cells, Colon[,c("cell50x50","scientific_name","state","point.dist","COG.dist","pred1970","ClimChange")],#
               by.x="CGRSNAME", by.y="cell50x50", all=T)
  
  ## Add in species traits
  MyDat<-merge(MyDat, species[,c("Scientific.Name","Mean.Mass","Generation.Length",
                                 "Clutch","Migratory.Status","rangeSize","HWI","MigDist",
                                 "MAES.habitat","Policy","HabitatBreadth","dietBreadth",
                                 "Montane","Persecuted","TipLabel")],
               by.x="scientific_name", by.y="Scientific.Name")#, all.x=T)
  MyDat<-droplevels(MyDat)
  
  ## For each species want the current extent of of favourable land cover and change in extent of favourable land cover
  landAsc<-data.frame(MAES=levels(as.factor(MyDat$MAES.habitat) ) )
  landAsc$ESAgrp<-c("Croplnd","Croplnd","Grsslnd","Shrblnd","Wtr_bds","Sprsly_","Urban","Wetland","Forest")
  
  MyDat$CurrEoFLC<-0
  MyDat$ChangeEoFLC<-0
  MyDat<-MyDat[complete.cases(MyDat$MAES.habitat),]
  
  for (x in 1:length(MyDat$CurrEoFLC)){
    print(x)
    LC<-landAsc[landAsc$MAES==as.character(MyDat[x,"MAES.habitat"]),"ESAgrp"]
    MyDat$ChangeEoFLC[x]<-MyDat[x,paste0(LC,"_change")]
    MyDat$CurrEoFLC[x]<-MyDat[x,paste0(LC,"_1992")]
  }
  
  ## Format final data frame for modelling
  MyDat$animal<-MyDat$TipLabel
  data.model<-MyDat[,c("scientific_name","CGRSNAME","lon","lat","point.dist","rangeSize","pred1970","COG.dist",#"code",
                       "ClimChange","NAME_FAO","Migratory.Status", "alt_range","IUCN",#"alt_mean",,
                       "Generation.Length","Mean.Mass","Clutch","Policy","CurrEoFLC","ChangeEoFLC", "HWI","MigDist",
                       "Shannon_change","Shannon_1992","totalConfidence","HabitatBreadth","dietBreadth","Persecuted", 
                       "MAES.habitat","Montane","state","animal")]
  
  ## Subset down to cells with similar sampling effor between the two years
  data.model<-data.model[data.model$totalConfidence==2,]
  formula <- as.formula(state~ClimChange+point.dist+CurrEoFLC+Mean.Mass+Policy+Montane+COG.dist+
                          MigDist+Generation.Length+alt_range+IUCN+HWI+
                          HabitatBreadth+dietBreadth+Clutch+
                          ChangeEoFLC+pred1970+Persecuted+
                          Shannon_1992+rangeSize)
  data.model<-data.model[complete.cases(data.model$scientific_name),]
  data.model<-data.model[complete.cases(data.model$HabitatBreadth),]
  data.model<-data.model[complete.cases(data.model$ClimChange),]
  data.model<-droplevels(data.model)
  missSp<-unique(data.model[is.na(data.model$animal),"scientific_name"])
  data.model$animal<-ifelse(is.na(data.model$animal), data.model$scientific_name, data.model$animal)
  data.model$animal<-gsub(" ","_",data.model$animal)
  data.model<-droplevels(data.model)
  data.model$Montane<-as.factor(as.character(data.model$Montane))
  data.model$Policy<-as.factor(as.character(data.model$Policy))
  data.model$Persecuted<-as.factor(as.character(data.model$Persecuted))
  
  ## Standardised variables
  cont<-c("ClimChange","Mean.Mass","point.dist","Generation.Length","COG.dist",
          "alt_range","pred1970","IUCN","MigDist","HWI",
          "Clutch","CurrEoFLC","ChangeEoFLC","Shannon_change","Shannon_1992",
          "rangeSize","HabitatBreadth","dietBreadth")
  for ( i in cont){data.model[,i]<-scale(as.numeric(as.character(data.model[,i])))}
  
  ## create folder for output
  results.dir.path.sub<-paste0(results.dir.path.main,"\\Final")
  dir.create(results.dir.path.sub)
  
  ## Set up parallel
  clust <- makeCluster(10)
  clusterExport(clust, varlist = c("trlist","missSp","data.model","formula","results.dir.path.sub","DatI","iNum"))
  clusterEvalQ(clust, library(phytools))
  clusterEvalQ(clust, library(MCMCglmm))
  
  system.time(parLapply(clust, 1:10,function(x){
    # x<-1
    
    ConsensusTree<-trlist[[x]]
    
    ## ~~~~~~ FORMAT THE CONSENSUS TREE
    ## Add missing species
    for(k in missSp){
      ConsensusTree<-add.species.to.genus(ConsensusTree,k)
    }
    
    ## trim the consensus tree
    tr<-drop.tip(ConsensusTree,ConsensusTree$tip.label[-match(unique(data.model$animal),ConsensusTree$tip.label)])
    # Prune tree to just include species of interest
    
    # Remove node names from phylogeny
    tr$node.label <- NULL
    # MCMCglmm needs ultrametric trees
    if(is.ultrametric(tr)==FALSE){
      tr <- chronos(tr)
      class(tr) <- "phylo"
    }
    
    # check again
    is.ultrametric(tr)
    
    ## Define the covariance matrix
    inv.phylo<-inverseA(tr,nodes="TIPS",scale=TRUE)
    
    ## ~~~~~~~~~~~~~~~~~~~ Model Fitting M1 ~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## params for sample size etc.
    THIN <- 200
    BURNIN <- 20000
    NITT <-220000
    
    # PRIORS- Model M1 is a simple phylogenetic model without intraspecific correlation structure
    priorpr.m1 <- list(R = list(V = 1, nu = 1),#0.002),
                       G = list(G1 = list(V = 1, nu = 0.002),
                                G2 = list(V = 1, nu = 0.002),
                                G3 = list(V = 1, nu = 0.002),
                                G4 = list(V = 1, nu = 0.002)))
    
    #Then you can run your global model: 
    system.time(mmF<- MCMCglmm(formula,#+X2010.IUCN.Red.List.category,
                               random=~animal+idv(NAME_FAO)+idv(scientific_name)+idv(CGRSNAME),#
                               family="categorical",ginverse=list(animal=inv.phylo$Ainv),prior=priorpr.m1,
                               data=data.model,nitt=NITT,burnin=BURNIN,thin=THIN, singular.ok = T))# 
    
    ## ~~~~~~~~~~~~~~~~~~~ Model Performance M1 ~~~~~~~~~~~~~~~~~~~~~~~ ##
    # Calculation of the variance in fitted values
    mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))
    # MCMCglmm - marginal
    vmVarF<-numeric(dim(mmF$Sol)[1])
    for(i in 1:length(vmVarF)){
      Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
      vmVarF[i]<-Var}
    R2m<-vmVarF/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3]+mmF$VCV[,4]+mmF$VCV[,5])
    margR2M1<-mean(R2m);margR2M1
  
    mod<-list(mod1=mmF,MarginalR2M1=margR2M1)
    ##Save block models
    save(mod, file=paste(results.dir.path.sub,"\\Colonisation_MCMCglmm_iteration_",iNum,"_phylo_",x,".rda",sep=""),compress="bzip2")
  }))
  
}


plot(M1)
library(coda)
gelman.plot(x=mcmc.list(M1$Sol, M2$Sol))
gelman.diag(mcmc.list(M1$Sol, M2$Sol))
 
library(BayesianTools)
correlationPlot(data.frame(mmF$Sol)[,1:10])
correlationPlot(data.frame(mmF$Sol)[,11:19])
correlationPlot(data.frame(mmF$Sol))
