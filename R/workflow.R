


.getBestComb <- function(species,preds,n=1,comb,bgn=1000) {
  # The following two lines are not needed if we sure the input is presence-only with "species" column
  # but to make sure, I just re-generated the SpatVector!!
  .sp <- vect(crds(species)) # make a SpatVector with coordinates of the presence locations
  .sp$species <- 1 # add species column
  
  .auc <- rep(0,length(comb))
  #------
  for (i in seq_along(comb)) {
    .w <- comb[[i]] # which variables/layers combination in item i of the comb list?
    .prs <- preds[[.w]] # selected variable combinations!
    d <- sdmData(species~., .sp, .prs, bg=list(method='gRandom',n=bgn))
    m <- sdm(species~., d, methods = c('glmp','brt','rf','svm','cart','maxent','bioclim.dismo','mlp','mda'),
             n=n,replication='boot')
    e <- .getEval(m,n='species',setting=list(method='weighted',stat='auc')) # from model -> get Evaluation based on Ensemble
    .auc[i] <- e@statistics$AUC
  }
  
  cat('\n The top 10 combinations (sorted): ',paste(order(.auc,decreasing = T),collapse = ', '))
  
  .auc
}
#-------


# The core function in Workflow (can be used in a more efficient parallelisation!)
.coreWF <- function(i,...) {
  .w <- comb[[i]] # which variables/layers combination in item i of the comb list?
  .prs <- preds[[.w]] # selected variable combinations!
  d <- sdmData(species~., .sp, .prs, bg=list(method = 'gRandom',n=bgn))
  m <- sdm(species~., d, methods = c('glmp','brt','rf','svm','cart','maxent','bioclim.dismo','mlp','mda'),
           n=n,replication='boot')
  #----
  if (path == '.') write.sdm(m, filename = paste0('M_',i,'_',spn,'_',filename,'.sdm'))
  else write.sdm(m, filename = paste0(path,'/M_',i,'_',spn,'_',filename,'.sdm'))
  #----------
  # Save predictions:
  if (path == '.') .pr <- predict(m, .prs, filename = paste0('predicts/pr_COMB-',i,'_',spn,'_',filename,'.tif'))
  else .pr <- predict(m, .prs, filename = paste0(path,'/predicts/pr_COMB-',i,'_',spn,'_',filename,'.tif'))
  
  #------
  
  if (path == '.') en <- ensemble(m, .pr, setting=list(method='weighted',stat='auc'),filename = paste0('ensembles/en_COMB-',i,'_',spn,'_',filename,'.tif'))
  else en <- ensemble(m, .pr, setting=list(method='weighted',stat='auc'),filename = paste0(path,'/ensembles/en_COMB-',i,'_',spn,'_',filename,'.tif'))
  #---------
  rm(m); gc()
  TRUE
}
#---------
# error-handling of WF:
.wf_ER <- function(i,...) {
  e <- try(.coreWF(i,...),silent=TRUE)
  if (inherits(e,'try-error')) FALSE
  else TRUE
}
#--------------



runWF <- function(species,preds,n=10,comb,bgn=1000,path='.',filename,ncore=8,spn) {
  if (!is.null(ncore)) {
    ncore <- c(parallel::detectCores()-2,ncore)
    ncore <- ncore[which.min(ncore)]
  }
  
  
  # The following two lines are not needed if we sure the input is presence-only with "species" column
  # but to make sure, I just re-generated the SpatVector!!
  .sp <- vect(crds(species)) # make a SpatVector with coordinates of the presence locations
  .sp$species <- 1 # add species column
  
  if (!is.null(ncore) && ncore > 1) {
    # we need to save the files, and read them in parallel sessions (revise the code if you already have files)
    writeVector(.sp,'_temp_species_file.shp',overwrite=T)
    writeRaster(preds,filename='_temp_predictors_file.tif',overwrite=T)
    
    require(parallel)
    cl <- makeCluster(ncore)
    clusterEvalQ(cl, {
      library(terra)
      library(sdm)
      getmethodNames()
      .sp <- vect('_temp_species_file.shp')
      preds <- rast('_temp_predictors_file.tif')
      NULL
    })
    clusterExport(cl, c("path","filename","spn","comb","n","bgn",".coreWF"))
    
    o <- parLapply(cl,1:length(comb),.wf_ER)
    stopCluster(cl)
    gc()
    
    return(unlist(o))
    
  } else {
    sapply(cl,1:length(comb),.wf_ER)
  }
  
}

