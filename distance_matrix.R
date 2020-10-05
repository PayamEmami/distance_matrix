dmatrix.dissimilarity.fast <- function(x, nblocks = 5,
threads = 2, ffuse = T, SampleSize = 250, use = "pairwise.complete.obs",pathSave="partialCors/",fileuse=F, ...){
  library(tidyverse)
  require(doFuture)
  require(foreach)
  require(ff)
  require(psych)
  require(Matrix)
  library(data.table)
  # initialize parallel framework
  NCOL <- ncol(x)
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}
  ## preallocate square matrix of dimension ncol(x)
  if(!fileuse)
  {
    if (!ffuse) {corMAT <- matrix(nrow = NCOL, ncol = NCOL)}
    if (ffuse) {corMAT<-ff(vmode = "single", dim = c(NCOL, NCOL))}
  }else{
    if(!dir.exists(pathSave))
    {
      dir.create(pathSave)
    }
  }
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  column<-colnames(x)
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  plan(tweak(multiprocess,workers=threads))
  registerDoFuture()
  results <- foreach(i = 1:nrow(COMBS),.packages = c("ff","psych","data.table","Matrix")) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    COR <- cor(x[, G1, with = F], x[, G2, with = F], use = use,method = "pearson")
    if(file){
      B1 <- column[G1]
      B2 <- column[G2]
      write.table(melt(COR),file = gzfile(paste(pathSave,i,".text",sep = "")),quote = F,row.names = F,col.names = F,sep = "\t")
      flush.console() 
    }else{
      corMAT[G1, G2] <- COR
      corMAT[G2, G1] <- t(COR)
    }
    COR <- NULL
   # COR <- NULL
    B1 <- NULL
    B2 <- NULL
  }
  # 
  # stopCluster(Cluster)
  # gc()
  if(!fileuse)
  {
    if (ffuse) {corMAT <- 1 - as.ram(corMAT)
    }else{
      corMAT <- 1 - corMAT
    }
  }else{
    corMAT <-
      list.files(path =pathSave,full.names = T)   %>%
      map_df(~read_delim(gzfile(.),delim = "\t",col_names = F))  %>% dcast(.,formula = X1~X2, value.var='X3')
    rownames(corMAT)<-corMAT[,1]
    corMAT<-corMAT[,-1]
    corMAT<-corMAT[column,column]
    corMAT[lower.tri(corMAT)] <- t(corMAT)[lower.tri(corMAT)]
    corMAT <- 1 - corMAT
  }
  return(corMAT)
}
