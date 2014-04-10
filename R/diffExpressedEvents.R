.readAndPrepareData <- function(countsData,conditions) {
  ###################################################
  ### code chunk number 1: Read data
  ###################################################
  sortedconditions <- sort(conditions)
  n <- length(unique(sortedconditions))
  nr <- rle(sortedconditions)$lengths
  sortedindex <- order(conditions)+2
  namesData <- c("ID","Length",rep(NA,length(conditions)))
  for (k in 1:nr[1]){
    namesData[2+k] <- paste(sortedconditions[k],"_r",k,sep="",collapse="")
  }
  for (i in 2:n) {
    for (j in 1:nr[n]) {
      namesData[2+cumsum(nr)[n-1]+j] <- paste(sortedconditions[cumsum(nr)[n-1]+j],"_r",j,sep="",collapse="")
    }
  }
  countsData[,-(1:2)] = countsData[,sortedindex]
  colnames(countsData) <- namesData
  options(warn=-1) # suppress the warning for the users
  countsData$Path <- gl( 2, 1, dim(countsData)[1], labels = c("UP", "LP"))

  ###################################################
  ### code chunk number 2: Normalisation
  ###################################################
  # Normalisation with DESeq
  conds <- c()
  for( i in 1:n ) {
    for( j in 1:nr[i] ) {
      conds <- c( conds,paste( "Cond", i, sep = "",collapse = "") )
    }
  } 
  cds <- newCountDataSet( countsData[ ,3:(3+length(conds)-1)], conds ) # create object
  cdsSF <- estimateSizeFactors(cds)
  sizeFactors( cdsSF )
  shouldWeNormalise=sum(is.na(sizeFactors(cdsSF))) < 1
  dim <- dim(countsData)[2]
  countsData[ ,(dim+1):(dim+length(conds)) ] <- round(counts(cdsSF, normalized=shouldWeNormalise))
  colnames(countsData)[(dim+1):(dim+length(conds))] <- paste(namesData[3:(3+sum(nr)-1)],"_Norm",sep="")
  return(list(countsData, conds, dim, n, nr, sortedconditions))
}

.eventtable <- function(df,startPosColumn4Counts, endPosCol4Counts){
  eventTab = data.frame(ID=rep(as.factor(df['ID']), endPosCol4Counts-startPosColumn4Counts+1),
  cond=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[2]}))),
  counts=as.numeric(df[startPosColumn4Counts:endPosCol4Counts]),
  path=as.factor(unlist(lapply(strsplit(names(df)[startPosColumn4Counts:endPosCol4Counts],"_"),FUN=function(d){d[1]}))),
  row.names = NULL)
  return(eventTab)
}

.fitNBglmModelsDSSPhi <- function(eventdata, phiDSS, phiDSScond, phiGlobal, nbAll){
    # S: simple, A: additive, I : interaction models
    # Poisson model  
  nbglmS0 <- negbin(counts~cond, data=eventdata, random=~1, fixpar=list(3,0))
  nbglmA0 <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4,0))
  nbglmI0 <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5,0))  

  nbAnov0 <- anova(nbglmS0, nbglmA0, nbglmI0)
  nbAIC0 <- c(AIC(nbglmS0,k = log(nbAll))@istats$AIC, AIC(nbglmA0,k = log(nbAll))@istats$AIC, AIC(nbglmI0,k = log(nbAll))@istats$AIC)
    # singular.hessian:  true when fitting provided a singular hessian, indicating an overparamaterized model.
  nbSingHes0 <- c(nbglmS0@singular.hessian, nbglmA0@singular.hessian, nbglmI0@singular.hessian)
    # code: ‘code’ An integer (returned by ‘optim’) indicating why the optimization process terminated.
  nbCode0 <- c(nbglmS0@code, nbglmA0@code, nbglmI0@code)  
  
    # binomial negative model, with global phi
  nbglmSgb <- negbin(counts~cond, data=eventdata, random=~1, fixpar=list(3,phiGlobal))
  nbglmAgb <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4,phiGlobal))
  nbglmIgb <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5,phiGlobal)) 
  
  nbAnovgb <- anova(nbglmSgb, nbglmAgb, nbglmIgb)
  nbAICgb <- c(AIC(nbglmSgb,k = log(nbAll))@istats$AIC, AIC(nbglmAgb,k = log(nbAll))@istats$AIC, AIC(nbglmIgb,k = log(nbAll))@istats$AIC)
    # the BIC in fact, since we use k = log(nobs)
  nbSingHesgb <- c(nbglmSgb@singular.hessian, nbglmAgb@singular.hessian, nbglmIgb@singular.hessian)
  nbCodegb <- c(nbglmSgb@code, nbglmAgb@code, nbglmIgb@code)
    
    # binomial negative model, with phi DSS
  nbglmS <- negbin(counts~cond, data=eventdata, random=~1, fixpar=list(3, phiDSS))
  nbglmA <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4, phiDSS))
  nbglmI <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5, phiDSS))
  
  nbAnov <- anova(nbglmS, nbglmA, nbglmI)
  nbAIC <- c(AIC(nbglmS,k = log(nbAll))@istats$AIC, AIC(nbglmA,k = log(nbAll))@istats$AIC, AIC(nbglmI,k = log(nbAll))@istats$AIC)
  nbSingHes <- c(nbglmS@singular.hessian, nbglmA@singular.hessian, nbglmI@singular.hessian)
  nbCode <- c(nbglmS@code, nbglmA@code, nbglmI@code)
  
    # binomial negative model, with phi DSS, conditionally to the expression mean  
  nbglmScond <- negbin(counts~cond, data=eventdata, random=~1, fixpar=list(3, phiDSScond))
  nbglmAcond <- negbin(counts~cond + path, data=eventdata, random=~1, fixpar=list(4, phiDSScond))
  nbglmIcond <- negbin(counts~cond * path, data=eventdata, random=~1, fixpar=list(5, phiDSScond))

  nbAnovcond <- anova(nbglmScond, nbglmAcond, nbglmIcond)
  nbAICcond <- c(AIC(nbglmScond,k = log(nbAll))@istats$AIC, AIC(nbglmAcond,k = log(nbAll))@istats$AIC, AIC(nbglmIcond,k = log(nbAll))@istats$AIC)
  nbSingHescond <- c(nbglmScond@singular.hessian, nbglmAcond@singular.hessian, nbglmIcond@singular.hessian) 
  nbCodecond <- c(nbglmScond@code, nbglmAcond@code, nbglmIcond@code) 
  
  rslts <- c(nbAnov0@anova.table$'P(> Chi2)'[2:3],
  nbAnovgb@anova.table$'P(> Chi2)'[2:3],
  nbAnov@anova.table$'P(> Chi2)'[2:3],
  nbAnovcond@anova.table$'P(> Chi2)'[2:3],
             nbAIC0,
             nbAICgb, 
             nbAIC,
             nbAICcond,
             
             nbCode0,
             nbCodegb,
             nbCode,
             nbCodecond,
             
             nbSingHes0,
             nbSingHesgb,
             nbSingHes,
             nbSingHescond)
  return(rslts)  
}

kissplice2counts <- function(fileName) {
  toConvert <- file(fileName, open = "r")
  line <- readLines(toConvert)
  index <- 1
  indexNames <- 1
  if (substr(line[index], start = 0, stop = 1) == '>') {
    len <- length(strsplit(line[index], "|", fixed = TRUE)[[1]])-6
  }
  events.mat <- matrix(NA,length(line)/2,len+2)
  events.names <- rep(NA,length(line)/2)
  firstLineChar <- substr(line[index], start = 0, stop = 1)
  if (firstLineChar == '>') {#checks if the line contains the header beginning with ">", not the sequence
    while (index <= length(line)) {
      nbCharLine <- length(strsplit(line[index],"|")[[1]]) # number of characters in the line
      lineInfos <- substr(line[index],start = 2,stop =nbCharLine) #the line without ">" that we want to avoid
      lineSplit <- strsplit(lineInfos, "|", fixed=TRUE)[[1]] # gets pieces of information separated by "|" in KisSplice format
      #example of lineSplit :  
      # lineSplit[1] : "bcc_3929"
      #          [2] : "Cycle_0"
      #          [3] : "Type_1"     
      #          [4] : "upper_path_length_90" 
      #          [5] : "C1_1" (first condition)
      #          [5+len] : "Cn_5" (last condition)              
      #          [5+len+1] : rank_0.90267" 
      eventName <- paste(lineSplit[1],lineSplit[2],sep="_") # concatenates the two first elements of lineSplit, ie the event name
      events.names[indexNames] <- eventName
      lengthInfo <- lineSplit[4]
      events.mat[indexNames,1] <- as.numeric(strsplit(lengthInfo,"_")[[1]][4]) #fills the first column of the matrix with length info
      lineCounts <- lineSplit[5:(5+len)] # gets every condition of the line and its associated count
      events.mat[indexNames,2:dim(events.mat)[2]] <- as.numeric(lapply(lineCounts, function(x) strsplit(x, "_")[[1]][2])) #fills the matrix others columns with the counts of the conditions of the line 
    index <- index + 2
    indexNames <- indexNames + 1
    }
  }
  class(events.mat) <- "numeric"
  events.df <- as.data.frame(events.mat)
  events.df <- data.frame(events.names,events.df)
  return (events.df)
}

qualityControl <- function(countsData,conditions) {

  ###################################################
  ### code chunk number 1: Read and prepare data
  ###################################################
  list_Data <-.readAndPrepareData(countsData,conditions)
  countsData <- list_Data[[1]]
  conds <- list_Data[[2]]
  dim <- list_Data[[3]]
  n <- list_Data[[4]]
  nr <- list_Data[[5]]

  ###################################################
  ### code chunk number 2: fig_hclust_norm
  ###################################################
  plot(hclust(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))])),"ward"))
  par(ask=TRUE)

  ###################################################
  ### code chunk number 3: replicates
  ###################################################
  heatmap(as.matrix(as.dist(1-cor(countsData[ ,(dim+1):(dim+length(conds))]))), margins = c(10,10))

  ###################################################
  ### code chunk number 4: intra-group and inter-group-variance
  ###################################################
  ## Mean and variance over all conditions and replicates (normalized counts!) 
  countsData$mn <- apply(countsData[ ,(dim+1):(dim+length(conds))], 1, mean)
  countsData$var <- apply(countsData[ ,(dim+1):(dim+length(conds))], 1, var)
  ## correction term
  nbAll <- sum(nr) # number of all observations in all groups
  countsData$ct  <- apply(countsData[ ,(dim+1):(dim+length(conds))], 1, sum)^2/nbAll
  ## sum of squares between groups
  countsData$ss <- apply(countsData[ ,(dim+1):(dim+length(conds)/n)],1,sum)^2/nr[1] + 
  apply(countsData[ ,((dim+1)+length(conds)/n):(dim+length(conds))],1,sum)^2/nr[2] 
  ## substract the correction term from the SS and divide by the degrees of 
  df <- 1 # freedom(groups); here: df=2-1=1
  countsData$varInter <- (countsData$ss - countsData$ct)/df
  # intra-variability 
  countsData$varC1 <- apply( countsData[ ,(dim+1):(dim+nr[1])], 1, var )
  countsData$varC2 <- apply( countsData[ ,((dim+1)+nr[1]):(dim+nr[2]+nr[1])], 1, var )
  countsData$varIntra <- apply( data.frame(countsData$varC1, countsData$varC2), 1, mean )

  ###################################################
  ### code chunk number 5: intra-vs-inter
  ###################################################
  plot( x = countsData$varIntra, y = countsData$varInter, xlab = "Intra-variability", ylab = "Inter-variability", las = 1, log = "xy")
  abline( a = 0, b = 1, col = 2, lty = 2, lwd = 2 )
}

diffExpressedVariants <- function(countsData,conditions) {

  ###################################################
  ### code chunk number 1: Read and prepare data
  ###################################################
  list_Data <-.readAndPrepareData(countsData,conditions)
  countsData <- list_Data[[1]]
  n <- list_Data[[4]]
  nr <- list_Data[[5]]
  sortedconditions <- list_Data[[6]]

  ##################################################
  ## code chunk number 2: event-list
  ##################################################
  # reduce data frame to the interesting columns
  nbAll <- sum(nr)
  dataPart <- countsData[ ,c(1:2,which(grepl("_Norm",names(countsData)))) ] # -> when not want to filter data 
  dataPart$Path <- gl( 2,1,dim(countsData)[1],labels=c("UP","LP") ) 
  dataPart2 <- cbind( dataPart[ seq( 1, dim( dataPart )[1], 2 ), ], dataPart[ seq( 2, dim(dataPart)[1] , 2 ), grepl("Norm", names( dataPart) ) ] )
  names(dataPart2)[3:(3+nbAll-1)] <- paste( "UP",names(dataPart2)[3:(3+nbAll-1)], sep="_" )
  names(dataPart2)[(3+nbAll+1):(3+2*nbAll+1-1)] <- paste( "LP", names(dataPart2)[(3+nbAll+1):(3+2*nbAll+1-1)], sep="_" )
  dataPart2 <- dataPart2[ ,c(-(3+nbAll))]
  if (anyDuplicated(dataPart2[ ,1]) > 0) {
    dataPart2 <- dataPart2[!duplicated(as.character(dataPart2[ ,1])), ]
  }
  rownames(dataPart2)=as.character(dataPart2[ ,1])

  # create list for the complete data set
  allEventtables  <- apply(dataPart2,1,.eventtable, startPosColumn4Counts = which(grepl("UP",names(dataPart2)))[1],endPosCol4Counts = ncol(dataPart2))

  ###################################################
  ### code chunk number 3: DSS dispersion estimation
  ###################################################
  dataNormCountsEvent <- as.matrix(dataPart2[ ,3:ncol(dataPart2)]) # the counts matrix
  colnames(dataNormCountsEvent) <- 1:ncol(dataNormCountsEvent)
  designs <- rep(c(1:(n*2)),c(nr,nr)) # the design matrix
  dispData <- newSeqCountSet(dataNormCountsEvent, as.data.frame(designs))
  dispData <-estDispersion(dispData)
  dispDataMeanCond <- newSeqCountSet(dataNormCountsEvent, as.data.frame(designs))
  # if (n > 2) {
  #     dispDataMeanCond <- estDispersion(dispData)
  #     } else {
  #       dispDataMeanCond <- estDispersion(dispData,trend=T)
  #     }
  dispDataMeanCond <- estDispersion(dispData)

  ###################################################
  ### code chunk number 4: variance - mean - Event level1
  ###################################################
  # compute mean and variance per Event (instead of per allele)
  event.mean.variance.df <- as.data.frame(cbind(apply(dataPart2[ ,which(grepl("_Norm",names(dataPart2)))],1,mean),apply(dataPart2[ ,which(grepl("_Norm",names(dataPart2)))],1,var)))
  names(event.mean.variance.df) <- c("Mean", "Variance")
  rownames(event.mean.variance.df) <- as.character(dataPart2[ ,1])
  # estimate the dispersion parameter D of the Quasi-Poisson distribution
  lm.D <- lm(event.mean.variance.df$Variance ~ event.mean.variance.df$Mean-1)
  # coef(lm.D)
  ## estimate the overdispersion parameter theta of the NB distritution
  modelNB <- Variance ~ Mean + 1/theta * Mean^2
  nls.modelNB <- nls(modelNB, data=event.mean.variance.df, start=list(theta=100))
  # coef(nls.modelNB)
  phi <- 1/coef(nls.modelNB) # to be used as fixed parameter later on
  #compute model fit
  # log plot, so exclude 0 from the x values
  x <- c(seq(0.1,1,0.1), seq(2,5000,1) )
  yQP <- x*coef(lm.D)
  yNB <- x + 1/coef(nls.modelNB) * x^2

  ###################################################
  ### code chunk number 5: Eventlevel3
  ###################################################
  plot(event.mean.variance.df$Mean, event.mean.variance.df$Variance, 
     xlab="Mean Event count", 
     ylab="Variance Event count",
     log="xy", las=1)
  abline(a=0, b=1, col=2, lwd=2)
  lines(x,yQP,col=3, lwd=2)
  lines(x,yNB,col=6, lwd=2)
  legend("topleft", c("Poisson","Quasi-Poisson", "Negative Binomial"), 
   text.col=c(2,3,6), box.lty=0);

  ###################################################
  ### code chunk number 6: pALLGlobalPhi.glm.nb
  ###################################################
  pALLGlobalPhi.glm.nb=data.frame(t(rep(NA,44)))
  for (i in 1:length(allEventtables)) {
    pALLGlobalPhi.glm.nb[i, ] = try(.fitNBglmModelsDSSPhi(allEventtables[[i]],dispersion(dispData)[i],dispersion(dispDataMeanCond)[i], phi, nbAll) ,silent=T)
  }

  ###################################################
  ### code chunk number 7: excl_errors
  ###################################################
  nonsing.events = which(!grepl("Error",pALLGlobalPhi.glm.nb[ ,1]))
  pALLGlobalPhi.glm.nb.nonsing=data.frame(t(rep(NA,44)))
  j=1
  for (i in nonsing.events) {
    pALLGlobalPhi.glm.nb.nonsing[j, ] = as.numeric(pALLGlobalPhi.glm.nb[i, ])
    j=j+1
  }
  pALLGlobalPhi.glm.nb=pALLGlobalPhi.glm.nb.nonsing
  colnames(pALLGlobalPhi.glm.nb) <- c("(0)A vs S","(0)I vs A",
                                    "(gb)A vs S","(gb)I vs A",
                                    "A vs S","I vs A",
                                    "(c)A vs S","(c)I vs A",
                                    
                                    "(0)bicS","(0)bicA","(0)bicI",
                                    "(gb)bicS","(gb)bicA","(gb)bicI",
                                    "bicS","bicA","bicI",
                                    "(c)bicS","(c)bicA","(c)bicI",
                                    
                                    "(0)codeS","(0)codeA","(0)codeI",
                                    "(gb)codeS","(gb)codeA","(gb)codeI", 
                                    "codeS","codeA","codeI",
                                    "(c)codeS","(c)codeA","(c)codeI",
                                    
                                    "(0)shS","(0)shA","(0)shI",
                                    "(gb)shS","(gb)shA","(gb)shI", 
                                    "shS","shA","shI",
                                    "(c)shS","(c)shA","(c)shI")
  rownames(pALLGlobalPhi.glm.nb) <- dataPart2[nonsing.events,1]
  pALLGlobalPhi.glm.nb = pALLGlobalPhi.glm.nb[!is.na(pALLGlobalPhi.glm.nb[ ,1]), ]

  ###################################################
  ### code chunk number 10: best model
  ###################################################
  bestmodel.table.n = apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)
  bestmodel.table = bestmodel.table.n
  bestmodel.table[bestmodel.table==1]="Poisson"
  bestmodel.table[bestmodel.table==2]="NB, global phi"
  bestmodel.table[bestmodel.table==3]="NB, DSS phi"
  bestmodel.table[bestmodel.table==4]="NB, cond DSS phi"
  bestmodel.singhes = c()
  for (i in 1:length(bestmodel.table.n)) {
    bestmodel.singhes[i] = c(pALLGlobalPhi.glm.nb[i,c(35,38,41,44)])[bestmodel.table.n[i]]
  }
  bestmodel.singhes = unlist(bestmodel.singhes)
  bestmodel = table(bestmodel.table)
  bestmodel2 = table(bestmodel.table,bestmodel.singhes)
  colnames(bestmodel2) = c("F","T")
  bestmodel3 = table(bestmodel.table,apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum),dnn=c("best model","singular hessian"))
  colnames(bestmodel3) = paste(colnames(bestmodel3), "Models")

  ###################################################
  ### code chunk number 11: glmnet
  ###################################################

  pALLGlobalPhi.glm.nb.glmnet = pALLGlobalPhi.glm.nb
  pALLGlobalPhi.glm.nb.glmnet$glmnet.pval = 1
  pALLGlobalPhi.glm.nb.glmnet$glmnet.code = 0
  singhes0 = which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min) == 1)# Event pour lesquels le modèle poissonien est plus adapté
  for (i in singhes0) {
    Xinter   = model.matrix(~cond*path,data= allEventtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]); 
    outinter = glmnet(Xinter,allEventtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]$counts,family="poisson",lambda=1e-4,alpha=0)
	  Xprinc   = model.matrix(~path+cond,data= allEventtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]); 
	  outprinc = glmnet(Xprinc,allEventtables[[which(rownames(dataPart2)==rownames(pALLGlobalPhi.glm.nb)[i])]]$counts,family="poisson",lambda=1e-4,alpha=0)
	  Pv       = 1-pchisq(deviance(outprinc) - deviance(outinter),df=1)
	  pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[i] = Pv
	  pALLGlobalPhi.glm.nb.glmnet$glmnet.code[i] = outinter$jerr
  }

  ###################################################
  ###################################################

  pALLGlobalPhi.glm.nb$final.pval.a.ia = 1
  i <- 1 #Poisson model
	li.singhes <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li.singhes] = pALLGlobalPhi.glm.nb.glmnet$glmnet.pval[li.singhes]
	
  i  <- 2 # negative binomial model, with global phi
	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,4]
	
  i  <- 3 # negative binomial model, phi estimated with DSS
	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,6]

  i  <- 4 # negative binomial model, phi estimated with DSS, conditionally to the expression mean
	li <- which(apply(pALLGlobalPhi.glm.nb[ ,c(11,14,17,20)],1,which.min)==i & apply(pALLGlobalPhi.glm.nb[ ,c(35,38,41,44)],1,sum)==0)
	pALLGlobalPhi.glm.nb$final.pval.a.ia[li] <- pALLGlobalPhi.glm.nb[li,8]
	
  pALLGlobalPhi.glm.nb$final.padj.a.ia <- p.adjust(pALLGlobalPhi.glm.nb$final.pval.a.ia, method="fdr")

  signifEvents <- cbind(dataPart2[nonsing.events, ],pALLGlobalPhi.glm.nb$final.padj.a.ia )[ pALLGlobalPhi.glm.nb$final.padj.a.ia <= 0.05, ]
  # sorting by deltaPSI / deltaF
  finalDelta <- NA
  #1st condition
  PSI.replicat1 <- c()
  for (j1 in 1:nr[1]) {
    nameUp <- paste('UP_',sortedconditions[j1], '_r', j1, '_Norm', sep='')
    nameLow <- paste('LP_', sortedconditions[j1],'_r', j1, '_Norm', sep='') 
    tmp.psi <- signifEvents[, nameUp] / (signifEvents[,nameUp] + signifEvents[ ,nameLow])
    PSI.replicat1 <- cbind( PSI.replicat1, tmp.psi)
    PSIcondI <- apply(PSI.replicat1, MARGIN=1, mean, na.rm = T)
    for (k1 in (2:n)) { # = 2nd condition
      PSI.replicat2 <- c()
      for (l in (1:nr[n])) {
        nameUp <- paste('UP_',sortedconditions[cumsum(nr)[n-1]+l], '_r', l, '_Norm', sep='')
        nameLow <- paste('LP_', sortedconditions[cumsum(nr)[n-1]+l],'_r', l, '_Norm', sep='') 
        tmp.psi <- signifEvents[, nameUp] / (signifEvents[,nameUp] + signifEvents[ ,nameLow])
        PSI.replicat2 <- cbind( PSI.replicat2, tmp.psi)
      }
      PSIcondJ <- apply(PSI.replicat2, MARGIN=1, mean, na.rm = T)
      #deltaPSI for cond i,j
      deltaPSIij <- abs( PSIcondJ- PSIcondI ) 
      finalDelta <- apply( cbind( finalDelta, deltaPSIij) , MARGIN = 1, max, na.rm = T)
    }
  }
  for (i in 2:(n-1)) {
    #1st condition
    PSI.replicat1 <- c()
    for (k2 in (1:nr[n-1])) {
      nameUp <- paste('UP_',sortedconditions[cumsum(nr)[n-1]+k2], '_r', k2, '_Norm', sep='')
      nameLow <- paste('LP_', sortedconditions[cumsum(nr)[n-1]+k2],'_r', k2, '_Norm', sep='') 
      tmp.psi <- signifEvents[, nameUp] / (signifEvents[,nameUp] + signifEvents[ ,nameLow])
      PSI.replicat1 <- cbind( PSI.replicat1, tmp.psi)
    }
    PSIcondI <- apply(PSI.replicat1, MARGIN=1, mean, na.rm = T)
    for (k3 in (i+1):n) {
      PSI.replicat2 <- c()
      for (l in 1:nr[n]) {
        nameUp <- paste('UP_',sortedconditions[cumsum(nr)[n-1]+l], '_r', l, '_Norm', sep='')
        nameLow <- paste('LP_', sortedconditions[cumsum(nr)[n-1]+l],'_r', l, '_Norm', sep='') 
        tmp.psi <- signifEvents[, nameUp] / (signifEvents[,nameUp] + signifEvents[ ,nameLow])
        PSI.replicat2 <- cbind( PSI.replicat2, tmp.psi)
      }
      PSIcondJ <- apply(PSI.replicat2, MARGIN=1, mean, na.rm = T)
      # deltaPSI for cond i,j
      deltaPSIij <- abs( PSIcondJ- PSIcondI ) 
      finalDelta <- apply( cbind( finalDelta, deltaPSIij) , MARGIN = 1, max, na.rm = T)
    }
  }

  signifEvents <- cbind(signifEvents, finalDelta)# adding DeltaPsi/f to the final table
  colnames(signifEvents)[length(colnames(signifEvents))] <- 'Deltaf/DeltaPSI'# renaming last columns
  signifEvents.sorted <- signifEvents[ order( finalDelta, decreasing = T), ]
  return(signifEvents.sorted)
}