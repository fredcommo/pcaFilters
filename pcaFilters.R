####################################################
# Contains the 3 functions to use to perform a PCA-based filtering
# eset is a (p,n) matrix with p genes by rows, and n samples by columns
# pcaGenes <- prcomp(eset)
# score <- pcaTrace1.1(eset, pcaGenes)
# info <- pcaInfo(score)
# select <- pcaSelect(score, p = 0.05) # returns the indices corresponding to 5% of the information
####################################################
require(parallel)
require(multicore)

.computeInform <- function(D, a, Df){
  # Returns probes in the rejected region: D<=q(a)
  cat('\talpha:', a)
  alpha = 10^(-a)
  Pmax = 1 - alpha
  Q <- qchisq(p = Pmax, df = Df)
  return(which(D <= Q))
}

.computeTrace <- function(Data){
  # Returns the trace of the cov. matrix
  acpTest <- prcomp(Data)
  return(sum(acpTest$sdev^2))
}

pcaTrace1.1 <- function(Data, PCA, Dim = 2:3, weight = 0.25, Plot = TRUE,...){
  require(parallel)
  # Data: the original data set
  # PCA: the PCA = prcomp(Data), so a pca on probes.
  # Dim: the dimensions to consider on PCA$x
  # Plot the information curve, and returns a list:
    # m = nrow(Data)
    # n = ncol(Data)
    # PCdim = ncol(X)
    # Dist = D, the vector of distances
    # Score = Score, the Trace scores according to alpha, as a data set
    # lModel = model, the Richard model parameters.
  X = scale(PCA$x[,Dim])
  X <- as.data.frame(X)
  cat('Matrix size:', dim(PCA$x), '\tUsed dims:', Dim, '\n')
  D <- rowSums(X^2)
  #a.values <- 10^seq(-4, log10(15), len = 9)
  a.values <-  quantile(D, probs = seq(.01, .99, len = 10))
  
  cat('Testing quantiles... ')
  Score <- mclapply(a.values, function(a){
    inform <- which(D<=a) #.computeInform(D, a, ncol(X))
    nInform <- length(inform)
    if(nInform>=2){
      tmpTrace <- .computeTrace(t(Data[inform,]))
      tmpScore <- c(aValues = a, nProbes = nInform, Trace = tmpTrace)
      #cat('\t#Rejected:', length(inform), '\tTrace:', tmpTrace, '\n')
      tmpScore
      }
    }, mc.cores = 4)
  
  Score <- as.data.frame(do.call(rbind, Score))
  tmpTrace <- .computeTrace(t(Data))
  Score <- rbind(Score, cbind(aValues = max(D), nProbes = nrow(Data), Trace = tmpTrace))
  Score <- rbind(Score, cbind(aValues = max(D)*10, nProbes = nrow(Data), Trace = tmpTrace))
  rownames(Score) <- seq(1, nrow(Score))
  Score <- as.data.frame(Score)
  #cat('alpha:', max(D), '\t#rejected:', nrow(Data), '\tTrace:', tmpTrace, '\n')
  
  x <- as.numeric(log10(Score$aValues))
  y <- as.numeric(Score$Trace)
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]; x <- x[-na.index]
    }
  y = y/(max(y)*1.01)*100
  if(any(y <=0)) y[y<=0] <- 1e-3
  model <- .Richard.w5PL.v2(x, y, w = weight, Plot = Plot, add.points = TRUE,
                           xlab = expression(-Log10(quantile)), ylab = 'Information (%)',...)
  cat('Done.\n')
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}

pcaInfo <- function(pcaScore){
  # pcaScore: pcaTrace output
  # Returns the information table containing the number of informative probes given th proportion of information
  Fmax = pcaScore$lModel$top
  Fb = pcaScore$lModel$bottom
  xc = pcaScore$lModel$xc
  b = pcaScore$lModel$scal
  d = pcaScore$lModel$d
  informTable <- lapply(seq(0.05, 1, by = 0.05), function(p){
    yTarg = Fb + (Fmax - Fb)*p
    xTarg = xc - b*log10(((Fmax-Fb)/(yTarg - Fb))^(1/d) - 1)
    aTarg <- 10^(xTarg)
#     alpha = 10^(-aTarg)
#     Pmax = 1 - alpha
#     Q <- qchisq(p = Pmax, df = pcaScore$PCdim)  
#    inform <- which(pcaScore$Dist >= Q)
    inform <- which(pcaScore$Dist >= aTarg)
    Inform <- length(inform)
    NonInform <- pcaScore$m - Inform
    cbind(Prop = p, Inform = Inform, nonInform = NonInform,
          propInform = round(Inform/pcaScore$m, 4))
  })
  return(as.data.frame(do.call(rbind,informTable)))
}

pcaSelect <- function(pcaScore, p = 0.05){
  # pcaScore: pcaTrace output
  # p: the proprtion of information required
  # Returns the indices of the slected features, according to p.
  if(p==0) return(1:pcaScore$m)
  Fmax = pcaScore$lModel$top
  Fb = pcaScore$lModel$bottom
  xc = pcaScore$lModel$xc
  b = pcaScore$lModel$scal
  d = pcaScore$lModel$d
  
#  yTarg = (Fmax + Fb)*p
#  xTarg = xmid - b*log(((Fmax - Fb)/(yTarg - Fb))^(1/d) - 1)
  yTarg = Fb + (Fmax - Fb)*p
  xTarg = xc - b*log10(((Fmax-Fb)/(yTarg - Fb))^(1/d) - 1)
  
  aTarg <- 10^(xTarg)
#   alpha = 10^(-aTarg)
#   Pmax = 1 - alpha
#   Q <- qchisq(p = Pmax, df = pcaScore$PCdim)
#  inform <- which(pcaScore$Dist >= Q)
  inform <- which(pcaScore$Dist >= aTarg)
  cat('Removed:', pcaScore$m-length(inform),'\tSelected:', length(inform), '\n')
  return(inform)
}


.Richard.w5PL.v2 <- function(x, y, w = 0.25, Plot = FALSE, add.points = FALSE,
                             add.intercept = TRUE, pcol = "royalblue1", add.line = FALSE,
                             lcol = "navy", tan.col = "purple",...){ # Xlim = range(x), Ylim = range(y),
  
  # x    		: x-axis values
  # y				: y-axis values
  # w = 0.25			: weights coefficient
  # Plot = F			: if results have to be visualized
  # add.points = F		: add points on the plot
  # add.intercept = T	: add the intercept point on plot
  # pcol = "royalblue1"	: define points color	
  # add.line = F		: add the tangente line
  # lcol = "navy"		: define the regression color
  # tan.col = "purple"	: define the tangente line color
  # Xlim = range(x)		: the x-axis range
  # Ylim = range(y)		: the y-axis range
  # Title = ""		: to add a plot title
  
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]
    x <- x[-na.index]
  }
  
  # Fonction logistique 5PL
  Richard <- function(x, Fb, Fmax, b, c, d){
    y <- Fb + (Fmax - Fb)/(1 + 10^(-(x-c)/b))^d
    return(y)
  }
  
  # Fonction sce (somme carr? r?sidus) avec pond?rations
  sce.5P <- function(param, xobs, yobs, w) {
    Fb <- 0 #param[1]
    Fmax <- param[2]
    b <- param[3]
    c <- param[4]
    d <- param[5]
    ytheo <- Richard(xobs, Fb, Fmax, b, c, d)
    sq.res <- (yobs - ytheo)^2
    weights <- 1/sq.res^w
    return(sum(weights*sq.res))
  }
  
  # Fonction sce (somme carr? r?sidus) avec pond?rations
  sce.5P.diag <- function(yobs, ytheo, w) {
    sq.res <- (yobs - ytheo)^2
    weights <- 1/sq.res^w
    return(weights)
  }
  
  # initialisation des parametres
  Fb.ini = min(y)
  Fmax.ini = max(y)	#*1.05
  c.ini = (max(x) + min(x))/2
  z <- (y)/(Fmax.ini - y)
  if (any(abs(z)==Inf)) z[abs(z)==Inf] <- NA
  b.ini = coef(lm(x~log(z)))[2]
  # b.ini = 1					
  d.ini = 1
  init <- c(Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)
  
  # Estimation du modele
  best<-nlm(f = sce.5P, p = init, xobs = x, yobs = y, w = w)
  
  # R?cup?ration des param?tres
  Fb <- best$estimate[1]
  Fmax <- best$estimate[2]
  b <-best$estimate[3]
  xc <- best$estimate[4]
  d <- best$estimate[5]
  
  # Diagnostic de r?gression
  yfit <- Richard(x, Fb, Fmax, b, xc, d)
  weights <- sce.5P.diag(y, yfit, w)
  lm.test <- lm(yfit~y)	#, weights = weights)
  r.sq <- summary(lm.test)$adj.r.squared
  lm.slope <-coef(lm.test)[2]
  p.slope <- summary(lm.test)$coefficients[2,4]
  
  # Estimation des valeurs pour graphique
  newx <- seq(min(x), max(x), length=100)						
  newy <- Richard(newx, Fb, Fmax, b, xc, d)
  
  # coordonn?es du pt d'inflexion
  Xflex = xc + b*log10(d)
  Yflex = Fb + (Fmax-Fb)*(d/(1 + d))^d
  
  # coordonn?es du pt rep = 0.5
  Y50 = (max(newy) + min(newy))/2
  X50 = xc - b*log10(((Fmax-Fb)/(Y50 - Fb))^(1/d) - 1)
  
  # pente au pt d'inflexion
  B = (Fmax - Fb)/b*log(10)*(d/(d+1))^(d+1)
  
  A = Yflex  - B*(Xflex)
  
  # pente finale
  x.ini <- x[1]
  x.end <- x[length(x)]
  y.ini <- Richard(x.ini, Fb, Fmax, b, xc, d)
  y.end <- Richard(x.end, Fb, Fmax, b, xc, d)
  Bf = (d*(Fmax-Fb)/b)*10^(-1/b*(x.end-xc))*(1+10^(-1/b*(x.end-xc)))^(-d-1)
  Af = y.end - Bf*x.end
  
  # x.intercept
  Cy0 = -A/B
  
  # Repr?sentations graphiques
  if(Plot){
    plot(y~x, type = 'n',...)
    lines(newy~newx, col = pcol,...)
    lines(I(A+B*x)~x, col = tan.col, lwd = 2)
    abline(h = 0, lwd = 1, lty = 3, col = 'grey75')
    if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
    
    if(add.points) points(y~x, col = pcol,...)
    
    if(add.line){
      lines(newy~newx, col = lcol, lwd = 2)
      lines(I(A+B*x)~x, col = tan.col, lwd = 1)
      if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
    }
  }
  return(list(bottom = Fb, top = Fmax, xc = xc, scal = b, d = d,
              Xflex = Xflex, Yflex = Yflex, slope = B, x.intercept = Cy0, Yini = y.ini, Yend = y.end, end.slope = Bf,
              lm.rsq = r.sq, lm.slope = lm.slope, p.value = p.slope, xfit = newx, yfit = newy))
}

pcaPerf <- function(eset, pcaProbes, Condition, trueList, Start = 1, End = 5,
                     threshold = 1e-2, mcCores = detectCores()/2){
  # From PC1 to PC4, adding up to PC5
  V <- apply(eset, 1, sd, na.rm = TRUE)
  if(!is.factor(trueList)) trueList <- factor(trueList)
  perfTable <- c()
  for (i in Start:(End-1)){
    tmpTable <- lapply(seq((i+1), End), function(j){
                  cat('Testin from', i, 'to', j, '\n')
                  tmpScore <- pcaTrace1.1(eset, pcaProbes, Dim = i:j, Plot = FALSE)
                  cat('Running multiple testing...\t')
                  tmpPerf <- mclapply(c(0.05, seq(0.1, 0.9, by = 0.1)), function(p){
                                select <- pcaSelect(tmpScore, p)
                                if(length(select)>1){
                                  Removed <- table(trueList[-select])/table(trueList)
                                  bestRemove <- .dimTest(eset, select, Condition, threshold, V, Remove = TRUE)
                                  #bestReplace <- .dimTest(eset, select, Condition, threshold, V, Remove = FALSE)
                                  tabRemove <- table(trueList[select][bestRemove])/table(trueList)
                                  #tabReplace <- table(trueList[bestReplace])/table(trueList)
                                  return(cbind(first = i, last = j, p = p,
                                           nSelect = length(select),
                                           nRemoved = nrow(eset) - length(select),
                                           randRemoved = Removed[names(Removed) == 'random'],
                                           signifRemoved = Removed[names(Removed) == 'signif'],
                                           signatureRemove = length(bestRemove),
                                           randInTest = tabRemove[names(tabRemove) == 'random'],
                                           signifInTest = tabRemove[names(tabRemove) == 'signif']))
                                        #   signatureReplace = length(bestReplace),
                                        #   randInReplaceTest = tabReplace[names(tabReplace) == 'random'],
                                        #   signifInReplaceTest = tabReplace[names(tabReplace) == 'signif']))
                                  }
                                }, mc.cores = mcCores)
                  return(do.call(rbind, tmpPerf))
                  })
    cat('Done.\n')
    tmpTable <- as.data.frame(do.call(rbind, tmpTable))
    perfTable <- rbind(perfTable, tmpTable)
    }
  return(perfTable)
}

varPerf <- function(eset, probs , Condition, trueList,
                    threshold = 1e-3, mcCores = detectCores()/2){
  S <- apply(eset, 1, sd, na.rm = TRUE)
  tmpPerf <- mclapply(probs, function(p){
    q <- quantile(S, probs = p)
    select <- which(S > q)
    if(length(select)>1){
      Removed <- table(trueList[-select])/table(trueList)
      bestRemove <- .dimTest(eset, select, Condition, threshold, S, Remove = TRUE)
      #bestReplace <- .dimTest(eset, select, Condition, threshold, V, Remove = FALSE)
      tabRemove <- table(trueList[select][bestRemove])/table(trueList)
      #tabReplace <- table(trueList[bestReplace])/table(trueList)
      return(cbind(prob = p, q = q,
                   nSelect = length(select),
                   nRemoved = nrow(eset) - length(select),
                   randRemoved = Removed[names(Removed) == 'random'],
                   signifRemoved = Removed[names(Removed) == 'signif'],
                   signature = length(bestRemove),
                   randInTest = tabRemove[names(tabRemove) == 'random'],
                   signifInTest = tabRemove[names(tabRemove) == 'signif']
             #      signatureReplace = length(bestReplace),
             #      randInReplaceTest = tabReplace[names(tabReplace) == 'random'],
             #      signifInReplaceTest = tabReplace[names(tabReplace) == 'pseudo'])
             ))
    }
  }, mc.cores = mcCores)
  return(as.data.frame(do.call(rbind, tmpPerf)))
}

.dimTest <- function(Data, select, Condition, threshold, V, Remove){
  require(multtest)
  if(Remove)
    Data <- Data[select,]
  else{
    n <- nrow(Data) - length(select)
    Data[-select,] <- matrix(rnorm(n*ncol(Data), 0, median(V)), n, ncol(Data))
    }  
  cat('Running test on', nrow(Data),'features...\n')
  Test <- mt.maxT(Data, classlabel=Condition, B = 5000)
  Best <- Test$index[Test$adjp < threshold]
  cat('Signature:', length(Best), 'of', nrow(Data), 'at p<', threshold,'\n')
  return(Best)
}

compareScores <- function(eset, pcaScore, Grp, P = seq(.05, .9, by = .05)){
  require(multtest)
  # Compare perfs PCA-filter Vs. Var-filter by selecting the same number of features.
  
  S <- apply(eset, 1, sd, na.rm = TRUE)
  # ROC curve PCA-filter Vs. Var-filter
  output <- lapply(P, function(p){
    select <- pcaSelect(score, p)
    q <- quantile(S, probs = 1-length(select)/nrow(eset))
    idxVar <- which(S>=q)
    cat('p:', p, '\tPCA select:', length(select), '\tVar select:', length(idxVar), '\n')
    mtPCA <- mt.maxT(eset[select,], classlabel=Grp, B = 5000)
    mtVar <- mt.maxT(eset[idxVar,], classlabel=Grp, B = 5000)
    
    remPCA <- nrow(eset) - length(select)
    randRemPCA <- sum(grepl('random', rownames(eset)[-select]))
    signifRemPCA <- sum(grepl('signif', rownames(eset)[-select]))
    signifPCA <- sum(mtPCA$adjp < 1e-3 & grepl('signif', rownames(mtPCA)))
    
    remVar <- nrow(eset) - length(idxVar)
    randRemVar <- sum(grepl('random', rownames(eset)[-idxVar]))
    signifRemVar <- sum(grepl('signif', rownames(eset)[-idxVar]))
    signifVar <- sum(mtVar$adjp < 1e-3 & grepl('signif', rownames(mtVar)))
    
    cbind(prop = p,
          remPCA = remPCA, randRemPCA = randRemPCA, signifRemPCA = signifRemPCA, signifPCA = signifPCA,
          remVar = remVar, randRemVar = randRemVar, signifRemVar = signifRemVar, signifVar = signifVar)
  })
  return(as.data.frame(do.call(rbind, output)))
}

plotCompare <- function(compareScore, nRandom, nSignif,
                        type = c('random', 'signif', 'score')){
  type <- match.arg(type)
  switch(type,
         random = {y1 <- compareScore$randRemPCA/nRandom;
                   y2 <- compareScore$randRemVar/nRandom;
                   title <- 'Random probes removed'},
         signif = {y1 <- compareScore$signifPCA/nSignif;
                   y2 <- compareScore$signifVar/nSignif;
                   title <- 'Significant probes detected'},
         score = {y1 <- compareScore$randRemPCA/nRandom*compareScore$signifPCA/nSignif;
                  y2 <- compareScore$randRemVar/nRandom*compareScore$signifVar/nSignif;
                  title <- 'Combined score'})

  P <- compareScore$prop
  plot(P, y1, type = 'l', lwd = 5, col = 'blue3',
       xlab = 'Quantiles', ylab = 'Proportion', ylim = range(0,1), main = title)
  points(P, y1, lwd = 5, cex = 1.5, col = 'steelblue3')
  lines(P, y2, type = 'l', lwd = 5, col = 'red3')
  points(P, y2, lwd = 5, cex = 1.5, col = 'indianred3')
  legend('topleft', title = 'Filtered', legend = c('by PCA', 'by variance'), lwd = 3,
         col = c('blue3', 'red3'), cex = 2, bty = 'n')
} 

plotROC <- function(compareScores, nRandom, nSignif){
  sensitPCA <- compareScore$randRemPCA/nRandom
  specPCA <- compareScore$signifRemPCA/nSignif
  sensitVar <- compareScore$randRemVar/nRandom
  specVar <- compareScore$signifRemVar/nSignif
  plot(specPCA, sensitPCA, col = 'blue3', lwd = 5, type = 'l', ylim = range(0, 1),
       xlab = 'Significant removed', ylab = 'Random removed', main = 'Removed probes')
  lines(specVar, sensitVar, col = 'red3', lwd = 5)
  abline(0, 1, lwd = 5, lty = 3, col = 'grey30')
  legend('bottomright', title = 'Filtered', legend = c('by PCA', 'by variance'), lwd = 3,
         col = c('blue3', 'red3'), cex = 2, bty = 'n')
}
  
plotScore <- function(perf, type = c('random', 'signif', 'score')){
  type <- match.arg(type)
  switch(type,
         random = {Y <- perf$randRemoved; title = 'Random probes removed'},
         signif = {Y <- perf$signifRemoved; title = 'Significant probes removed'},
         score = {Y <- perf$randRemoved*(1-perf$signifRemoved); title = 'Performance score'})
  plot(c(0,1), c(0,1.2), type= 'n', xlab = 'Filter value', ylab = 'Score', main = title)
  leg <- cols <- ltys <- c()
  for(i in unique(perf$first)){
    k = 1
    for(j in (i+1):max(perf$last)){
      idx <- which(perf$first == i & perf$last== j)
      x <- perf$p[idx]
      y <- Y[idx]
      lines(x, y, lwd = 3, col = i+1, lty = k)
      leg <- c(leg, paste0('PC', i, ' to PC', j))
      cols <- c(cols, i+1)
      ltys <- c(ltys, k)
      k = k + 1
    }
  }
  legend('topright', legend = leg, col = cols, lty = ltys, lwd = 3, cex = 1.5, ncol = i-1, bty = 'n')
}

plotVarPerf <- function(varPerf){
  x <- varPerf$prob
  y1 <- varPerf$randRemoved #*(1 - varPerf$signifRemoved)
  y2 <- varPerf$signifRemoved #*(1 - varPerf$randInTest)
  plot(x, y1, ylim = range(0, 1), type = 'l', col = 'blue3', lwd = 3,
       xlab = 'Quantiles', ylab = 'Proportion', las = 1)
  lines(x, y2, col = 'red3', lwd = 3)
  lines(x, y1*(1-y2), col = 'green3', lwd = 3)
  legend('topleft', legend = c('Ran. Removed', 'Signif removed', 'Score'),
         lwd = 3, col = c('blue3', 'red3', 'green3'), cex = 1, bty = 'n')
}

visualizeProbes <- function(eset, pcaProbes, p = 0.1, dim1 = 2, dim2 = 3){
  # Visualize PCA selection depending on axes
  cols <- rep('grey75', nrow(eset))
  cols[grep('random', rownames(eset))] <- rgb(0,0.7,0.5,0.5) #green
  cols[grep('signif[^_rand]', rownames(eset))] <- 'red3'
  cols[grep('signif_rand', rownames(eset))] <- 'cyan'
  
#   cols <- ifelse(grepl('random', rownames(eset)), rgb(0,0.7,0.5,0.5),
#                 ifelse(grepl('signif', rownames(eset)), rgb(.8,0,0,.5), rgb(.8,.8,.8,.5)))
  cexs <- rep(1, nrow(eset))
  score <- pcaTrace1.1(eset, pcaProbes, Dim = dim1:dim2, Plot = FALSE)
  select <- pcaSelect(score, p)
  cols[-select] <- rgb(0,0,.8,.5)
  cexs[-select] <- 1.25
  pairs(pcaProbes$x[,1:3], col = cols, cex = cexs,
            main =paste('Removals using PC:', dim1, 'to', dim2))
}
  
visualizeSamples <- function(eset, pcaScore, filterValues = c(0.05, 0.1, 0.2),...){
  par(mfrow = c(length(filterValues), 2), mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
  #redDot = 10^(1/score$lModel$x.intercept)
  lapply(filterValues, function(p){
    cat('p:', p, '\t')
    select <- pcaSelect(pcaScore, p)
    n <- length(select)
    pcaS <- prcomp(t(eset[select,]))
    pcaR <- prcomp(t(eset[-select,]))
    plotPCA(pcaS, main = paste(n, 'informative probes'),...)
    plotPCA(pcaR, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
            main = paste(nrow(eset)-n,'rejected probes'),...)  
    })
  par(op)  
}

plotFilter <- function(eset, selectList,...){
  M <- apply(eset, 1, mean, na.rm = TRUE)
  S <- apply(eset, 1, sd, na.rm = TRUE)
  par(mfrow = c(length(selectList), 3), mar = c(5, 5, 4, 2.5), cex.main = 1.5, cex.lab = 1.25, cex.axis = 1)
  lapply(1:length(selectList), function(i){
    select <- selectList[[i]]
    smoothScatter(M, log10(S), colramp = colorRampPalette(c("white", "blue4")),
                main = 'Filtered probes')
    points(M[-select], log10(S[-select]), pch = 19, cex = 0.2, col = 'grey30')
    legend('bottomright', legend = 'Removed probes', pch = 19, cex = 1.25, col = 'grey30', bty = 'n')
    pcaS <- prcomp(t(eset[select,]))
    pcaR <- prcomp(t(eset[-select,]))
    plotPCA(pcaR, xlim = range(pcaS$x[,1], pcaR$x[,1]), ylim = range(pcaS$x[,2], pcaR$x[,2]),
          main = paste(nrow(eset)-length(select), 'Rejected probes'),...)
    plotPCA(pcaS, xlim = range(pcaS$x[,1], pcaR$x[,1]), ylim = range(pcaS$x[,2], pcaR$x[,2]),
            main = paste(length(select), 'Selected probes'),...)}
    )
  par(op)
}

pcaFilter <- function(eset, pcaScore, filterValues = 0.05,...){
  #redDot = 10^(1/score$lModel$x.intercept)
  select <- lapply(filterValues, function(p){
    cat('p:', p, '\t')
    pcaSelect(pcaScore, p)
  })
  plotFilter(eset, select,...)  
}

randomFilter <- function(eset, n,...){
  M <- apply(eset, 1, mean, na.rm = TRUE)
  S <- apply(eset, 1, sd, na.rm = TRUE)
  select <- sample(1:nrow(eset), n)
  plotFilter(eset, list(select),...)
}

varFilter <- function(eset, n,...){
  M <- apply(eset, 1, mean, na.rm = TRUE)
  S <- apply(eset, 1, sd, na.rm = TRUE)
  Q <- quantile(S, probs = 1-n/nrow(eset))
  highVarList <- lapply(Q, function(q){lowVar <- which(S > q)})
  plotFilter(eset, highVarList,...)
}

# Multiple testing
mtTest <- function(eset, grp, thresh = 1e-3, B = 5000){
  require(multtest)
  n <- nrow(eset)
  Test <- mt.maxT(eset, classlabel = grp, B = B)
  Test <- Test[order(Test$index),]
  best <- Test$index[Test$adjp<thresh]
  nRandom <- sum(grepl('random', rownames(eset)))
  foundRandom <- sum(grepl('random', rownames(eset)[best]))
  nSignif <- sum(grepl('signif', rownames(eset)))
  foundSignif <- sum(grepl('signif', rownames(eset)[best]))
  cat('\nThreshold adj. p-value:', thresh, '\n')
  cat('#Signature:', length(best), 'of', n, '\tprop:', round(length(best)/n, 4)*100, '\n')
  cat('#Random:', foundRandom, 'of', nRandom, '\tprop:', round(foundRandom/nRandom, 4)*100, '\n')
  cat('#Signif:', foundSignif, 'of', nSignif, '\tprop:', round(foundSignif/nSignif, 4)*100, '\n')
  return(Test)
}

.modelGlm <- function(x, y){
  x <- as.numeric(x)
  y <- c(0,1)[y]
  model <- glm(y ~ x, family = binomial)
  pvalue <- summary(model)$coefficients[2,4]
  return(pvalue)
}

glmProfile <- function(x, y, probeName){
  x <- as.numeric(x)
  y <- c(0,1)[y]
  pTtest <- t.test(x~y)$p.value
  model <- glm(y ~ x, family = binomial)
  pGlm <- summary(model)$coefficients[2,4]
  newx <- seq(min(x), max(x), len = 100)
  fit <- predict(model, newdata = list(x = newx), type = 'response')
  plot(x, y, cex = 1.25, col = Cols, xlim = range(quantile(x, probs = c(0.02, 1))),
       xlab = 'Log10(int)', ylab = 'Probability',
       main = paste0('probe ', probeName, ' (t.test p: ', signif(pTtest, 3),')'))
  lines(newx, fit, lwd = 3)
  legend('right', legend = paste('logist. reg.\np:',signif(pGlm,3)), cex = 1.25, bty = 'n')
  k = k+1
  axis(side = 4, at = c(0,1), labels = c('Adk', 'Sq'))
}

GMmodel <- function(x, G = 1:9, resamp = length(x),...){
  x <- x[sample(1:length(x), min(length(x), resamp))]
  d <- density(x)
  model <- Mclust(x, G = G)
  n = length(x); nG = model$G; m <- model$parameters$mean; p <- model$parameters$pro
  s <- sqrt(model$parameters$variance$sigmasq); if(length(s)==1) s <- rep(s, nG)
  plot(d, lwd = 2, ...)
  for (i in 1:nG){
    xtmp <- seq(min(d$x), max(d$x), len = 1000)
    dtmp <- dnorm(xtmp, m[i], s[i])
    polygon(xtmp, dtmp*p[i], col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
  }  
}

saveFilters <- function(eset, score, P = c(.05, .10, .20, .50)){
  S <- apply(eset, 1, sd, na.rm = TRUE)
  output <- lapply(P, function(p){
    pcaS <- pcaSelect(score, p)
    q <- quantile(S, 1-length(pcaS)/nrow(eset))
    varR <- which(S < q)
    list(pcaR = rownames(eset)[-pcaS], varR = rownames(eset)[varR])
  })
  names(output) <- paste0('F', P)
  return(output)
}

