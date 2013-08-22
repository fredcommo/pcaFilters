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

simulEset <- function(eset){
  M <- apply(eset, 1, mean, na.rm = TRUE)
  S <- apply(eset, 1, sd, na.rm = TRUE)
  Random <- generateRandom(eset, nrow(eset))
  SignifGrps <- generateGrps(M, S, ncol(eset), nrow(eset)*.01,
                                 nGrp = sample(2:4, 1), minP = 0.5, maxP = 0.9)
  signifGrps <- SignifGrps$Data
  Grps <- SignifGrps$grps
  colnames(Random) <- colnames(signifGrps) <- colnames(eset)
  return(rbind(Random, signifGrps))
}

multiTrace <- function(eset, pcaProbes){
  score <- lapply(1:3, function(i){
    tmp <- lapply((i+1):4, function(j){
      cat('From', i, 'to', j, '\n')
      model <- pcaTraceQ(eset, pcaProbes, Dim = i:j, Normalize = FALSE, Plot = FALSE)
      cbind(x = model$lModel$xfit, y = model$lModel$yfit)
    }
    )
    names(tmp) <- seq(i+1, i+length(tmp))
    tmp
  }
  )
  plot(c(-4,3),c(0, 30e3), type = 'n', xlab = 'Log10(quantiles)', ylab = 'Information')
  k = 1
  l <- list()
  for(i in 1:length(score)){
    tmp <- score[[i]]
    for(j in 1:length(tmp)){
      lines(tmp[[j]][,1], tmp[[j]][,2], col = k, lwd = 3)
      k = k + 1
      l[[k]] <- paste('from', i, 'to', names(tmp)[[j]])
    }
  }
  legend('bottomright', title = 'PCs', legend = l[-1], lwd = 3, col = 1:k)
}
  
pcaPerf <- function(eset, pcaProbes, Condition, trueList, Start = 1, End = 5,
                     threshold = 1e-2, mcCores = detectCores()/2){
  # From PC1 to PC4, adding up to PC5
  V <- apply(eset, 1, sd, na.rm = TRUE)
  if(!is.factor(trueList)) trueList <- factor(trueList)
  perfTable <- c()
  for (i in Start:(End-1)){
    cat('Testing from PC', i, '\n')
    tmpTable <- lapply(seq((i+1), End), function(j){
                  cat('Testin from', i, 'to', j, '\n')
                  tmpScore <- pcaTraceQ(eset, pcaProbes, Dim = i:j, Plot = FALSE)
                  # cat('Running multiple testing...\t')
                  tmpPerf <- mclapply(c(0.05, seq(0.1, 0.9, by = 0.1)), function(p){
                                select <- pcaSelectQ(tmpScore, p)
                                if(length(select)>1){
                                  Removed <- table(trueList[-select])/table(trueList)
                                  randRemoved = Removed[names(Removed) == 'random']
                                  signifRemoved = Removed[names(Removed) == 'signif']                                  #bestRemove <- .dimTest(eset, select, Condition, threshold, V, Remove = TRUE)
                                  #bestReplace <- .dimTest(eset, select, Condition, threshold, V, Remove = FALSE)
                                  #tabRemove <- table(trueList[select][bestRemove])/table(trueList)
                                  #tabReplace <- table(trueList[bestReplace])/table(trueList)
                                  return(cbind(first = i, last = j, p = p,
                                               nSelect = length(select),
                                               nRemoved = nrow(eset) - length(select),
                                               randRemoved = randRemoved,
                                               signifRemoved = signifRemoved,
                                               score = randRemoved*(1-signifRemoved)
                                               ))
                                        #   signatureRemove = length(bestRemove),
                                        #   randInTest = tabRemove[names(tabRemove) == 'random'],
                                        #   signifInTest = tabRemove[names(tabRemove) == 'signif']))
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
      #bestRemove <- .dimTest(eset, select, Condition, threshold, S, Remove = TRUE)
      #bestReplace <- .dimTest(eset, select, Condition, threshold, V, Remove = FALSE)
      #tabRemove <- table(trueList[select][bestRemove])/table(trueList)
      #tabReplace <- table(trueList[bestReplace])/table(trueList)
      return(cbind(prob = p, q = q,
                   nSelect = length(select),
                   nRemoved = nrow(eset) - length(select),
                   randRemoved = Removed[names(Removed) == 'random'],
                   signifRemoved = Removed[names(Removed) == 'signif']))
              #     signature = length(bestRemove),
              #     randInTest = tabRemove[names(tabRemove) == 'random'],
              #     signifInTest = tabRemove[names(tabRemove) == 'signif']
             #      signatureReplace = length(bestReplace),
             #      randInReplaceTest = tabReplace[names(tabReplace) == 'random'],
             #      signifInReplaceTest = tabReplace[names(tabReplace) == 'pseudo'])
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

compareScores <- function(eset, pcaScore, Grp, P = seq(.05, .9, by = .05), 
                          type = c('useDist', 'useChi2'), mcCores = detectCores()/2){
  require(multtest)
  # Compare perfs PCA-filter Vs. Var-filter by selecting the same number of features.
  type <- match.arg(type)
  cat(type, '\n')
  switch(type,
          useDist = {Select <- pcaSelectD},
          useChi2 = {Select <- pcaSelectQ})
  S <- apply(eset, 1, sd, na.rm = TRUE)
  # ROC curve PCA-filter Vs. Var-filter
  output <- lapply(P, function(p){
    #select <- pcaSelectQ(pcaScore, p)
    select <- Select(pcaScore, p)
#    if(length(select) == 0) select <- seq(1, nrow(eset))
    q <- quantile(S, probs = 1-length(select)/nrow(eset))
    idxVar <- which(S>q)
    cat('p:', p, '\t') #PCA select:', length(select), '\tVar select:', length(idxVar), '\n')
  #  mtPCA <- mt.maxT(eset[select,], classlabel=Grp, B = 5000)
  #  mtVar <- mt.maxT(eset[idxVar,], classlabel=Grp, B = 5000)
    
    remPCA <- nrow(eset) - length(select)
    if(length(select) > 0){
      randRemPCA <- sum(grepl('random', rownames(eset)[-select]))
      signifRemPCA <- sum(grepl('signif', rownames(eset)[-select]))
  #  signifPCA <- sum(mtPCA$adjp < 1e-3 & grepl('signif', rownames(mtPCA)))
    }
    else{
      randRemPCA <- sum(grepl('random', rownames(eset)))
      signifRemPCA <- sum(grepl('signif', rownames(eset)))      
    }
    
    remVar <- nrow(eset) - length(idxVar)
    if(length(select) > 0){
      randRemVar <- sum(grepl('random', rownames(eset)[-idxVar]))
      signifRemVar <- sum(grepl('signif', rownames(eset)[-idxVar]))
   # signifVar <- sum(mtVar$adjp < 1e-3 & grepl('signif', rownames(mtVar)))
    }
    else{
      randRemVar <- sum(grepl('random', rownames(eset)))
      signifRemVar <- sum(grepl('signif', rownames(eset)))
    }
    cbind(prop = p,
          remPCA = remPCA, randRemPCA = randRemPCA, signifRemPCA = signifRemPCA, #signifPCA = signifPCA,
          remVar = remVar, randRemVar = randRemVar, signifRemVar = signifRemVar)#, signifVar = signifVar)
  })#, mc.cores = mcCores)
  return(as.data.frame(do.call(rbind, output)))
}

plotRemoved <- function(eset, score = NULL, p, removed = NULL, Grps){
  myPalette <- c('darkblue', 'orangered', 'seagreen3', 'purple', 'steelblue1')
  par(mfrow = c(2, 2), cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, 
      las = 1, mar = c(5, 6.5, 4, 2))
  if(!is.null(score)){
    x <- score$lModel$xfit
    y <- score$lModel$yfit
    target <- .invRichardQ(score$lModel, p)
    xTarg <- target$xTarg
    yTarg <- target$yTarg
    plot(x, y, type = 'l', lwd = 3, col = 'blue',
         xlab = 'Log10(quantiles)', ylab = 'Information', main = 'Information curve')
    segments(x0 = c(-50, xTarg), y0 = c(yTarg, -50),
             x1 = c(xTarg, xTarg), y1 = c(yTarg, yTarg), lwd = 3, col = 'red')
    select <- pcaSelectQ(score, p)
    removed <- seq(1, nrow(eset))[-select]
    text(min(x)+1.5, yTarg, labels = paste0(p*100, '%'), cex = 1.5, font = 2, pos = 3)
  }
  M <- apply(eset, 1, mean, na.rm = TRUE)
  S <- apply(eset, 1, sd, na.rm = TRUE)
  smoothScatter(M, log10(S), colramp = colorRampPalette(c("white", "blue4")),
                xlab = 'Means', ylab = 'Log10(sDev)', main = 'Probes location')
  points(M[removed], log10(S[removed]), pch = 1, cex = 0.5, col = rgb(.7,.7,.7,.1))
  #  points(srM, log10(srS), pch = 19, cex = 0.75, col = 'cyan')
  points(M[-removed], log10(S[-removed]), pch = 19, cex = 0.75, col = rgb(1,0,0.2,.1))
  #   legend('bottomright', legend = c('original', 'random added', 'random grps added', 'main grps added'),
  #          pch = 19, col = c('purple', 'grey70', 'cyan', 'red3'), cex = 1.5, ncol = 2, bty = 'n')
  l1 <- paste0('Rejected: n = ', length(removed))
  l2 <- paste0('Selected: n = ', nrow(eset) - length(removed))
  legend('bottomright', legend = c(l1, l2),
         pch = 19, col = c('grey70', 'red3'), cex = 1.5, bty = 'n')
  pcaS <- prcomp(t(eset[-removed,]))
  pcaR <- prcomp(t(eset[removed,]))
  plotPCA(pcaS, col = myPalette[Grps], main = 'PCA on selected probes', pch = 19,
          xlim = range(pcaS$x[,1], pcaR$x[,1]),  ylim = range(pcaS$x[,2], pcaR$x[,2]))
  plotPCA(pcaR, col = myPalette[Grps], main = 'PCA on rejected probes', pch = 19,
          xlim = range(pcaS$x[,1], pcaR$x[,1]),  ylim = range(pcaS$x[,2], pcaR$x[,2]))
  par(op)
}

plotCompare <- function(compareScore, nRandom, nSignif,
                        type = c('random', 'signif', 'score')){
  type <- match.arg(type)
  switch(type,
         random = {y1 <- compareScore$randRemPCA/nRandom;
                   y2 <- compareScore$randRemVar/nRandom;
                   title <- 'Random probes removed';
                   legPos <- 'bottomright'},
         signif = {y1 <- compareScore$signifRemPCA/nSignif;
                   y2 <- compareScore$signifRemVar/nSignif;
                   title <- 'Significant probes removed';
                   legPos <- 'topleft'},
         score = {y1 <- compareScore$randRemPCA/nRandom*(1-compareScore$signifRemPCA/nSignif);
                  y2 <- compareScore$randRemVar/nRandom*(1-compareScore$signifRemVar/nSignif);
                  title <- 'Combined score';
                  legPos <- 'topright'})

  P <- compareScore$prop
  plot(P, y1, type = 'l', lwd = 5, col = 'blue3',
       xlab = 'Quantiles', ylab = 'Proportion', ylim = range(0,1), main = title)
  points(P, y1, lwd = 5, cex = 1.5, col = 'steelblue3')
  lines(P, y2, type = 'l', lwd = 5, col = 'red3')
  points(P, y2, lwd = 5, cex = 1.5, col = 'indianred3')
  legend(legPos, title = 'Filtered', legend = c('by PCA', 'by variance'), lwd = 3,
         col = c('blue3', 'red3'), cex = 2, bty = 'n')
} 

plotROC <- function(compareScores, nRandom, nSignif){
  sensitPCA <- c(compareScore$randRemPCA/nRandom, 1)
  specPCA <- c(compareScore$signifRemPCA/nSignif, 1)
  aucPCA <- AUC(specPCA, sensitPCA)
  sensitVar <- c(compareScore$randRemVar/nRandom, 1)
  specVar <- c(compareScore$signifRemVar/nSignif, 1)
  aucVar <- AUC(specVar, sensitVar)
  p <- wilcox.test(sensitPCA, sensitVar)$p.value

  plot(specPCA, sensitPCA, col = 'blue3', lwd = 5, type = 'l', ylim = range(0, 1),
       xlab = 'Significant removed', ylab = 'Random removed', main = 'Removed probes')
  lines(specVar, sensitVar, col = 'red3', lwd = 5)
  abline(0, 1, lwd = 5, lty = 3, col = 'grey30')
  l1 <- paste0('by PCA (AUC = ', round(aucPCA, 3), ')')
  l2 <- paste0('by variance (AUC = ', round(aucVar, 3), ')')
  l3 <- paste0('Wilcoxon: p = ', round(p, 3))
  legend('bottomright', title = 'Filtered', legend = c(l1, l2, l3), lwd = 3,
         col = c('blue3', 'red3', NA), cex = 2, bty = 'n')
}

AUC <- function(x, y){
  # Using trapezoide estimation
  auc <- lapply(2:length(x), function(i){
    dx <- abs(x[i]-x[i-1])
    dy <- abs(y[i]-y[i-1])
    y[i-1]*dx +1/2*dy*dx
  })
  return(do.call(sum, auc))
}

plotScore <- function(perf, type = c('random', 'signif', 'score')){
  type <- match.arg(type)
  switch(type,
         random = {Y <- perf$randRemoved; title = 'Random probes removed'},
         signif = {Y <- perf$signifRemoved; title = 'Significant probes removed'},
         score = {Y <- perf$score; title = 'Performance score'})
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
  score <- pcaTraceQ(eset, pcaProbes, Dim = dim1:dim2, Plot = FALSE)
  select <- pcaSelectQ(score, p)
  cols[-select] <- rgb(0,0,.8,.5)
  cexs[-select] <- 1.25
  pairs(pcaProbes$x[,1:5], col = cols, cex = cexs,
            main =paste('Removals using PC:', dim1, 'to', dim2))
}
  
visualizeSamples <- function(eset, pcaScore, filterValues = c(.025, .05, .1, .2, .5),...){
  par(mfrow = c(length(filterValues), 2), las = 1, mar = c(5, 6.5, 4, 2),
      cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
  #redDot = 10^(1/score$lModel$x.intercept)
  lapply(filterValues, function(p){
    cat('p:', p, '\t')
    select <- pcaSelectQ(pcaScore, p)
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
    points(M[-select], log10(S[-select]), pch = 19, cex = 0.2, col = 'grey90')
    legend('topright', legend = 'Removed probes', pch = 1, cex = 1.25,
           col = 'darkblue', border = 'darkblue', bty = 'n')
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
    pcaSelectQ(pcaScore, p)
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


radPairs <- function(data, p, idx, Cex = 1, Pcol = 'grey', Lcol = 'red3',...){
  
  # data: a data set (or matrix) of observations by row, and variates by columns.
  # smooth.fun: what smoothing function to use. "LM" performs and add a standard linear regression (default), "RLM" performs a robust linear regression, "SMOOTH" use a kernel smoothing.
  # cor.method: what method to compute the correlation coefficients. Default use 'Spearman'
  # Pcol: the points color.
  # Lcol: the line color.
  # CexCor: a parameter to adjust the size of the text in the upper panel.
  # Title: if TRUE (default) a subtitle is added indicating the smoothing and correlation methods used.
  #...: other optional graphic parameters.
  
  radius <- function(x, y, p, idx, Pcol = Pcol, Lcol = Lcol,...){
    D2 <- x[idx]^2+y[idx]^2
    q <- quantile(D2, p)
    theta <- seq(0, 90, by = 1)
    x1 <- sqrt(q)*cos(theta)
    y1 <- sqrt(q)*sin(theta)
    points(x, y, col = Pcol,...)
    lines(x1, y1, col = Lcol,...)
  }
  
  panel.cor <- function(x, y, p, idx, digits = 3, Cex = Cex){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    D2 <- x^2+y^2
    q <- quantile(D2[idx], p)
    txt1 <- format(c(q, 0.123456789), digits=digits)[1]
    txt1 <- paste('q95(D2):', txt1)
    txt2 <- paste('n=', sum(D2[-idx]<=q), sep="")
    txt3 <- paste('ratio=', round(sum(D2[idx]<=q)/sum(D2[-idx]<=q), 3), sep="")
    #text(0.5, 0.5, paste(txt1, '\n', txt2), cex = Cex)
    legend('center', legend = c(txt1, txt2, txt3), cex = Cex,
           lwd = c(2, -1, -1), col = c(Lcol, NA, NA), bty = 'n')
  }
  
  pairs(data,
        panel=function(x,y)radius(x, y, p, idx, Pcol, Lcol,...),
        upper.panel=function(x, y)panel.cor(x, y, p, idx, Cex = Cex)
        #upper.panel = NULL
  )
  #if(subTitle) title(sub = paste("Lower panel : Smoothing by", smooth.type, " / Upper panel : Correlation by", cor.type), ...)
}

roundness <- function(Data, Theta = seq(0, 2*pi, len = 500), minN = 50, Plot = FALSE){
  f <- function(x, slope){return(slope*x)}
  between <- function(y, y1, y2){which(y>=y1 & y<=y2)}
  x <- as.numeric(Data[,1])
  y <- as.numeric(Data[,2])
  
  Rtheta <- lapply(Theta, function(theta){
    eps = 0.0
    idx = NA
    while(length(idx) < minN){
      eps = eps + .025
      y1 <- f(x, sin(theta-eps)/cos(theta-eps))
      y2 <- f(x, sin(theta+eps)/cos(theta+eps))
{
        if(theta>pi/2 & theta<3*pi/2) idx <- between(y, y2, y1)
        else idx <- between(y, y1, y2)
      }
    }
    D <- sqrt(rowSums(Data[idx,]^2))
    return(max(D, na.rm = TRUE))
  }
  )
  
  Rtheta <- do.call(c, Rtheta)
  N <- sum(!is.na(Rtheta))
  Rhat <- 1/N*sum(Rtheta, na.rm = TRUE)
  ahat <- 2/N*sum(Rtheta*cos(Theta), na.rm = TRUE)
  bhat <- 2/N*sum(Rtheta*sin(Theta), na.rm = TRUE)
  Delta <- Rtheta - (Rhat + ahat*cos(Theta) + bhat*sin(Theta))
  if(Plot){
    plot(Data, cex = .5, asp = 1)
    points(Rtheta*cos(Theta), Rtheta*sin(Theta), col = 'blue')
    lines(Rhat*cos(Theta), Rhat*sin(Theta), col = 'red', lwd = 2)
  }
  return(Delta)
}

Classifier <- function(eset, select, resp, Family, alpha){
  .getFit <- function(train, test, resp, Family, alpha){
    model <- cv.glmnet(x = t(train), y = resp, family = Family, alpha = alpha, maxit = 1e4)
    return(predict(model, newx = t(test), s="lambda.min", type = 'class'))
  }
  .getSpe <- function(fit, obs){
    #  fit <- ifelse(1/(1+exp(-fit))<.5, levels(obs)[1], levels(obs)[2])
    res <- table(fit, obs)
    return(sum(diag(res))/sum(res))    
  }
  
  samp <- sample(1:ncol(eset), ncol(eset)*.6)
  
  fit1 <- try(.getFit(eset[select, samp], eset[select, -samp], resp[samp], Family, alpha), silent = TRUE)
  if(class(fit1)!="try-error") c1 <- .getSpe(fit1, resp[-samp])
  else c1 <- NA
  
  fit2 <- try(.getFit(eset[-select, samp], eset[-select, -samp], resp[samp], Family, alpha), silent = TRUE)
  if(class(fit2)!="try-error") c2 <- .getSpe(fit2, resp[-samp])
  else c2 <- NA
  
  return(c(c1 = c1, c2 = c2))
}

plotErr <- function(err, Fact,...){
  .pMean <- function(x){n <- length(x);return(prod(x)^(1/n))}
  .pSd <- function(x){n <- length(x);p <- .pMean(x); return(sqrt(p*(1-p)/n))}
  P <- unique(Fact)
  M <- as.vector(by(err, Fact, .pMean))
  S <- as.vector(by(err, Fact, .pSd))
  points(P, M,...); lines(P, M, ...)
  sapply(1:length(S), function(i,...){
    m <- M[i]; s <- S[i]; p <- P[i]
    segments(p, m-s, p, m+s)
    points(c(p, p), c(m-s, m+s), pch = '_')
#     segments(x0 = c(p, p-.02, p-.02), y0 = c(m-2*s, m-2*s, m+2*s),
#              x1 = c(p, p+.02, p+.02), y1 = c(m+2*s, m-2*s, m+2*s),...)
  })
}

saveFilters <- function(eset, score, P = c(.025, .05, .10, .20, .50)){
  S <- apply(eset, 1, sd, na.rm = TRUE)
  output <- lapply(P, function(p){
    pcaS <- pcaSelectQ(score, p)
    q <- quantile(S, 1-length(pcaS)/nrow(eset))
    varR <- which(S < q)
    list(pcaR = rownames(eset)[-pcaS], varR = rownames(eset)[varR])
  })
  names(output) <- paste0('F', P)
  output$var50 <- rownames(eset)[S<quantile(S, .5)]
  return(output)
}

