require(parallel)
require(multicore)

.computeBounds <- function(M, S, By = 0.2){
  # M: vector of means
  # S: vector of variances (or sd)
  # By: size of intervals
  mCuts <- cut(M, breaks = seq(min(M, na.rm = TRUE), max(M, na.rm = TRUE), by = By)) #
  labs <- levels(mCuts)
  mBounds <- cbind.data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                              upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  mBounds$lower[1] <- min(M, na.rm = TRUE)
  mBounds$upper[nrow(mBounds)] <- max(M, na.rm = TRUE)  
  mBounds <- cbind.data.frame(lower = mBounds$lower, med = (mBounds$upper + mBounds$lower)/2, upper = mBounds$upper)
  sBounds <- lapply(seq(1, nrow(mBounds)), function(x){
    index <- as.numeric(which(M >= mBounds$lower[x] & M < mBounds$upper[x]))
    if(length(index)<1) cbind(S = NA, M = mBounds$med[x], Ml = mBounds$lower[x], Mu = mBounds$upper[x])
    else cbind(S = S[index], M = mBounds$med[x], Ml = mBounds$lower[x], Mu = mBounds$upper[x])
  })
  sBounds <- do.call(rbind, sBounds)
  rownames(sBounds) <- seq(1, nrow(sBounds))
  return(as.data.frame(sBounds))
}

.dfunc <- function(x){
  sbar <- sd(x, na.rm = TRUE)
  mbar <- mean(x, na.rm = TRUE)
  return(1/(sbar*sqrt(2*pi))*exp(-1/2*((x-mbar)/sbar)^2))
}

.pickM <- function(dM){
  m = as.numeric(sample(dM$x, 1, prob = dM$y))
  return(m)
}

.pickS <- function(MStable, m){
  j <- which(MStable$Ml <= m & MStable$Mu > m)
  #  if(length(j)<1) {s <- min(dS$x, na.rm = TRUE); cat('min S used for m:', m, '\n')}
  if(length(j)<1) {s <- min(MStable$S, na.rm = TRUE); cat('min S used for m:', m, '\n')}
  else {
    tmpS <- MStable$S[j]
    if(length(tmpS)>2)
      s <- sample(tmpS, 1, prob = pchisq(tmpS, df = 1, ncp = abs(m), lower.tail = FALSE))
    else s <- sample(tmpS, 1)
  }
  return(s)
}

generateRandom <- function(Data, n, mcCores = 1){
  require(parallel)
  require(multicore)
  cat('Resampling', n, 'rows...\t')
  # Returns a (n,p) matrix of premutation values from p randow Data rows
  Samp <- Data[sample(1:nrow(Data), n, replace = TRUE),]
  output <- mclapply(1:n, function(i){
    return(sample(as.numeric(Samp[i,])))
    }, mc.cores = mcCores)
  output <- do.call(rbind, output)
  rownames(output) <- paste0('random', seq(1, nrow(output)))
  colnames(output) <- paste0('sample', seq(1, ncol(output)))
  cat(nrow(output), 'done\n')
  return(output)
}
# 
# 
# .generateGrps <- function(Mvector, Svector, n, p, grps = NULL, nGrp = 2,
#                           minP = 0.3, maxP = 0.8, mcCores = 1){
#   require(parallel)
#   require(multicore)
#   # Provide vector of means and coresspondong Sdev
#   # n, p : number of samples (provides if grps if not provided), number of probes, resp.
#   # minP, maxP : min and max proportion of samples in a grp for which a specif probe is generated.
#   # Returns a (p, n) matrix of probes specific to grps and the vector of grps as a list.
#   
#   # if grps not provided, generate random grps according to nGrps
#   if(is.null(grps))
#     grps <- factor(rbinom(n, nGrp - 1, 0.5), labels = LETTERS[1:nGrp])
#   else if(!is.factor(grps)) grps <- as.factor(grps)
#   if(length(grps) != n) n <- length(grps)
#   nGrps <- nlevels(grps)
#   dM <- density(Mvector)
#   MStable <- .computeBounds(Mvector, Svector, By = 0.5) #Improve this func to avoid NAs in the S column.
#   if(any(is.na(MStable$S))) MStable <- MStable[!is.na(MStable$S),]
#   output <- mclapply(seq(1, p),
#                    function(x){ if(x%%ceiling(p/10) == 0) cat(x, '\t')
#                                 m <- max(min(MStable$M),.pickM(dM))
#                                 if(is.na(m)) stop('m is NA')
#                                 grp <- sample(levels(grps), 1)
#                                 s <- .pickS(MStable, m)
#                                 values <- rnorm(n, m, s)
#                                 # with respect to grps but some samples only, according to minP/maxP
#                                 idx <- sample(which(grps == grp))
#                                 rbi <- rbinom(length(idx), 1, prob = sample(seq(minP, maxP, by = 0.05), 1))
#                                 idx <- idx[rbi == 1]
#                                 change <- sample(c('inc', 'dec'), 1)
#                                 if(change == 'inc')
#                                   newmu <- min(max(MStable$M)*0.95, m*sample(seq(1.5, 2.5, by = 0.2), 1))
#                                 else
#                                   newmu <- max(min(MStable$M)*1.05, m/sample(seq(1.5, 2.5, by = 0.2), 1))
#                                 s <- .pickS(MStable, newmu)
#                                 values[idx] <- rnorm(length(idx), newmu, s)
#                                 if(any(is.na(values))) stop('s', s, 'is NA for m:', m)                                      
#                                 return(values)
#                    }, mc.cores = mcCores)
#   output <- do.call(rbind, output)
#   rownames(output) <- paste0('signif', seq(1, nrow(output)))
#   colnames(output) <- paste0('sample', seq(1, ncol(output)))
#   return(list(Data = output, Grps = grps))
# }

# .generateGrps2
# vector or means using .pickM(dM, n): output = n means
# vector on s using .pickS(MStable, m): output = s sdev with respect to m
# values matrix(p, n), each row is rnorm(n, mi, si)
# grps matrix(p, n), each row is 0/1 vector of group assignment
# for each p in values, assign new m/s where grp matrix = 1
# return values

generateGrps <- function(Mvector, Svector, n, p, grps = NULL, nGrp = 3,
                           minP = 0.5, maxP = 0.9, mcCores = detectCores()/2){
  p <- round(p)
  if(is.null(grps))
    grps <- factor(rbinom(n, nGrp - 1, c(.75, .25)), labels = LETTERS[1:nGrp])
  else if(!is.factor(grps)) grps <- as.factor(grps)
  if(length(grps) != n) n <- length(grps)
  cat('Generating', p, 'probes with respect to', levels(grps),'...')

  nGrps <- nlevels(grps)
  dM <- density(Mvector)
  MStable <- .computeBounds(Mvector, Svector, By = 0.5) #Improve this func to avoid NAs in the S column.
  if(any(is.na(MStable$S))) MStable <- MStable[!is.na(MStable$S),]

  cat('\nMean values...')
  meanValues <- mclapply(1:p, function(i){max(min(MStable$M),.pickM(dM))}, mc.cores = mcCores)
  meanValues <- do.call(c, meanValues)

  cat('\nsDev values...')
  sdValues <- mclapply(meanValues, function(mu){.pickS(MStable, mu)}, mc.cores = mcCores)
  sdValues <- do.call(c, sdValues)

  cat('\nProbes values...')
  values <- mclapply(1:p, function(i){rnorm(n, meanValues[i], sdValues[i])}, mc.cores = mcCores)
  values <- do.call(rbind, values)

  cat('\nAssignment matrix (Ig) according to groups...')
  Ig <- mclapply(1:p, function(i){
    grp <- sample(levels(grps), 1)
    rbi <- rbinom(sum(grps == grp), 1, prob = sample(seq(minP, maxP, by = 0.05), 1))
    tmp <- ifelse(grps != grp, 0, rbi)
    }, mc.cores = mcCores)
  Ig <- do.call(rbind, Ig)
    
  cat('\nIncrease/decrease factors...')
#  xgRange <- c(seq(-6, -3, by = .01), seq(1, 4, by = .01))
  xgRange <- seq(0, 4, by = .01)
  xg <- sample(xgRange, p, replace = TRUE)
  newMu <- meanValues + (1 + xg)
  if(any(newMu < min(MStable$M))){
    nLo <- sum(newMu < min(MStable$M))
    newMu[newMu < min(MStable$M)] <- quantile(MStable$M, .25)*(1 + rnorm(nLo, 0, .1))
    }
  if(any(newMu > max(MStable$M))){
    nHi <- sum(newMu > max(MStable$M))
    newMu[newMu > max(MStable$M)] <- quantile(MStable$M, .75)*(1 + rnorm(nHi, 0, .1))
    }
  newSd <- mclapply(newMu, function(mu){.pickS(MStable, mu)}, mc.cores = mcCores)
  newSd <- do.call(c, newSd)
#  newSd <- abs(2*log(newMu))   

  cat('\nFinal probes values...')
#   grpValues <- lapply(1:p, function(i){
#    cat(i,':\tFrom', meanValues[i], '\tto', newMu[i], '\tRatio:', 2^(newMu[i]-meanValues[i]),'\n')
#    values[i, Ig[i,]==1] <- rnorm(sum(Ig[i,]==1), newMu[i], newSd[i])
#    values[i,]
#  })
  
  grpValues <- mclapply(1:p, function(i){
    values[i, Ig[i,]==1] <- rnorm(sum(Ig[i,]==1), newMu[i], newSd[i])
    values[i,]
  }, mc.cores = mcCores)
  grpValues <- as.data.frame(do.call(rbind, grpValues))
  rownames(grpValues) <- paste0('signif', seq(1, nrow(grpValues)))
  colnames(grpValues) <- paste0('sample', seq(1, ncol(grpValues)))
  cat('\nDone.\n')
  return(list(Data = grpValues, grps = grps))
}

visualizeContruct <- function(eset){
  # Visualize the construction
  # S vs. M distribution
  oriProbes <- eset[-grep('random|signif', rownames(eset)),]
  M <- apply(oriProbes, 1, mean, na.rm = TRUE)
  S <- apply(oriProbes, 1, sd, na.rm = TRUE)
  tab <- .computeBounds(M, S)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
  boxplot(log10(tab$S) ~ factor(round(tab$M,2)), 
          outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
          xlab = 'Means', ylab = 'Log10(Sdev)', main = 'MS plot')
  
  # Visualize the added probes
  rProbes <- eset[grep('random', rownames(eset)),]
#  srProbes <- eset[grep('signif_rand', rownames(eset)),]
  sProbes <- eset[grep('signif[^_rand]', rownames(eset)),]
  rM <- apply(rProbes, 1, mean, na.rm = TRUE)
  rS <- apply(rProbes, 1, sd, na.rm = TRUE)
#  srM <- apply(srProbes, 1, mean, na.rm = TRUE)
#  srS <- apply(srProbes, 1, sd, na.rm = TRUE)
  sM <- apply(sProbes, 1, mean, na.rm = TRUE)
  sS <- apply(sProbes, 1, sd, na.rm = TRUE)
  
  smoothScatter(M, log10(S), colramp = colorRampPalette(c("white", "blue4")),
                xlab = 'Means', ylab = 'Log10(sDev)', main = 'Added probes')
  points(rM, log10(rS), pch = 19, cex = 0.5, col = 'grey75')
#  points(srM, log10(srS), pch = 19, cex = 0.75, col = 'cyan')
  points(sM, log10(sS), pch = 19, cex = 0.75, col = 'red3')
#   legend('bottomright', legend = c('original', 'random added', 'random grps added', 'main grps added'),
#          pch = 19, col = c('purple', 'grey70', 'cyan', 'red3'), cex = 1.5, ncol = 2, bty = 'n')
  legend('bottomright', legend = c('original', 'random added', 'grps added'),
         pch = 19, col = c('purple', 'grey70', 'red3'), cex = 1.5, bty = 'n')
  par(op)
}

WelchDf <- function(eset, grps){
  grps <- as.factor(grps)
  n1 <- apply(eset[, grps == levels(grps)[1]], 1, length)
  v1 <- apply(eset[, grps == levels(grps)[1]], 1, var)
  n2 <- apply(eset[, grps == levels(grps)[2]], 1, length)
  v2 <- apply(eset[, grps == levels(grps)[2]], 1, var)
  dfw <- (v1/n1+v2/n2)^2/(v1^2/(n1^2*(n1-1)) + v2^2/(n2^2*(n2-1)))
  return(dfw)
}