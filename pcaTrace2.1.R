pcaTrace2.1 <- function(Data, PCA, Dim = 2:3, nv = 2,...){
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/Richard.w5PL.v2.R')
  #X = scale(PCA$x[,Dim])
  X = scale(PCA$u[,Dim])
  X <- as.data.frame(X)
  D <- rowSums(X^2)
  
  cat('Testing\n')
  Score <- lapply(1e-5*3^seq(0, 12), function(a){
    alpha = 10^(-a)
    Pmax = 1 - alpha
    Q <- qchisq(p = Pmax, df = ncol(X))
    inform <- which(D <= Q)
    lInf <- length(inform)
    tmpScore <- NA
    if(lInf>=nv & lInf<nrow(Data)){
      subData <- Data[inform,]
      svdTest <- fast.svd(t(subData))
      sValues <- svdTest$d[1:nv]
      sumValues <- sum(sqrt(sValues))
      tmpScore <- c(aValues = a, nProbes = lInf, D = sValues, sumD = sumValues)
      cat('alpha:', a, '\t#probes:', lInf, '\tsDev:', sumValues, '\n')
      tmpScore
      }
    }
  )
  
  Score <- as.data.frame(do.call(rbind, Score))
#  plot(seq(1,nrow(Score)), Score$sumD, type = 'l')
  
#   plot(c(1,nrow(Score)), range(Score[,3:(3+nv-1)]), type = 'n')
#   for(i in 3:(3+nv-1))
#     lines(seq(1, nrow(Score)), Score[,i])
  
  # Full matrix
  svdTest <- fast.svd(t(Data))
  sValues <- svdTest$d[1:nv]
  sumValues <- sum(sqrt(sValues))
  first.a <- Score$aValues[1]
  last.a <- Score$aValues[nrow(Score)]
  tmpScore <- cbind(aValues = last.a*3^seq(1, 2),
                    nProbes = rep(nrow(Data), 2),
                    D = matrix(rep(sValues, 2), 2, byrow = TRUE),
                    sumD = rep(sumValues, 2))
#  tmpStart <- cbind(aValues = first.a/3, nProbes = 0, D = t(rep(1e-3, nv)), sumValues = 0)
#  colnames(tmpScore) <- colnames(tmpStart) <- colnames(Score)
#  Score <- rbind(tmpStart, Score, tmpScore)
  colnames(tmpScore) <- colnames(Score)
  Score <- rbind(Score, tmpScore)
  cat('alpha:', last.a*3, '\tsDev:', sumValues, '\n')
  
  cat('\n')
  rownames(Score) <- seq(1, nrow(Score))
  Score <- as.data.frame(Score)

  x <- as.numeric(log10(Score$aValues))
  y <- as.numeric(Score$sumD)
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]
    x <- x[-na.index]
  }
  y = (y - min(y)*0.99)/(max(y)*1.01 - min(y)*0.99)
  if(any(y <=0)) y[y<=0] <- 1e-3
  model <- Richard.w5PL.v2(x, y, w = 0.25, Plot = TRUE, add.points = TRUE,
                            xlab = expression(-Log10(aValues)), ylab = 'Information (%)')#,...)
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}
