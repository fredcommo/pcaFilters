pcaTrace1.1 <- function(Data, PCA, Dim = 2:3,...){
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/Richard.w5PL.v2.R')
  X = scale(PCA$x[,Dim])
  X <- as.data.frame(X)
  D <- rowSums(X^2)
#  a.values <- 10^(seq(log10(1e-4), log10(1e-1), len = 12))
  a = 0.1
  for(i in 1:8) a <- c(a, a[length(a)]*2)
  a.values <- c(1e-4, 1e-3, 1e-2, a)
  
  cat('Testing\n')
  Score <- lapply(a.values, function(a){
    cat('\talpha:', a)
    alpha = 10^(-a)
    Pmax = 1 - alpha
    Q <- qchisq(p = Pmax, df = ncol(X))
    inform <- which(D <= Q)
    lInf <- length(inform)
    
    if(lInf>=2 & lInf<nrow(Data)){
      subData <- Data[inform,]
      acpTest <- prcomp(t(subData))
      tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
      tmpScore <- c(aValues = a, nProbes = lInf, Trace = tmpTrace)
      cat('\t#Probes:', lInf, '\tTrace:', tmpTrace, '\n')
      tmpScore
      }
    }
  )
  
  Score <- as.data.frame(do.call(rbind, Score))
  
  # Full matrix
  acpTest <- prcomp(t(Data))
  last.a <- Score$aValues[nrow(Score)]
  tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
  tmpScore <- cbind(aValues = c(last.a*2, last.a*4, last.a*8),
                    nProbes = rep(nrow(Data), 3),
                    Trace = rep(tmpTrace, 3))
  Score <- rbind(Score, tmpScore)
  
  cat('\n')
  rownames(Score) <- seq(1, nrow(Score))
  Score <- as.data.frame(Score)

  x <- as.numeric(log10(Score$aValues))
  y <- as.numeric(Score$Trace)
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]
    x <- x[-na.index]
  }
  y = y/(max(y)*1.01)*100
  if(any(y <=0)) y[y<=0] <- 1e-3
  model <- Richard.w5PL.v2(x, y, w = 0.25, Plot = TRUE, add.points = TRUE,
                            xlab = expression(-Log10(aValues)), ylab = 'Information (%)',...)  
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}
