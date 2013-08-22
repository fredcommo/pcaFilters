pcaTrace <- function(Data, PCA, Dim = 2:3,...){
  source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/Richard.w5PL.v2.R')
  X = scale(PCA$x[,Dim])
  X <- as.data.frame(X)
  D <- rowSums(X^2)
  a.values <- 10^(seq(log10(1e-4), log10(32), len = 10))
  
  useTrace = T; Score <- c()
  cat('Testing a = ')
#  plot(X)
  Score <- lapply(1:(length(a.values)-1), function(i){
    cat('\t', a.values[i])
    if(i == 1) Pmin = 0
    else Pmin = 1 - 10^(-a.values[i])
    Pmax = 1 - 10^(-a.values[i+1])
    Q <- qchisq(p = c(Pmin, Pmax), df = ncol(X))
    inform <- which(D>=Q[1] & D < Q[2])
    lInf <- length(inform)
    points(X[inform,1], X[inform,2], col = i, pch = 0.1)
    #tmpScore <- NA    
    if(lInf>1){
      tmpNselect <- length(inform)
      subData <- Data[inform,]
      tmpS2 <- var(as.vector(t(subData)))
      tmpScore <- c(aValues = a.values[i], nProbes = tmpNselect, S2 = tmpS2)
      
      # calculer les variances sur les axes1 et 2 de la nouvelle acp ?
      if(useTrace){
        acpTest <- prcomp(t(subData))
        tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
        tmpScore <- c(tmpScore, Trace = tmpTrace)
        }
      tmpScore
      }
    }
  )
  Score <- as.data.frame(do.call(rbind, Score))
    
  cat('\n')
  rownames(Score) <- seq(1, nrow(Score))
  Score <- as.data.frame(Score)

  x <- as.numeric(log10(Score$aValues))
  y <- as.numeric(cumsum(Score$Trace))
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]
    x <- x[-na.index]
  }
  y = y/max(y)*100
  if(any(y<=0)) y[y<=0] <- 1e-3
  if(any(y>=100)) y[y>=100] <- 99.999
  model <- Richard.w5PL.v2(x, y, w = 0.25, Plot = TRUE, add.points = TRUE,
                            xlab = expression(-Log10(aValues)), ylab = 'Information (%)')#,...)  
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}
