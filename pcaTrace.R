pcaTrace <- function(Data, PCA, Dim = 2:3,...){
  require(corpcor)
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/Richard.w5PL.v2.R')
  X = scale(PCA$x[,Dim])
  X <- as.data.frame(X)
  D <- rowSums(X^2)
  a = 0.1
  for(i in 1:8) a <- c(a, a[length(a)]*2)
  a.values <- c(1e-4, 1e-3, 1e-2, a)
  #par(mfrow = c(3,2))
  #graphValues <- c(1e-3, 0.1, 0.8)
  useTrace = T; Score <- c()
  cat('Testing\n')
  for(a in a.values){
    cat('\talpha: ', a)
    alpha = 10^(-a)
    Pmax = 1 - alpha
    Q <- qchisq(p = Pmax, df = ncol(X))
    inform <- which(D <= Q)
    lInf <- length(inform)
    tmpScore <- NA
#    if(a %in% graphValues)
#      plot(X, col = ifelse(D<=Q, 'red3', 'grey'), cex = 0.1, main = length(inform))
    
    if(lInf>1 & lInf<nrow(Data)){
      tmpNselect <- length(inform)
      subData <- Data[inform,]
      tmpS2 <- var(as.vector(t(subData)))
      tmpScore <- c(aValues = a, nProbes = tmpNselect, S2 = tmpS2)
      
      # calculer les variances sur les axes1 et 2 de la nouvelle acp ?
      if(useTrace){
        acpTest <- prcomp(t(subData))
        tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
        tmpScore <- c(tmpScore, Trace = tmpTrace)
#         if(a %in% graphValues){
#           PC3 <- acpTest$x[,3]
#           mp = min(PC3)
#           Mp = max(PC3)
#           pcex = 2*(PC3-mp)/(Mp - mp) + 0.5
#           plot(acpTest$x, cex = pcex, xlim = range(-100, 100),
#                pch = 19, col = rgb(0.5, 0, 0.5, 0.5), main = paste('Trace: ', round(tmpTrace, 2)))
#         }
        cat('\t#Probes:', lInf, '\tTrace:', tmpTrace, '\n')
      }
      Score <- rbind(Score, tmpScore)
      }
  }
#  par(mfrow = c(1,1))
  
  # Full matrix
  acpTest <- prcomp(t(Data))
  tmpS2 <- var(as.vector(t(Data)))
  tmpScore <- c(aValues = a*2, nProbes = nrow(Data), S2 = tmpS2)
  tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
  tmpScore <- c(tmpScore, Trace = tmpTrace)
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
  y = y/max(y)*100
  if(any(y <=0)) y[y<=0] <- 1e-3
  model <- Richard.w5PL.v2(x, y, w = 0.25, Plot = TRUE, add.points = TRUE,
                            xlab = expression(-Log10(aValues)), ylab = 'Information (%)',...)  
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}
