
.computeInformD <- function(D, q){
  # Returns probes in the rejected region: D<=q(a)
  cat('\tquantile:', q)
  return(which(D <= q))
}

.computeTraceD <- function(Data){
  # Returns the trace of the cov. matrix
  acpTest <- prcomp(Data)
  return(sum(acpTest$sdev^2))
}

pcaTraceD <- function(Data, PCA, Dim = 2:3, weight = 0.25, Plot = TRUE,...){
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
  avalues <-  quantile(D, probs = seq(.01, .99, len = 10))
  
  cat('Testing quantiles... ')
  Score <- mclapply(avalues, function(a){
    inform <- which(D<=a) #.computeInform(D, a, ncol(X))
    nInform <- length(inform)
    if(nInform>=2){
      tmpTrace <- .computeTraceD(t(Data[inform,]))
      tmpScore <- c(aValues = a, nProbes = nInform, Trace = tmpTrace)
      #cat('\t#Rejected:', length(inform), '\tTrace:', tmpTrace, '\n')
      tmpScore
    }
  }, mc.cores = 4)
  
  Score <- as.data.frame(do.call(rbind, Score))
  tmpTrace <- .computeTraceD(t(Data))
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
  model <- .RichardD(x, y, w = weight, Plot = Plot, add.points = TRUE,
                            xlab = expression(-Log10(quantile)), ylab = 'Information (%)',...)
  cat('Done.\n')
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}

pcaInfoD <- function(pcaScore){
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

pcaSelectD <- function(pcaScore, p = 0.05){
  # pcaScore: pcaTrace output
  # p: the proprtion of information required
  # Returns the indices of the slected features, according to p.
  if(p==0) return(1:pcaScore$m)
#   Fmax = pcaScore$lModel$top
#   Fb = pcaScore$lModel$bottom
#   xc = pcaScore$lModel$xc
#   b = pcaScore$lModel$scal
#   d = pcaScore$lModel$d
#   
#   #  yTarg = (Fmax + Fb)*p
#   #  xTarg = xmid - b*log(((Fmax - Fb)/(yTarg - Fb))^(1/d) - 1)
#   yTarg = Fb + (Fmax - Fb)*p
#   xTarg = xc - b*log10(((Fmax-Fb)/(yTarg - Fb))^(1/d) - 1)
  target <- .invRichardD(pcaScore$lModel, p)  
  xTarg <- target$xTarg
  aTarg <- 10^(xTarg)
  #   alpha = 10^(-aTarg)
  #   Pmax = 1 - alpha
  #   Q <- qchisq(p = Pmax, df = pcaScore$PCdim)
  #  inform <- which(pcaScore$Dist >= Q)
  inform <- which(pcaScore$Dist >= aTarg)
  cat('Removed:', pcaScore$m-length(inform),'\tSelected:', length(inform), '\n')
  return(inform)
}

.RichardD <- function(x, y, w = 0.25, Plot = FALSE, add.points = FALSE,
                             add.intercept = TRUE, pcol = "royalblue1", add.line = FALSE,
                             lcol = "navy", tan.col = "purple",...){ # Xlim = range(x), Ylim = range(y),
  
  # x      	: x-axis values
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

.invRichardD <- function(model, target){
  # target : proportion of signal to reach.
  Fmax = model$top
  Fb = model$bottom
  xc = model$xc
  b = model$scal
  d = model$d
  yTarg = Fb + (Fmax - Fb)*target
  xTarg = xc - b*log10(((Fmax-Fb)/(yTarg - Fb))^(1/d) - 1)
  return(list(xTarg = xTarg, yTarg = yTarg))
}