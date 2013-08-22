pcaInfo <- function(pcaScore){
  Fmax = pcaScore$lModel$top
  Fb = pcaScore$lModel$bottom
  xmid = pcaScore$lModel$xmid
  b = pcaScore$lModel$scal
  d = pcaScore$lModel$d
  informTable <- c()
  for(i in seq(0.05, 1, by = 0.05)){
    yTarg = (Fmax + Fb)*i
    xTarg = xmid - b*log(((Fmax - Fb)/(yTarg - Fb))^(1/d) - 1)
    
    aTarg <- 10^(xTarg)
    alpha = 10^(-aTarg)
    Pmax = 1 - alpha
    Q <- qchisq(p = Pmax, df = pcaScore$PCdim)  
    inform <- which(pcaScore$Dist >= Q)
    nInform <- length(inform)
    nNonInform <- pcaScore$m - nInform
    informTable <- rbind(informTable,
                         c(Prop = i, Inform = nInform, nonInform = nNonInform, propInform = nInform/pcaScore$m))
  }
  return(as.data.frame(informTable))
}
