pcaSelect <- function(pcaScore, p = 0.5){
  Fmax = pcaScore$lModel$top
  Fb = pcaScore$lModel$bottom
  xmid = pcaScore$lModel$xmid
  b = pcaScore$lModel$scal
  d = pcaScore$lModel$d
  
  yTarg = (Fmax + Fb)*p
  xTarg = xmid - b*log(((Fmax - Fb)/(yTarg - Fb))^(1/d) - 1)
  
  aTarg <- 10^(xTarg)
  alpha = 10^(-aTarg)
  Pmax = 1 - alpha
  Q <- qchisq(p = Pmax, df = pcaScore$PCdim)
  inform <- which(pcaScore$Dist >= Q)
  return(inform)
}
