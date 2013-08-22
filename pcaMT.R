# Special

#source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaFilterD.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaFilterQ.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaFilterHelpers.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/simulationFunctions.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/graph3D.8.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/heatmap.3.R')
require(multtest)
require(glmnet)
require(VennDiagram)
require(mclust)
require(graphics)
require(parallel)
require(snow)
require(multicore)
require(corpcor)
require(synapseClient)

op <- par(no.readonly = TRUE)


#########################
mt1 <- function(p, score, eset){
  cat('p:', p, '\t')
  select <- pcaSelectQ(score, p)
  mt <- try(mt.maxT(eset[select,], classlabel=Grps, B = 2500), silent = TRUE)
  nSignif <- sum(mt$adjp<=1e-3)
  trueSignif <- sum(grepl('signif', rownames(eset)))
  n <- sum(mt$adjp<=1e-3 & grepl('signif', rownames(mt)))
  return(c(n = length(select), p1 = nSignif/nrow(eset), p2 = nSignif/nrow(mt), p3 = n/trueSignif))
}

mt2 <- function(p, S, eset){
  cat('p:', p, '\t')
  select <- which(S > quantile(S, p))
  mt <- try(mt.maxT(eset[select,], classlabel=Grps, B = 2500), silent = TRUE)
  nSignif <- sum(mt$adjp<=1e-3)
  trueSignif <- sum(grepl('signif', rownames(eset)))
  n <- sum(mt$adjp<=1e-3 & grepl('signif', rownames(mt)))
  return(c(n = length(select), p1 = nSignif/nrow(eset), p2 = nSignif/nrow(mt), p3 = n/trueSignif))
}
#########################

#parentId <- 'syn2009712'; Design <- 'lesion' #'mfolfox6' # Tsuji/CRC
parentId <- 'syn1930301'; Design <- 'TumourType' # Lee/lung Adk vs Sq

e <- synGet(parentId)

query <- synapseQuery(paste0('select name, id from entity where parentId=="', parentId, '"'))
query
synId <- query$entity.id[grep('^GSE[0-9]+.rds', query$entity.name)]
synId
RDS <- synGet(synId)

# Use readRDS instead for parentId = 'syn1930301' (kimData)
#Data <- get(load(RDS@filePath))
Data <- readRDS(RDS@filePath)

# Getting data and samples
eset <- Data$eset
samples <- Data$samples
oriM <- apply(eset, 1, mean, na.rm = TRUE)
oriS <- apply(eset, 1, sd, na.rm = TRUE)
Grps <- samples[,which(colnames(samples) == Design)]

output <- mclapply(1:20, function(i){
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaFilterQ.R')
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaFilterHelpers.R')
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/simulationFunctions.R')
  cat('Run:', i, '\t')
  SignifGrps <- generateGrps(oriM, oriS, ncol(Data$eset), nrow(Data$eset)*.01,
                               grps = Grps, minP = 0.5, maxP = 0.9)
  signifGrps <- SignifGrps$Data
  Grps <- SignifGrps$grps

  # colnames(Random) <- colnames(signifGrps) <- colnames(signifRand) <- colnames(eset)
  colnames(signifGrps) <- colnames(Data$eset)
  eset <- rbind(Data$eset, signifGrps)

  pcaProbes <- prcomp(eset)
  score <- pcaTraceQ(eset, pcaProbes, Plot = FALSE)
  cat('Running test on PCA...\n')
  P <- seq(0.5, .90, len = 3)
  pow1 <- lapply(P, function(p) {return(mt1(p, score, eset))})
  pow1 <- cbind.data.frame(p = P, do.call(rbind, pow1))

  cat('Running test on Var...\n')
  S <- apply(eset, 1, sd, na.rm = TRUE)
  P <- 1 - pow1$n/nrow(eset)
  pow2 <- lapply(P, function(p) {return(mt2(p, S, eset))})
  pow2 <- cbind.data.frame(p = P, do.call(rbind, pow2))
  return(list(pow1 = pow1, pow2 = pow2))
  cat('\n')
  }, mc.cores = 2)
pow1 <- do.call(rbind, lapply(1:length(output),function(i) output[[i]]$pow1))
pow2 <- do.call(rbind, lapply(1:length(output),function(i) output[[i]]$pow2))


plot(factor(pow1$p), pow1$p1, border = rgb(0, .2, .8, .25), outcol = 'blue', outpch = '+',
     names = round(unique(pow1$p), 2))
boxplot(pow2$p1 ~ factor(pow1$p), border = rgb(0.8, .2, 0, .25), outcol = 'red', outpch = '+',
        add = TRUE, names = NA)


select <- pcaSelectQ(score, .85)

pcaSamples <- prcomp(t(eset[select,]))
plot(pcaSamples$x, col = c(1:4)[Grps])

Cols <- c('blue', 'red', 'green', 'yellow')[Grps]
heatmap.3(eset[select,], density.info = 'none', trace = 'none',
          col = colorpanel(100, "darkblue", "grey95", "yellow"),
          scale = 'row', ColSideColors = Cols)


