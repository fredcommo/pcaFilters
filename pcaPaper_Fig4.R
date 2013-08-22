# PCA filetring: special script for Fig4

source('/Users/fredcommo/Documents/MyProjects/ProjetACP/ScriptsPCAfilter/pcaFilterQ.R')
source('/Users/fredcommo/Documents/MyProjects/ProjetACP/ScriptsPCAfilter/pcaFilterHelpers.R')
source('/Users/fredcommo/Documents/MyProjects/ProjetACP/ScriptsPCAfilter/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/ProjetACP/ScriptsPCAfilter/simulationFunctions.R')
source('/Users/fredcommo/Documents/MyProjects/ProjetACP/ScriptsPCAfilter/graph3D.8.R')
source('/Users/fredcommo/Documents/MyProjects/ProjetACP/ScriptsPCAfilter/heatmap.3.R')
require(multtest)
require(glmnet)
require(VennDiagram)
require(mclust)
require(graphics)
require(parallel)
require(multicore)
require(corpcor)
require(synapseClient)

op <- par(no.readonly = TRUE)


#ent <- synGet('syn1898672')
#source(ent@filePath)

mainPath <- paste0('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_results_v3_Q')
parentId <- 'syn1930301'; Design <- 'TumourType' # Lee
# parentId <- 'syn1939680'; Design <- 'disease_state' # Kruner
# parentId <- 'syn1935134'; Design <- 'tissue_type' # Hou
# parentId <- 'syn1939658'; Design <- 'NONE' # Smith
# parentId <- 'syn1939654'; Design <- 'subtype' # Shlicker
# parentId <- 'syn1939656'; Design <- 'NONE' # Jorissen
# parentId <- 'syn1939669'; Design <- 'Sample_characteristics_ch1' # Richardson
# parentId <- 'syn1939676'; Design <- 'gg' # Metzer
# parentId <- 'syn1939671'; Design <- 'ER' # Lu
# parentId <- 'syn1939678'; Design <- 'Sample_source_name_ch1' # CLarke
# parentId <- 'syn1930524'; Design <- 'Condition' # Platinium

e <- synGet(parentId)
resDir <- file.path(mainPath, paste0(propertyValue(e, 'name'), '_results_v3'))
dir.create(resDir)
setwd(resDir)

query <- synapseQuery(paste0('select name, id from entity where parentId=="', parentId, '"'))
query
synId <- query$entity.id[grep('^syn[0-9](.*)+.rds', query$entity.name)]
synId

RDS <- synGet(synId)
# Use readRDS instead for parentId = 'syn1930301' (kimData)
{
  if(parentId == 'syn1930301') Data <- readRDS(RDS@filePath)
  else Data <- get(load(RDS@filePath))
}

# Getting data and samples
eset <- Data$eset
samples <- Data$samples
#samples <- read.csv('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_Jorissen_GSE14333/GSE14333_samples_FC.txt',
#                     header = T, sep = '\t')
if(Design!='NONE')
  Grps <- samples[,which(colnames(samples)==Design)]
colnames(eset) <- gsub('\\..*|_.*', '', colnames(eset))
all(colnames(eset) == as.character(samples$Sample_geo_accession))

# Original M & S values
oriM <- apply(eset, 1, mean, na.rm = TRUE)
oriS <- apply(eset, 1, sd, na.rm = TRUE)

# Generate random probes and signif probes
Random <- generateRandom(eset, floor(nrow(eset)*.1))

# Adding significant probes to main groups
{
  if(Design == 'NONE')
    SignifGrps <- generateGrps(oriM, oriS, ncol(eset), nrow(eset)*.01,
                               nGrp = 2, minP = 0.5, maxP = 0.9)
  else
    SignifGrps <- generateGrps(oriM, oriS, ncol(eset), nrow(eset)*.01,
                               grps = Grps, minP = 0.5, maxP = 0.9)
}
signifGrps <- SignifGrps$Data
Grps <- SignifGrps$grps
colnames(Random) <- colnames(signifGrps) <- colnames(eset)
eset <- rbind(eset, Random, signifGrps)

# PCA Filtering
pcaProbes <- prcomp(eset)
probeCols <- rep(rgb(.6, 0.3, .7, .1), nrow(eset))
probeCols[grep('random', rownames(eset))] <- rgb(.7, .7, .7, 1)
pcaProbes <- prcomp(eset)
score <- pcaTraceQ(eset, pcaProbes, Plot = FALSE)
pcaIni <- prcomp(t(eset))
Xlim <- range(pcaIni$x[,1]); Ylim <- range(pcaIni$x[,2])

M <- matrix(c(1:8, rep(9, 4)), 3, 4, byrow = T)
nf <- layout(M, c(1, 1, 1, 1), c(.25, .25, .5))
#layout.show(nf)

P <- c(.05, .1, .25, .5)
# part 1-4
par(mar = c(2, 5, 2, 2), cex.lab = 1.5)
for(p in P){
  select <- pcaSelectQ(score, p)
  Cols <- rep('red', nrow(eset))
  Cols[select] <- 'grey80'
#  Title <- paste0('Q', p, ' - ', nrow(eset)-length(select), ' probes')
  plot(pcaProbes$x[,2:3], cex = .2, col = Cols, main = paste0('Q', p))
  legend('topleft', legend = paste(nrow(eset)-length(select), ' probes'),
         pch = 19, col = 'red', cex = 1.25, bty = 'n')
}

# part 5- 8
for(p in P){
  select <- pcaSelectQ(score, p)
  pcaSamples <- prcomp(t(eset[-select,]))
  plot(pcaSamples$x, col = c('darkblue', 'orangered')[Grps], xlim = Xlim, ylim = Ylim)
}

#part9
par(mar = c(5, 15, 2, 12), cex.lab = 2, cex.axis = 1.5, las = 1)
targets <- lapply(P, function(p) .invRichardQ(score$lModel, p))
targets <- as.data.frame(do.call(rbind, targets))
score <- pcaTraceQ(eset, pcaProbes, lwd = 4)
legend('topleft', legend = 'Information curve (5PL)', cex = 1.75, bty = 'n')
for(i in 1:nrow(targets)){
  points(targets$xTarg[i], targets$yTarg[i], cex = 2, col = 'red', lwd = 5)
  text(as.numeric(targets$xTarg[i])-.25, as.numeric(targets$yTarg[i])+5,
       labels = paste0(P[i]*100, '%'), pos = 3, cex = 1.75, font = 2)
}
par(op)

