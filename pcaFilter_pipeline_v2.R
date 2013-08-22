# PCA filetring

############################
updateWiki <- function(synId , Infos){
  cat('Updating the wiki...\n')
  folderWikiUri <- sprintf("/entity/%s/wiki", synId)
  folderWiki <- synRestGET(folderWikiUri)
  folderUpdateWikiUri <- sprintf("/entity/%s/wiki/%s", synId, folderWiki$id)
  #  folderWiki$attachmentFileHandleIds[length(folderWiki$attachmentFileHandleIds)+1] <- file@fileHandle$id
  
  geoUrl <- paste0('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', Infos$geoId)
  addWiki <- paste0(Infos$synId, '|',
                    Infos$tissue, '|',
                    Infos$contrib, '|',
                    Infos$platform, '|',
                    '[', Infos$geoId, '](', geoUrl, ')', '|',
                    Infos$nSamples, '|',
                    Infos$Design, '\n')
  folderWiki$markdown <- paste0(folderWiki$markdown, addWiki)
  folderWiki <- synRestPUT(folderUpdateWikiUri, folderWiki)
  cat('Done\n')
  #  unlink(synapseCacheDir(), recursive = TRUE, force = TRUE)
  #  unlink(tempdir(), recursive = TRUE, force = TRUE)
}
############################

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
require(multicore)
require(corpcor)
require(synapseClient)

op <- par(no.readonly = TRUE)


#ent <- synGet('syn1898672')
#source(ent@filePath)

mainPath <- paste0('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_results_v4_Q')
# parentId <- 'syn1930301'; Design <- 'TumourType' # Lee
# parentId <- 'syn1939680'; Design <- 'disease_state' # Kruner
# parentId <- 'syn1935134'; Design <- 'tissue_type' # Hou
# parentId <- 'syn1939658'; Design <- 'Survival' # Smith
# parentId <- 'syn1939654'; Design <- 'subtype' # Shlicker
 parentId <- 'syn1939656'; Design <- 'NONE' # Jorissen
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

if(parentId == 'syn1939658'){ # special case for Smith
  dfs <- as.numeric(gsub('dfs_time: ', '', samples[,20]))
  Grps <- ifelse(dfs>median(dfs, na.rm = TRUE), 'high', 'low')
}

if(parentId == 'syn1939656'){
  dfs <- gsub('(.*)+DFS_Time: ', '', samples$Location)
  dfs <- as.numeric(gsub(';(.*)+', '', dfs))
  Grps <- ifelse(dfs>median(dfs, na.rm = TRUE), 'high', 'low')
}

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



# Adding significant probes to random groups
# SignifRand <- generateGrps(oriM, oriS, ncol(eset), nrow(eset)*.01,
#                            nGrp = 3, minP = 0.5, maxP = 0.9)
# signifRand <- SignifRand$Data
# rownames(signifRand) <- paste0('signif_rand', seq(1,nrow(signifRand)))
# randGrps <- SignifRand$grps

#colnames(Random) <- colnames(signifGrps) <- colnames(signifRand) <- colnames(eset)
colnames(Random) <- colnames(signifGrps) <- colnames(eset)
eset <- rbind(eset, Random, signifGrps)

# Visualize the construction
png(paste0(synId, '_Construction.png'), width = 1500, height = 800)
visualizeContruct(eset)
dev.off()

# PCA Filtering
pcaProbes <- prcomp(eset)
probeCols <- rep(rgb(.6, 0.3, .7, .1), nrow(eset))
probeCols[grep('random', rownames(eset))] <- rgb(.7, .7, .7, 1)

# Compute the traces according to axes
png(paste0(synId, '_multiTrace.png'), width = 900, height = 600)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
multiTrace(eset, pcaProbes)
par(op)
dev.off()

# Visualize the radii
png(paste0(synId, '_radPairs.png'), width = 900, height = 900)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
radPairs(pcaProbes$x[,1:5], p = .95, idx = grep('random', rownames(eset)),
         Pcol = probeCols, Cex = 1.5, cex = .2)
par(op)
dev.off()

# Roundness
Theta <- seq(0, 2*pi, len = 1000)
randidx <- grep('random', rownames(eset))
R <- lapply(1:10, function(i) roundness(pcaProbes$x[randidx,i:(i+1)], Theta, minN = 50))
R <- as.data.frame(do.call(cbind, R))

png(paste0(synId, '_roudness.png'), width = 900, height = 600)
par(mfrow = c(3,2), mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
n = 6
for(i in 1:n){
  plot(Theta, R[,i], col = rgb(i/n, 0.2, 1-i/n, 0.1), ylim = range(R[,1:n], na.rm = TRUE),
       ylab = 'Roundness', main = paste('Axes', i, 'to', i+1))
  polygon(c(0,Theta,max(Theta)), c(0,R[,i],0), col = rgb(i/n, 0.2, 1-i/n, 0.1))
  legend('topright', legend = paste('Err.Sq.:', round(sum(R[,i]^2, na.rm = TRUE))), cex = 1.2, bty = 'n')
}
par(op)
dev.off()

png(paste0(synId, '_roudnessBoxplot.png'), width = 900, height = 600)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
boxplot(R, outcex = .1); abline(h = 0, col = 'red')
par(op)
dev.off()

# Compute the filtering performance depending on the PC axes:
trueList <- ifelse(grepl('signif', rownames(eset)), 'signif',
                   ifelse(grepl('random', rownames(eset)), 'random', 'original'))
trueList <- factor(trueList)
pcaScore <- pcaPerf(eset, pcaProbes, Grps, trueList, threshold = 1e-3)

# Visualize performances:
png(paste0(synId, '_VisualizePerf.png'), width = 1200, height = 1200)
par(mfrow = c(2, 2), cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, 
    las = 1, mar = c(5, 6.5, 4, 2))
plotScore(pcaScore)
plotScore(pcaScore, type = 'signif')
plotScore(pcaScore, type = 'score')
par(op)
dev.off()

# Keep the best dimension
pcaScore[which.max(pcaScore$score),]
score <- pcaTraceQ(eset, pcaProbes, Dim = 2:3)
for(p in c(.025, .05, .1, .2, .5)){
  png(paste0(synId, '_VisualizeFilter_', p, '.png'), width = 1200, height = 1200)
  plotRemoved(eset, score = score, p = p, Grps = Grps)
  dev.off()
}

Info <- pcaInfoQ(score); Info
png(paste0(synId, '_VisualizeSamples.png'), width = 1200, height = 1200)
visualizeSamples(eset, score, pch = 19, col = c('orangered', 'darkblue', 'green4', 'goldenrod3')[factor(Grps)])
par(op)
dev.off()

png(paste0(synId, '_VisualizeRemove_0..png'), width = 1200, height = 1200)
par(cex.main = 2.5, cex.lab = 2, cex.axis = 2)
visualizeProbes(eset, pcaProbes, p = 0.0)
dev.off()
png(paste0(synId, '_VisualizeRemove_0.025.png'), width = 1200, height = 1200)
par(cex.main = 2.5, cex.lab = 2, cex.axis = 2)
visualizeProbes(eset, pcaProbes, p = 0.025)
dev.off()
png(paste0(synId, '_VisualizeRemove_0.05.png'), width = 1200, height = 1200)
par(cex.main = 2.5, cex.lab = 2, cex.axis = 2)
visualizeProbes(eset, pcaProbes, p = 0.05)
dev.off()
png(paste0(synId, '_VisualizeRemove_0.1.png'), width = 1200, height = 1200)
visualizeProbes(eset, pcaProbes, p = 0.1)
dev.off()
png(paste0(synId, '_VisualizeRemove_0.2.png'), width = 1200, height = 1200)
visualizeProbes(eset, pcaProbes, p = 0.2)
par(op)
dev.off()

# Visualize PCA selection depending on axes
#visualizeProbes(eset, pcaProbes, p = .2)

# Compute the filtering performance depending on qth variance:
varScore <- varPerf(eset, seq(.1, .9, by = .1), Grps, trueList)
png(paste0(synId, '_VisualizeVarPerf.png'), width = 600, height = 600)
par(cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, 
    las = 1, mar = c(5, 6.5, 4, 2))
plotVarPerf(varScore)
par(op)
dev.off()

# ROC curve PCA-filter Vs. Var-filter
# Compare the filtering perf: S-filter quantile defined such that the same number of probes is filtered.
#compareScore <- compareScores(eset, score, Grp = Grps, P = seq(0, .95, by = 0.05))
newScore <- pcaTraceQ(eset, pcaProbes, Dim = 2:3)
compareScore <- compareScores(eset, newScore, Grp = Grps, P = seq(0, .95, by = 0.05), type = 'useChi2')
compareScore
nRandom <- sum(grepl('random', rownames(eset)))
nSignif <- sum(grepl('signif', rownames(eset)))

png(paste0(synId, '_VisualizeComparePerf.png'), width = 1300, height = 1200)
par(mfrow = c(2,2),mar = c(5, 6, 4, 2.5), cex.main = 2.5, cex.lab = 2, cex.axis = 1.5,
    las = 1)
plotCompare(compareScore, nRandom, nSignif)
plotCompare(compareScore, nRandom, nSignif, type = 'signif')
plotCompare(compareScore, nRandom, nSignif, type = 'score')
plotROC(compareScore, nRandom, nSignif)
par(op)
dev.off()


# Visualize selection S/M densities
M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)
select <- pcaSelectQ(score, .025)
q <- quantile(S, probs = 1 - length(select)/nrow(eset))

png(paste0(synId, '_VisualizeCVDistribution.png'), width = 1300, height = 1200)
par(mar = c(5, 6, 4, 2.5), cex.main = 2.5, cex.lab = 2.25, cex.axis = 2, las = 1)
GMmodel(log(S/M), resamp = 1e5, xlab = expression(Log10(S/M)),
        main = 'Coefficient of variation')
SM <- log(S/M)
SMpca <- ifelse(seq(1, length(SM)) %in% select, SM, max(SM)+1)
SMvar <- ifelse(seq(1, length(SM)) %in% which(S>q), SM, max(SM)+1)
hist(SMpca, nclass = 100, add = TRUE, prob = T, border = rgb(0,0,1,0.5), col = rgb(0,0,1,0.25))
hist(SMvar, nclass = 100, add = TRUE, prob = T, border = rgb(1,0,0,0.5), col = rgb(1,0,0,0.25))
lines(density(SMpca), lwd = 6, col = 'blue3')
lines(density(SMvar), lwd = 6, col = 'red3')
legend('topleft', title = 'Filtered', legend = c('by PCA', 'by Var'),
       lwd = 5, cex = 2, col = c('blue3', 'red3'), bty = 'n')
par(op)
dev.off()

output <- saveFilters(eset, score)
Results <- list(nSamples = ncol(eset),
                Design = Grps,
                pcaPerf = pcaScore,
                varPerf = varScore,
                compareScore = compareScore,
                rejectList = output)
#savePath <- '/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_results/'
filePath <- paste0(mainPath, '/',propertyValue(e, 'name'), '_Results_SurvMed_Q.rds')
save(Results, file = filePath)
cat('Pushing to synapse...\n')
file <- File(filePath, parentId = parentId)
file <- synStore(file,
                 activityName = 'PCA-based Filtering',
                 used = list(list(entity = synGet('syn1960969'), wasExecuted = TRUE),
                             list(entity = RDS, wasExecuted = FALSE))
                 )

# Name <- unlist(strsplit(propertyValue(e, 'name'), '_'))
# Design <-lapply(1:nlevels(Grps), function(i) levels(Grps)[i])
# Design <- do.call(paste, Design)
# Design
# Infos <- cbind.data.frame(synId = synId, tissue = Name[1], contrib = Name[2], geoId = Name[3],
#                           platform = 'hgu133plus2', nSamples = ncol(eset),
#                           Design = Design)
# updateWiki('syn1898587', Infos)
