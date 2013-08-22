# PCA pipeline2
# compute the frequences of excluded probes
require(synapseClient)
require(hgu133plus2.db)
require(VennDiagram)
# require(GO.db)
# require(preprocessCore)
source('~/Documents/MyProjects/FredScripts/ComBat3.R')
op <- par(no.readonly = TRUE)

getIntersect <- function(List, Type){
  idx <- which(tissues == Type)
  probes <- unlist(List)
  for(i in idx)
    probes <- intersect(probes, List[[i]])
  return(probes)
}

# Get a full matrix (Platinium data)
synId = 'syn1952964'
RDS <- synGet(synId)
Data <- get(load(RDS@filePath))
eset <- Data$eset
dim(eset)

#
setwd('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_results_v3_Q/')
fileList <- list.files()
fileList <- fileList[grep('.rds$', fileList)]
fileList <- fileList[-grep('Experiment', fileList)]
fileList <- fileList[-grep('Smith|Jorissen|Richardson', fileList)]
fileList
#
tissues <- lapply(fileList, function(f) unlist(strsplit(f, '_'))[1])
tissues <- do.call(c, tissues)
tissues
Names <- lapply(fileList, function(f) unlist(strsplit(f, '_'))[2])
Name <- do.call(c, Names)
Name
#
rejectList <- lapply(fileList, function(f){
  cat(f, '\n')
  result <- get(load(f))
  tmp <- result$rejectList$F0.1$pcaR
  tmp[-grep('random', tmp)]
})
names(rejectList) <- Names
fullList <- rownames(eset)

#
cat('Using PCA\n')
for(Type in c('Lung', 'Breast', 'CRC')){
  idx = which(tissues == Type)
  n <- length(idx)
  cat('Venn diagram:', Type, 'n:', n, '\n')
  nProbes <- length(unique(unlist(rejectList[idx])))
  venn.diagram(rejectList[idx], filename = paste0('./vennDiagrams/', Type, '_VennDiagram_F10.png'),
               fill = 2:(n+1), alpha = rep(0.1, n), cex = 1.5, cat.fontface = 4, lty = 1,
               main = paste0(Type, ': ', nProbes, ' rejected probes'), sub = '(% of total)',
               main.cex = 2, sub.cex = 1.75, useProp = TRUE, digits = 2)
  }

#
lungList <- getIntersect(rejectList, 'Lung'); length(lungList)
breastList <- getIntersect(rejectList, 'Breast'); length(breastList)
crcList <- getIntersect(rejectList, 'CRC'); length(crcList)
cat('Venn diagram all tissues\n')
nProbes <- length(unique(c(lungList, breastList, crcList)))
venn.diagram(list('Lung\n(n=3)' = lungList, 'Breast\n(n = 4)' = breastList, 'CRC(n = 3)' = crcList),
             filename = paste0('./vennDiagrams/', 'AllTissues_VennDiagram_F10.png'),
             main = paste0('All data sets: ', nProbes, ' rejected probes'),
             sub = '(% of total)', main.cex = 2, sub.cex = 1.75,
             fill = 2:4, alpha = rep(0.1, 3), useProp = TRUE, digits = 2,
             cex = 1.5, cat.fontface = 4, lty =1)

# On variance
cat('\nUsing VAR\n')
crit <- 'var50'
rejectList <- lapply(fileList, function(f){
  cat(f, '\n')
  result <- get(load(f))
  tmp <- result$rejectList$var50
  tmp[-grep('random', tmp)]
})
names(rejectList) <- Names
#
for(Type in c('Lung', 'Breast', 'CRC')){
  cat('Venn diagram:', Type, '\n')
  idx = which(tissues == Type)
  n <- length(idx)
  nProbes <- length(unique(unlist(rejectList[idx])))
  venn.diagram(rejectList[idx], fill = 2:(n+1),
               filename = paste0('./vennDiagrams/', Type, '_VennDiagram_', crit,'.png'),
               alpha = rep(0.1, n),
               cex = 1.5, cat.fontface = 4, lty = 1, fontfamily = 1,
               main = paste0(Type, ': ', nProbes, ' rejected probes'),
               sub = '(% of total)',
               main.cex = 2, sub.cex = 1.75,
               useProp = TRUE, digits = 2)
}
#
lungList <- getIntersect(rejectList, 'Lung'); length(lungList)
breastList <- getIntersect(rejectList, 'Breast'); length(breastList)
crcList <- getIntersect(rejectList, 'CRC'); length(crcList)
cat('Venn diagram all tissues\n')
nProbes <- length(unique(c(lungList, breastList, crcList)))
venn.diagram(list('Lung\n(n=3)' = lungList, 'Breast\n(n = 4)' = breastList, 'CRC(n = 3)' = crcList),
             filename = paste0('./vennDiagrams/', 'AllTissues_VennDiagram_', crit,'.png'),
             main = paste0('All data sets: ', nProbes, ' rejected probes'), 
             sub = '(% of total)', main.cex = 2, sub.cex = 1.75,
             fill = 2:4, alpha = rep(0.1, 3), useProp = TRUE, digits = 2,
             cex = 1.5, cat.fontface = 4, lty =1)

# Check annotations
# lungOnly <- setdiff(lungList, c(breastList, crcList))
# lungOnly <- select(hgu133plus2.db, lungOnly, c("SYMBOL", "GO"), "PROBEID")
# write.table(lungOnly, 'lungProbeList.txt', sep = '\t')
# lungTerms <- select(GO.db, lungOnly$GO, "TERM", "GOID")
# lungRejectList <- list(probeIds = lungOnly, GOterms = lungTerms)
# 
# # Chek levels
# lungOnly <- setdiff(lungList, c(breastList, crcList))
# breastOnly <- setdiff(breastList, c(lungList, crcList))
# crcOnly <- setdiff(crcList, c(breastList, lungList))
# common <- intersect(lungList, intersect(breastList, crcList))
# 
# synList <- c('syn1930301', 'syn1939680', 'syn1935134', 'syn1939658',
#              'syn1939654', 'syn1939656', 'syn1939669', 'syn1939676',
#              'syn1939671', 'syn1939678') #, 'syn1930524' : Platinium data
# rejectList <- lapply(synList, function(parentId){
#   e <- synGet(parentId)
#   e
#   query <- synapseQuery(paste0('select name, id from entity where parentId=="', parentId, '"'))
#   synId <- query$entity.id[grep('^syn[0-9]+(.*).rds', query$entity.name)]
#   RDS <- synGet(synId)
#   if(parentId == 'syn1930301')
#     Data <- readRDS(RDS@filePath)
#   else Data <- get(load(RDS@filePath))
# 
#   # Getting data
#   eset <- Data$eset
#   tmpList <- list(lung = eset[rownames(eset) %in% lungOnly,],
#                   breast = eset[rownames(eset) %in% breastOnly,],
#                   crc = eset[rownames(eset) %in% crcOnly,],
#                   common = eset[rownames(eset) %in% common,],
#                   Id = propertyValue(e, 'id'),
#                   Name = propertyValue(e, 'name'))
#   tmpList
#   }
# )
# 
# lungProbes <- breastProbes <- crcProbes <- commonProbes <- Batch <- labels <- c()
# for(tmp in rejectList){
#   cat(tmp$Name, '\t', tmp$Id, '\n')
#   cat(dim(tmp$lung), '\n')
#   cat(dim(tmp$breast), '\n')
#   cat(dim(tmp$crc), '\n')
#   cat(dim(tmp$common), '\n')
#   
#   lungProbes <- cbind(lungProbes, as.matrix(tmp$lung))
#   breastProbes <- cbind(breastProbes, as.matrix(tmp$breast))
#   crcProbes <- cbind(crcProbes, as.matrix(tmp$crc))
#   commonProbes <- cbind(commonProbes, as.matrix(tmp$common))
#   tmpBatch <- cbind(colnames(tmp$lung), tmp$Name)
#   Batch <- rbind(Batch, tmpBatch)
#   if(grepl('Lung', tmp$Name))
#     labels <- c(labels, c(rep('Lung', ncol(tmp$lung))))
#   if(grepl('Breast', tmp$Name))
#     labels <- c(labels, c(rep('Breast', ncol(tmp$breast))))
#   if(grepl('CRC', tmp$Name))
#     labels <- c(labels, c(rep('CRC', ncol(tmp$crc))))
#   }
# Batch <- as.data.frame(Batch)
# colnames(Batch) <- c('expIds', 'Batch')
# 
# X <- commonProbes
# 
# for(X in list(lungProbes, breastProbes, crcProbes, commonProbes)){
#   colnames(X) <- Batch$expIds
#   X <- ComBat3(X, Batch)
#   pcaX <- prcomp(t(X))
#   D2 <- rowSums(pcaX$x^2)
#   par(mfrow = c(2,2))
#   nCol = 10
#   plot(pcaX$x, col = c(1:3)[factor(labels)])
#   S <- sum(diag(var(pcaX$x[,1:nCol])))
#   legend('bottomright', legend = paste('Trace:', round(S, 2)), bty = 'n')
#   k = 1
#   for(tissue in c('Breast', 'CRC', 'Lung')){
#     tmp <- pcaX$x[labels == tissue,]
#     D2tmp <- rowSums(tmp^2)
#     Stmp <- sum(diag(var(tmp[,1:nCol])))
#     plot(tmp, col = k, xlim = range(pcaX$x[,1]), ylim = range(pcaX$x[,2]), main = tissue)
#     p <- pf(Stmp/S, nrow(tmp)-1, nrow(pcaX$x)-1)
#     leg1 <- paste('Trace:', round(Stmp, 2))
#     leg2 <- paste('p:', signif(p, 3) )
#     legend('bottomright', legend = c(leg1, leg2), bty = 'n')
#     k = k+1
#   }
#   par(op)
# }
