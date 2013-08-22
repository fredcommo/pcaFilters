# PCA filetring: compute the performance of a classifier, using rejected or selected probes
# 2013-07-30

mainPath <- paste0('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_results_v5_Q')

#ent <- synGet('syn1898672')
#source(ent@filePath)

# parentId <- 'syn1930301'; Design <- 'TumourType' # Lee/lung Adk vs Sq
# parentId <- 'syn1939680'; Design <- 'disease_state' # Kruner/lung Adk vs Sq
# parentId <- 'syn1935134'; Design <- 'tissue_type' # Hou/lung
# parentId <- 'syn1939658'; Design <- 'Survival' # Smith/CRC
# parentId <- 'syn1939654'; Design <- 'subtype' # Shlicker/CRC
# parentId <- 'syn1939656'; Design <- 'Survival' # Jorissen/CRC
# parentId <- 'syn1939669'; Design <- 'Sample_characteristics_ch1' # Richardson/breast
# parentId <- 'syn1939676'; Design <- 'gg' # Metzer/breast
# parentId <- 'syn1939671'; Design <- 'ER' # Lu/breast
# parentId <- 'syn1939678'; Design <- 'er_status' # Design <- 'source_name_ch1' # Clarke/breast Tum vs Normal
# parentId <- 'syn2009712'; Design <- 'mfolfox6' # Tsuji/CRC
# parentId <- 'syn2010042'; Design <- 'tissue' # Khamas/CRC
# parentId <- 'syn1930524'; Design <- 'Condition' # Platinium

computeClassif <- function(parenId, Design)
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

  e <- synGet(parentId)
  resDir <- file.path(mainPath, paste0(propertyValue(e, 'name'), '_results_v5'))
  dir.create(resDir)
  setwd(resDir)

  query <- synapseQuery(paste0('select name, id from entity where parentId=="', parentId, '"'))
  query
  synId <- query$entity.id[grep('^GSE[0-9]+.rds', query$entity.name)]
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

  if(parentId == 'syn2010042'){
    Grps <- ifelse(grepl('normal', Grps), 'Normal',
                   ifelse(grepl('cancer', Grps), 'Cancer', 'CellLines'))
    Grps <- as.factor(Grps)
  }
  table(Grps)

  N <- nrow(eset)
  oProbes <- prcomp(eset)
  oScore <- pcaTraceQ(eset, oProbes, Plot = FALSE)
  S <- apply(eset, 1, sd, na.rm = TRUE)
  {
    if(nlevels(Grps)> 2){resp <- Grps; Family = 'multinomial'}
    else {resp <- c(0, 1)[Grps]; Family = 'binomial'}
  }

  L<- 15; B <- 15; mcCores = 2
  P <- exp(seq(log(.001), log(.95), len = L))

  pcaErr <- mclapply(P, function(p){
    cat('p:', p, '\t')
    select <- pcaSelectQ(oScore, p); n <- length(select)
    if(n > 2 & N-n > 1000)
      tmp <- lapply(1:B, function(i)
      {c(p = p, sel = n, rem = N-n, Classifier(eset, select, resp, Family, alpha = .01))})
    else
    {tmp <- list(); tmp[[length(tmp)+1]] <- c(p = p, sel = n, rem = N-n, c1 = NA, c2 = NA)}
    do.call(rbind, tmp)
  }, mc.cores = min(B, mcCores))
  pcaErr <- as.data.frame(do.call(rbind, pcaErr))

  #P <- seq(.025, .95, len = L)
  Rem <- unique(pcaErr$rem)
  varErr <- mclapply(Rem, function(n){
    p <- n/N
    select <- which(S > quantile(S, p, na.rm = TRUE)); n <- length(select)
    cat('p:', p, '\t', 'removed:', N-n, '\n')
    if(n > 2 & N-n > 1000)
      tmp <- lapply(1:B, function(i)
      {c(p = p, sel = n, rem = N-n, Classifier(eset, select, resp, Family, alpha = .01))})
    else
    {tmp <- list(); tmp[[length(tmp)+1]] <- c(p = p, sel = n, rem = N-n, c1 = NA, c2 = NA)}
    do.call(rbind, tmp)
  }, mc.cores = min(B, mcCores))
  varErr <- as.data.frame(do.call(rbind, varErr))
  
  png(paste0(synId, '_ErrorRate1.png'), width = 1300, height = 600)
    par(mfrow = c(1, 2), las = 1, mar = c(5, 5, 4, 2), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
    Xlab <- 'Removed probes'; Ylab <- 'Classification'
    plot(c(0, N), c(0.2, 1), type = 'n', xlab = Xlab, ylab = Ylab, main = 'Filtered by PCA')
    plotErr(pcaErr$c1, pcaErr$rem, pch = 19, col = 'blue')
    plotErr(pcaErr$c2, pcaErr$rem, pch = 19, col = 'red')
    legend('bottomright', legend = c('removed', 'selected'), pch = 19, col = c('red', 'blue'), cex = 1.25, bty = 'n')

    plot(c(0, N), c(0.2, 1), type = 'n', xlab = Xlab, ylab = Ylab, main = 'Filtered by VAR')
    plotErr(varErr$c1, varErr$rem, pch = 19, col = 'blue')
    plotErr(varErr$c2, varErr$rem, pch = 19, col = 'red')
    legend('bottomright', legend = c('removed', 'selected'), pch = 19, col = c('red', 'blue'), cex = 1.25, bty = 'n')
    par(op)
  dev.off()

  png(paste0(synId, '_ErrorRate2.png'), width = 800, height = 600)
    par(las = 1, mar = c(5, 5, 4, 2), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
    plot(c(0, N), c(0.2, 1), type = 'n', xlab = Xlab, ylab = Ylab, main = 'Classification using removed probes')
    plotErr(pcaErr$c2, pcaErr$rem, pch = 19, cex = 1.5, col = 'blue')
    plotErr(varErr$c2, varErr$rem, pch = 19, cex = 1.5, col = 'red')
    w <- wilcox.test(pcaErr$c2, varErr$c2)$p.value
    legend('bottomright', legend = c('byPCA', 'byVar'), pch = 19, col = c('blue', 'red'), cex = 1.5, bty = 'n')
    legend('bottomleft', legend = paste('Wilcoxon: p =', signif(w, 3)), cex = 1.5, bty = 'n')
    par(op)
  dev.off()
}
