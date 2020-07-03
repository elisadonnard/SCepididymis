## modified inferNonExpressedGenes
inferNonExpressedGenes = function(scl){
  if(is(scl,'SoupChannelList')){
    for(i in seq_along(scl$channels)){
      message(sprintf("Inferring non-expressed genes for channel %s",names(scl$channels)[i]))
      scl$channels[[i]] = inferNonExpressedGenes(scl$channels[[i]])
    }
  }else if(is(scl,'SoupChannel')){
    #Construct expression fractions from table of counts
    rat = t(t(scl$toc)/scl$nUMIs)
    #Convert to useful format
    rat = as(rat,'dgTMatrix')
    # Fixed to add row and column names (AD 6/11/18)
    rownames(rat) <- rownames(scl$toc)
    colnames(rat) <- colnames(scl$toc)
    #Now convert to ratio
    rat@x = rat@x/scl$soupProfile[rat@i+1,'est']
    #Exclude the useless rows
    toKeep = unique(rat@i+1)
    toKeep = toKeep[scl$soupProfile[toKeep,'est']>0]
    rat = rat[toKeep,]
    #Get summary stats
    nCells = table(factor(rownames(rat)[rat@i+1],levels=rownames(rat)))
    lowCount = table(factor(rownames(rat)[rat@i[rat@x<1]+1],levels=rownames(rat)))
    tmp = split(rat@x,rownames(rat)[rat@i+1])[names(nCells)]
    extremity = sapply(tmp,function(e) sum(log10(e)**2)/length(e))[names(lowCount)]
    centrality = sapply(tmp,function(e) sum(1/(1+log10(e)**2))/length(e))[names(lowCount)]
    minFrac = log10(sapply(tmp,min))
    #Construct targets data frame
    targets = data.frame(nCells = as.integer(nCells),
                         lowCount = as.integer(lowCount),
                         lowFrac = as.integer(lowCount)/as.integer(nCells),
                         extremity = extremity,
                         centrality = centrality,
                         minFrac = minFrac,
                         isUseful = as.integer(lowCount)>0.1*ncol(rat))
    #Order by something
    targets = targets[order(targets$isUseful,targets$extremity,decreasing=TRUE),]
    scl$nonExpressedGenes = targets
  }else{
    stop("scl must be a SoupChannel or SoupChannelList object")
  }
  return(scl)
}

### modified adjustCounts
adjustCounts = function (scl, samp, pCut = 0.01, verbose = TRUE) 
{
  if (is(scl, "SoupChannel")) {
    sc = scl
    tmp = list(sc)
    names(tmp) = samp
    scl = SoupChannelList(tmp)
  }
  genes = rownames(scl$toc)
  cells = colnames(scl$toc)
  scl$toc = as(scl$toc, "dgTMatrix")
  # Fixed to add row and column names (ED 8/13/19)
  rownames(scl$toc) = genes
  colnames(scl$toc) = cells
  #
  rhos = unlist(lapply(scl$channels, function(e) e$rhos), 
                use.names = FALSE)
  names(rhos) = unlist(lapply(scl$channels, function(e) names(e$rhos)), 
                       use.names = FALSE)
  rhos = rhos[colnames(scl$toc)]
  if (verbose) 
    message("Calculating probability of each gene being soup")
  cMap = match(gsub("___.*", "", colnames(scl$toc)), colnames(scl$soupMatrix))
  p = pbinom(scl$toc@x - 1, scl$nUMIs[scl$toc@j + 1], 
             rhos[scl$toc@j + 1] * 
               scl$soupMatrix[cbind(scl$toc@i + 1, cMap[scl$toc@j + 1])], 
             lower.tail = FALSE)
  o = order(-(scl$toc@j + 1), p, decreasing = TRUE)
  if (verbose) 
    message("Calculating probability of the next count being soup")
  s = split(o, scl$toc@j[o] + 1)
  rTot = unlist(lapply(s, function(e) cumsum(scl$toc@x[e])), 
                use.names = FALSE)
  pSoup = pbinom(rTot - scl$toc@x[o] - 1, 
                 scl$nUMIs[scl$toc@j[o] + 1],
                 rhos[scl$toc@j[o] + 1], 
                 lower.tail = FALSE)
  if (verbose) 
    message("Filtering table of counts")
  pp = p[o] * pSoup
  w = which(pp < pCut)
  dropped = data.frame(cell = colnames(scl$toc)[scl$toc@j[o[-w]] + 
                                                  1], gene = rownames(scl$toc)[scl$toc@i[o[-w]] + 1], 
                       channel = colnames(scl$soupMatrix)[cMap[scl$toc@j[o[-w]] + 
                                                                 1]])
  if (verbose) {
    for (channel in colnames(scl$soupMatrix)) {
      message(sprintf("Most removed genes for channel %s are:", 
                      channel))
      x = sort(table(dropped$gene[dropped$channel == channel])/sum(colnames(scl$soupMatrix)[cMap] == 
                                                                     channel), decreasing = TRUE)
      print(x[seq_len(min(length(x), 100))])
    }
  }
  scl$atoc = sparseMatrix(i = scl$toc@i[o[w]] + 1, j = scl$toc@j[o[w]] + 
                            1, x = scl$toc@x[o[w]], dims = dim(scl$toc), dimnames = dimnames(scl$toc))

  scl$channels[[samp]]$atoc = as.matrix(scl$atoc)
  sc = scl$channels[[samp]]
  return(sc)
}


### modified to print gene name as title (ED)
plotSoupX = function (scl, geneSet, DR, dataType = c("binary", "counts", 
                                                     "expression", "ratio"), logData = TRUE, includePanels = c("Uncorrected", 
                                                                                                               "CorrectedExpression", "CorrectedCounts")) 
{
  dataType = match.arg(dataType)
  if (dataType == "binary") 
    logData = FALSE
  if (dataType == "ratio") 
    logData = TRUE
  if (is.null(scl$strainedExp) & is.null(scl$atoc)) 
    stop("No adjusted expression matrix found.")
  DR = as.data.frame(DR)
  if (ncol(DR) < 2) 
    stop("Need at least two reduced dimensions.")
  if (!(all(rownames(DR) %in% colnames(scl$toc)))) 
    stop("rownames of DR need to match column names of scl$toc")
  dfs = list()
  df = DR
  df$correction = "Uncorrected"
  if (dataType == "binary") {
    df$data = colSums(scl$toc[geneSet, rownames(df), drop = FALSE]) > 
      0
  }
  else if (dataType == "counts") {
    df$data = colSums(scl$toc[geneSet, rownames(df), drop = FALSE])
  }
  else {
    df$data = colSums(scl$toc[geneSet, rownames(df), drop = FALSE])/scl$nUMIs[rownames(df)]
  }
  if (dataType == "ratio") {
    df$data = scale(df$data)[, 1]
  }
  else {
    if (logData) 
      df$data = log10(df$data)
  }
  dfs[["raw"]] = df
  if (!is.null(scl$strainedExp)) {
    df = DR
    df$correction = "CorrectedExpression"
    if (dataType == "binary") {
      df$data = colSums(scl$strainedExp[geneSet, rownames(df), 
                                        drop = FALSE]) > 0
    }
    else if (dataType == "counts") {
      df$data = colSums(scl$strainedExp[geneSet, rownames(df), 
                                        drop = FALSE]) * 10000
    }
    else if (dataType == "ratio") {
      df$data = colSums(scl$strainedExp[geneSet, rownames(df), 
                                        drop = FALSE])/(colSums(scl$toc[geneSet, rownames(df), 
                                                                        drop = FALSE])/scl$nUMIs[rownames(df)])
    }
    else {
      df$data = colSums(scl$strainedExp[geneSet, rownames(df), 
                                        drop = FALSE])
    }
    if (logData) 
      df$data = log10(df$data)
    dfs[["correctedExpression"]] = df
  }
  if (!is.null(scl$strainedExp)) {
    df = DR
    df$correction = "CorrectedCounts"
    if (dataType == "binary") {
      df$data = colSums(scl$atoc[geneSet, rownames(df), 
                                 drop = FALSE]) > 0
    }
    else if (dataType == "counts") {
      df$data = colSums(scl$atoc[geneSet, rownames(df), 
                                 drop = FALSE])
    }
    else if (dataType == "ratio") {
      df$data = (colSums(scl$atoc[geneSet, rownames(df), 
                                  drop = FALSE])/colSums(scl$atoc[, rownames(df), 
                                                                  drop = FALSE]))/(colSums(scl$toc[geneSet, rownames(df), 
                                                                                                   drop = FALSE])/scl$nUMIs[rownames(df)])
    }
    else {
      df$data = colSums(scl$atoc[geneSet, rownames(df), 
                                 drop = FALSE])/colSums(scl$atoc[, rownames(df), 
                                                                 drop = FALSE])
    }
    if (logData) 
      df$data = log10(df$data)
    dfs[["correctedCounts"]] = df
  }
  dfs = do.call(rbind, dfs)
  dfs = dfs[dfs$correction %in% includePanels, ]
  lvls = includePanels
  dfs$correction = factor(dfs$correction, levels = lvls[lvls %in% 
                                                          dfs$correction])
  if (dataType == "ratio") {
    dfs$data[dfs$data < -1] = -1
    dfs$data[dfs$data > 1] = 1
  }
  gg = ggplot(dfs, aes(RD1, RD2)) + geom_point(aes(colour = data), 
                                               size = 0.5) + xlab("ReducedDim1") + ylab("ReducedDim2") + 
    labs(colour = "geneSet") + ggtitle(geneSet) + 
    facet_wrap(~correction)
  if (dataType == "ratio") 
    gg = gg + scale_colour_gradientn(colours = rainbow(50), 
                                     limits = c(-1, 1))
  gg$df = dfs
  gg
}
