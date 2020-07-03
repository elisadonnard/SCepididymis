library(SoupX)
library(SignallingSingleCell)
## load modified soupX functions:
source('soupX_funcs.R')
###################################################################################

## Read input UMI count files
for (s in c("Caput","Cauda","Corpus","Vasdef")) {
  samp = s
  # inFile = Sys.glob(paste("merged_",samp,".umi.distributions.txt", sep = ""))
  # x = data.table::fread(inFile,header=T,sep='\t')
  # x = as.data.frame(x)
  # rownames(x) = x[,1]
  # x = x[,-1]
  # x = as.matrix(x)
  # ## Save as R data file
  # rDataFile = paste("merged_",samp,".umi.distributions.rds", sep = "")
  # saveRDS(x,file=rDataFile)
  x = readRDS(paste("merged_",samp,".umi.distributions.rds", sep = ""))
  ###################################################################################

  ## Extract just the 'soup' cells and get real cells from the previous analysis
  dim(x)
  message(sprintf('total barcodes: %s',ncol(x)))
  nMin = 350 # minimum UMIs per cell
  nMax = 15000 # maximum UMIs per cell
  cSum = apply(x,2,sum)
  idx.allbc = names(which(cSum>0 & cSum<nMax))    # all barcodes
  idx.keep = names(which(cSum>nMin & cSum<nMax))    # valid cells
  
  ## Soup profile
  soup = x[,!(colnames(x)%in%idx.keep)]
  soupSum = apply(soup,2,sum)
  plot(density(log1p(soupSum)))
  abline(v=log1p(as.numeric(names(which.max(table(soupSum))))))
  
  ## Frequent number of UMIs in empty droplets
  soupcounts = as.numeric(names(which.max(table(soupSum))))
  total = sum(soup)
  upm = as.data.frame(rowSums(soup)/total)*1000000
  colnames(upm) = "upm_soup"
  plot_tsne_gene(epid, "Crisp1")
  upm = merge(upm, fData(epid),by=0)
  upm[,2:length(upm)] = log1p(upm[,2:length(upm)])

  ## Make a SoupChannel object:
  scl = SoupChannel(tod=x[,idx.allbc],
                    toc=x[,idx.keep],
                    channelName=samp, 
                    soupRange = c(0,nMin), 
                    keepDroplets = F)
  saveRDS(scl, paste("ed/batch/soupXfixed/",samp,"_scl_soupchannel.rds", sep = ""))
  scl = readRDS(paste("ed/batch/soupXfixed/",samp,"_scl_soupchannel.rds", sep = ""))
  
  ## Estimate contaminating soup fraction per cell
  scl = calculateContaminationFraction(scl, 
                                       nonExpressedGeneList = markerlist, 
                                       useToEst = cells, 
                                       cellGroups = pData(epid)[colnames(scl$toc),"celltype"])
  plotChannelContamination(scl)
  
  ## Plot cell level contamination fraction calculated by SoupX
  scl = interpolateCellContamination(scl)
  ggplot(as.data.frame(scl$rhos), aes(scl$rhos)) + 
    geom_histogram(binwidth = 0.01) +
    theme_bw()
  ggsave(paste("ed/batch/soupXfixed/",samp,"_rho_distribution.pdf", sep = ""))
  
  ### Re-Estimate contaminating soup fraction per cell using the frequent number of UMIs in empty droplets
  fixedrho = soupcounts/epid$UMI_sum
  names(fixedrho) = colnames(epid)
  epid$fixedrho = fixedrho
  plot_tsne_metadata(epid, color_by = "fixedrho")
  
  ## Add this new fraction to the Soup Channel object and replace the SoupX calculated
  scl$rhos= fixedrho[names(fixedrho)%in%colnames(scl$toc)]
  ggplot(as.data.frame(scl$rhos), aes(scl$rhos)) + 
    geom_histogram(binwidth = 0.01) +
    theme_bw()
  ggsave(paste("ed/batch/soupXfixed/",samp,"_fixedrho_distribution.pdf", sep = ""))
  
  ### Adjust counts to remove likely contaminating RNA from soup
  scl = strainCells(scl)
  scl = adjustCounts(scl, samp)
  rownames(scl$strainedExp) = rownames(scl$toc)
  colnames(scl$strainedExp) = colnames(scl$toc)
  rownames(scl$atoc) = rownames(scl$toc)
  colnames(scl$atoc) = colnames(scl$toc)
  saveRDS(scl, paste("ed/batch/soupXfixed/",samp,"_scl_soupXresult.rds", sep = ""))
  scl = readRDS(paste("ed/batch/soupXfixed/",samp,"_scl_soupXresult.rds", sep = ""))
  
  ### Plot number of genes expressed per cell before and after soupX
  genecount = as.data.frame(apply(scl$toc, 2, function(y) length(y[y>0])))
  colnames(genecount) = "genes"
  genecount$data = "raw"
  genecountX = as.data.frame(apply(as.matrix(scl$atoc), 2, function(y) length(y[y>0])))
  colnames(genecountX) = "genes"
  genecountX$data = "soupX"
  ggplot(rbind(genecount,genecountX), aes(genes, fill = data)) + 
    geom_density(alpha = 0.5) +
    scale_x_log10() +
    theme_bw()
  ggsave(paste("ed/batch/soupXfixed/",samp,"_genecount_distribution.pdf", sep = ""))
  
  ### Plot best corrected gene examples (highest expression in soup for each marker list used)
  pd = pData(epid)[,c("x","y")]
  colnames(pd) = c("RD1","RD2")
  allmarkers = unique(do.call(c, markerlist))
  allmarkerexp = upm$upm_soup[upm$Row.names%in%allmarkers]
  names(allmarkerexp) = upm$Row.names[upm$Row.names%in%allmarkers]
  ### top 20 corrected genes
  examplegenes = names(sort(allmarkerexp, decreasing = T)[1:20])
  pdf(paste("ed/batch/soupXfixed/",samp,"_example_genes.pdf", sep = ""),h=4,w=12)
  for (g in 1:length(examplegenes)) {
    print(plotSoupX(scl = scl, 
                    geneSet = examplegenes[g], 
                    DR = pd, 
                    dataType = "counts"))
  }
  dev.off()
}
