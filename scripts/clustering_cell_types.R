library(SignallingSingleCell)
#########################################################################
### Input data (corrected counts from soupX) and merge
#########################################################################
soupcnt = data.frame()
for (s in c("Caput","Cauda","Corpus","Vasdef")) {
  scl = readRDS(paste("ed/batch/soupXfixed/",s,"_scl_soupXresult.rds", sep = ""))
  cnt = as.matrix(scl$atoc)
  soupcnt = merge(soupcnt,as.data.frame(cnt),by=0, all = T)
  rownames(soupcnt) = soupcnt[,1]
  soupcnt = soupcnt[,-1]
}
dim(soupcnt)  
soupcnt[is.na(soupcnt)] = 0
#########################################################################
### Construct expression set
#########################################################################
soupcnt = construct_ex_sc(soupcnt)
soupcnt = soupcnt[rowSums(exprs(soupcnt))>0,]
soupcnt = soupcnt[,colSums(exprs(soupcnt))>0]
dim(soupcnt)
soupcnt$sample = stringr::str_split_fixed(colnames(soupcnt), pattern = "___", n = 2)[,1]
saveRDS(soupcnt, "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds") 
#########################################################################
### Select variable genes, reduce dimensions and cluster RAW counts
#########################################################################
soupcnt = readRDS("ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
vargenes = subset_genes(soupcnt, method = "CV", cutoff = 0.8, threshold = 0)
soupcnt = dim_reduce(soupcnt, genelist = vargenes, pre_reduce = "vPCA", nComp = 50, nVar = 0.85, scale = F)
soupcnt = cluster_sc(soupcnt, dimension = "2d", method = "density", num_clust = 30, xcol = "x", ycol = "y")
plot_tsne_metadata(soupcnt, color_by = "Cluster")
plot_tsne_metadata(soupcnt, color_by = "sample")
soupcnt$UMI_sum = colSums(exprs(soupcnt))
plot(density(log1p(soupcnt$UMI_sum)))
plot_tsne_metadata(soupcnt, color_by = "UMI_sum")
#########################################################################
### Normalize counts with scran
#########################################################################
genes = rowMeans(exprs(soupcnt))
meancut = quantile(genes, 0.8)
clusterlist = pData(soupcnt)$celltype
sce = SingleCellExperiment::SingleCellExperiment(list(
  counts = exprs(soupcnt)))
sce = scran::computeSumFactors(sce,
                               sizes=c(20, 40, 60, 80, 100),
                               positive=T,
                               clusters=clusterlist,
                               min.mean=meancut)
sf_cells = as.data.frame(sizeFactors(sce))
rownames(sf_cells) = colnames(sce)
colnames(sf_cells) = "sizefactor"
save(sf_cells, file="ed/batch/soupXfixed/corrected_counts/scran_size_factors.RData")
sfs <- sf_cells[sf_cells$sizefactor>0,]
gMean <- exp(mean(log(sfs)))
nLt <- length(sf_cells[sf_cells$sizefactor<(0.1*gMean),])
nGt <- length(sf_cells[sf_cells$sizefactor>(10.0*gMean),])
keepcells = rownames(sf_cells[sf_cells$sizefactor>(0.1*gMean) &
                                sf_cells$sizefactor<(10*gMean),,
                              drop=F])
pd = merge(pData(soupcnt)[,grep("sizefactor",colnames(pData(soupcnt)), invert = T)], as.data.frame(sf_cells), by=0)
rownames(pd) = pd[,1]
pd = pd[,-1]
pData(soupcnt) = pd[colnames(soupcnt),]
saveRDS(soupcnt, file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
norm_counts = sweep(exprs(soupcnt),
                    2,
                    pData(soupcnt)[colnames(soupcnt),"sizefactor"],
                    FUN='/')
norm_counts[!is.finite(norm_counts)] = 0
norm_counts = norm_counts[,keepcells]
norm_counts = as.matrix(norm_counts[rowSums(norm_counts)>0,
                                    colSums(norm_counts)>0])
dim(norm_counts)
normsoup = construct_ex_sc(norm_counts)
pData(normsoup) = pData(soupcnt)[colnames(normsoup),]
normsoup$normUMI_sum = colSums(exprs(normsoup))
plot(density(log1p(normsoup$normUMI_sum)))
saveRDS(normsoup, "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
#########################################################################
### Transform to gene fraction and recluster
#########################################################################
normsoup = readRDS("ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
fgc = apply(exprs(normsoup), 2, function(x) x/sum(x))
exprs(normsoup) = fgc
saveRDS(normsoup, "ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/norm_ex_sc_soupcnt.rds")
vargenes = subset_genes(normsoup, method = "CV", cutoff = 0.85, threshold = 0)
normsoup = dim_reduce(normsoup, genelist = vargenes, pre_reduce = "vPCA", nComp = 50, nVar = 0.95, tSNE_perp = 30, scale = F)
normsoup = cluster_sc(normsoup, dimension = "2d", method = "density", s=1, xcol = "x", ycol = "y")
plot_tsne_metadata(normsoup, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/tSNE_sample.pdf", h=5.5, w=5)
plot_tsne_metadata(normsoup, color_by = "Cluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/tSNE_Cluster.pdf", h=5.5, w=5)
plot_tsne_metadata(normsoup, color_by = "UMI_sum", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/tSNE_UMI_sum.pdf", h=5.5, w=5)
saveRDS(normsoup, "ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/norm_ex_sc_soupcnt.rds") 
#########################################################################
### Find marker genes (DE with edgeR)
#########################################################################
raw = readRDS("ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
normsoup = readRDS("ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/norm_ex_sc_soupcnt.rds")
findDEmarkers(raw,
              pd = pData(normsoup),
              DEgroup="Cluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir="ed/batch/soupXfixed/corrected_counts/results/markers/",
              minCells = 0.01)
#########################################################################
### Redefine cell types
#########################################################################
normsoup = readRDS(file = "ed/batch/soupXfixed/corrected_counts/gene_fraction_postScran/norm_ex_sc_soupcnt.rds")
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
normsoup$celltype[normsoup$Cluster%in%c("Cluster3", "Cluster11", "Cluster6", 
                                        "Cluster2", "Cluster17", "Cluster21",
                                        "Cluster5", "Cluster18", "Cluster12")] = "Principal"
normsoup$celltype[normsoup$Cluster%in%c("Cluster8")] = "Clear"
normsoup$celltype[normsoup$Cluster%in%c("Cluster10", "Cluster19")] = "Muscle"
normsoup$celltype[normsoup$Cluster%in%c("Cluster16", "Cluster14")] = "Immune"
normsoup$celltype[normsoup$Cluster%in%c("Cluster7", "Cluster4", 
                                        "Cluster15", "Cluster13")] = "Stromal"
normsoup$celltype[normsoup$Cluster%in%c("Cluster9", "Cluster20", "Cluster1")] = "Basal"
plot_tsne_metadata(normsoup, color_by = "celltype", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/tSNE_new_celltype.pdf", h=5.5, w=5)
saveRDS(normsoup, "results/norm_ex_sc_soupcnt.rds") 
#########################################################################
### DE markers for new broad cell types
#########################################################################
normsoup = readRDS(file = "results/norm_ex_sc_soupcnt.rds")
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
findDEmarkers(raw,
              pd = pData(normsoup),
              DEgroup="celltype",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/markers/",
              minCells = 0.01 )
#########################################################################
### Aggregate bulk for celltypes and plot heatmap
#########################################################################
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
normsoup = readRDS(file = "results/norm_ex_sc_soupcnt.rds")
raw = raw[,colnames(normsoup)]
pData(raw) = pData(normsoup)[colnames(raw),]
###
all_cluster_markers = scan("results/markers/all_cluster_markers.txt", what = character())
plot_heatmap(raw, all_cluster_markers, "bulk", cluster_by = "both",text_sizes = c(20, 10, 10, 10, 10, 5))
raw = calc_agg_bulk(raw, aggregate_by = c("sample"))
###
raw = calc_agg_bulk(raw, aggregate_by = c("Cluster"))
plot_heatmap(raw, all_cluster_markers, "bulk", cluster_by = "both",text_sizes = c(20, 10, 10, 10, 10, 5))
cellorder = c("Cluster15","Cluster4","Cluster13","Cluster7", #Stromal
              "Cluster10", "Cluster19", # Muscle
              "Cluster9", "Cluster1", "Cluster20", # Basal
              "Cluster14", "Cluster16", # Immune
              "Cluster8", # Clear
              "Cluster12","Cluster5","Cluster18", # Principal
              "Cluster3","Cluster6","Cluster11", # Principal
              "Cluster2","Cluster17","Cluster21") # Vas principal
colnames(fData(raw)) = stringr::str_split_fixed(colnames(fData(raw)), "_", n=2)[,1]
fData(raw) = fData(raw)[,cellorder]
colnames(fData(raw)) = paste0(colnames(fData(raw)),"_bulk")
cluster_result_12k = plot_heatmap(raw, 
                                  all_cluster_markers, 
                                  type = "bulk", 
                                  cluster_by = "row", 
                                  cluster_type = "kmeans", 
                                  text_sizes = c(20, 10, 10, 10, 10, 5),
                                  show_k = T,
                                  k = 13)
saveRDS(cluster_result_12k, "results/markers/cluster_result_12k.rds")
genelist = as.data.frame(cluster_result_12k[[2]]$cluster)
colnames(genelist) = "cluster" 
genelist$gene = rownames(genelist)
write.table(genelist, file = "results/markers/heatmap_clusters_all_cluster_markers.csv", row.names = F, sep = ",")
reordergenes = c(names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==1],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==10],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==13],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==4],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==3],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==12],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==9],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==6],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==5],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==11],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==7],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==8],
                names(cluster_result_12k[[2]]$cluster)[cluster_result_12k[[2]]$cluster==2])
# plot heatmap again but with pre-defined order
plot_heatmap(raw, reordergenes, 
             type = "bulk", cluster_by = F,
             text_sizes = c(20, 10, 10, 10, 10, 5),
             gene_names = F,
             title = "Fig2B")
ggsave("results/heatmap_all_21cluster_markers_12k.pdf", h=8, w=4)
#########################################################################
### Subclustering for each cell type
#########################################################################
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
normsoup = readRDS(file = "results/norm_ex_sc_soupcnt.rds")
raw = raw[,colnames(normsoup)]
pData(raw) = pData(normsoup)[colnames(raw),]
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
pData(normnotfraction) = pData(normsoup)[colnames(normnotfraction),]
#########################################################################
### Subclustering Clear cells
#########################################################################
normsoup_subClear<-subset_ex_sc(normsoup, variable = "celltype", select = c("Clear")) #348 cells
plot_tsne_metadata(normsoup_subClear, color_by = "sample", title = "Clear Cells")
pcgenes = scan("results/markers/markers_Principal_DEmarkers.tsv", what = character())
vargenes = subset_genes(normsoup_subClear, method = "CV", cutoff = 0.99, threshold = 0)
vargenes = setdiff(vargenes,pcgenes)
normsoup_subClear = dim_reduce(normsoup_subClear, genelist = vargenes, 
                               pre_reduce = "vPCA", nComp = 50, nVar = 0.7, 
                               tSNE_perp = 24, scale = F)
normsoup_subClear = cluster_sc(normsoup_subClear, dimension = "2d", 
                               method = "density", num_clust = 3, 
                               xcol = "x", ycol = "y")
plot_tsne_metadata(normsoup_subClear, color_by = "Cluster", size = 2, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Clear/tSNE_clear_Cluster.pdf", h=5.5, w=5)
plot_tsne_metadata(normsoup_subClear, color_by = "sample", size = 2, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Clear/tSNE_clear_sample.pdf", h=5.5, w=5)
normsoup_subClear$mergedCluster = normsoup_subClear$Cluster
saveRDS(normsoup_subClear, "results/Clear/norm_ex_sc_soupcnt_subClear.rds") 
normsoup_subClear = readRDS(file = "results/Clear/norm_ex_sc_soupcnt_subClear.rds")

findDEmarkers(raw,
              pd = pData(normsoup_subClear),
              DEgroup="Cluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Clear/markers/",
              minCells = 0.05)

### plot marker genes with norm values not fraction
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
normnotfraction_Clear = normnotfraction[,colnames(normsoup_subClear)]
pData(normnotfraction_Clear) = pData(normsoup_subClear)[colnames(normnotfraction_Clear),]
plot_tsne_gene(normnotfraction_Clear,"Iapp", log_scale = T)
plot_tsne_gene(normnotfraction_Clear,"Cp", log_scale = T)
plot_tsne_gene(normnotfraction_Clear,"Wfdc13", log_scale = T)


#AggBulk Clear and heatmap
raw_subClear = raw[,colnames(normsoup_subClear)]
pData(raw_subClear) = pData(normsoup_subClear)[colnames(raw_subClear),]
raw_subClear = calc_agg_bulk(raw_subClear, aggregate_by = c("Cluster"))
all_cluster_markers = scan("results/Clear/markers/clear_markers", what = character())
clearheatmap = plot_heatmap(raw_subClear, all_cluster_markers, "bulk", 
                            cluster_by = "row",
                            cluster_type = "kmeans",
                            k=4,
                            show_k = T,
                            text_sizes = c(20, 10, 10, 10, 10, 5))

ggsave("results/Clear/heatmap_Clear_bulk.pdf")
saveRDS(clearheatmap, file = "results/Clear/clearheatmap.rds")
clearheatmap = readRDS("results/Clear/clearheatmap.rds")
names(clearheatmap[[2]]$cluster[clearheatmap[[2]]$cluster=="2"])
geneorder = names(clearheatmap[[2]]$cluster)
plot_heatmap(normsoup_subClear, 
             geneorder, 
             "single_cell", 
             cluster_by = "col",
             facet_by = "Cluster",
             ceiling = 2,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Clear/heatmap_Clear_single_cell.pdf", h=30, w=20)


#########################################################################
### Subclustering Stromal cells
#########################################################################

normsoup_subStromal=subset_ex_sc(normsoup, variable = "celltype", select = c("Stromal")) #2222 cells
plot_tsne_metadata(normsoup_subStromal, color_by = "sample", title = "Stromal Cells")
pcgenes = scan("results/markers/markers_Principal_DEmarkers.tsv", what = character())
vargenes = subset_genes(normsoup_subStromal, method = "CV", cutoff = 0.97, threshold = 0)
vargenes = setdiff(vargenes,pcgenes)
normsoup_subStromal = dim_reduce(normsoup_subStromal, genelist = vargenes, 
                                 pre_reduce = "vPCA", nComp = 50, nVar = 0.7, 
                                 tSNE_perp = 100, scale = F)
normsoup_subStromal = cluster_sc(normsoup_subStromal, 
                                 dimension = "2d", 
                                 method = "density", 
                                 num_clust = 6,
                                 xcol = "x", ycol = "y")
plot_tsne_metadata(normsoup_subStromal, color_by = "Cluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Stromal/tSNE_Stromal_Cluster.pdf", h=5.5, w=5)
plot_tsne_metadata(normsoup_subStromal, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Stromal/tSNE_Stromal_sample.pdf", h=5.5, w=5)

normsoup_subStromal$mergedCluster = normsoup_subStromal$Cluster
normsoup_subStromal$mergedCluster[normsoup_subStromal$Cluster=="Cluster6"] = "Cluster1"
plot_tsne_metadata(normsoup_subStromal, color_by = "mergedCluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Stromal/tSNE_Stromal_mergedCluster.pdf", h=5.5, w=5)

saveRDS(normsoup_subStromal, "results/Stromal/norm_ex_sc_soupcnt_subStromal.rds") 
normsoup_subStromal = readRDS(file = "results/Stromal/norm_ex_sc_soupcnt_subStromal.rds")

findDEmarkers(raw,
              pd = pData(normsoup_subStromal),
              DEgroup="Cluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Stromal/markers/",
              minCells = 0.01 )

findDEmarkers(raw,
              pd = pData(normsoup_subStromal),
              DEgroup="mergedCluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Stromal/markers/merged/",
              minCells = 0.01 )

### plot marker genes with norm values not fraction
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
normnotfraction_Stromal = normnotfraction[,colnames(normsoup_subStromal)]
pData(normnotfraction_Stromal) = pData(normsoup_subStromal)[colnames(normnotfraction_Stromal),]
plot_tsne_gene(normnotfraction_Stromal,"Cxcl2", log_scale = T)
plot_tsne_gene(normnotfraction_Stromal,"Ccl2", log_scale = T)
plot_tsne_gene(normnotfraction_Stromal,"Tspan8", log_scale = T)

#AggBulk Stromal and heatmap
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
raw_subStromal = raw[,colnames(normsoup_subStromal)]
pData(raw_subStromal) = pData(normsoup_subStromal)[colnames(raw_subStromal),]
raw_subStromal = calc_agg_bulk(raw_subStromal, aggregate_by = c("Cluster"))
###
stromalheatmap = plot_heatmap(raw_subStromal, 
                              vargenes, 
                              "bulk", 
                              cluster_by = "row",
                              cluster_type = "kmeans",
                              k=7,
                              show_k = T,
                              text_sizes = c(20, 10, 10, 10, 10, 5))

ordergenes = names(stromalheatmap[[2]]$cluster)
### check single cell heatmap for clusters
plot_heatmap(normsoup_subStromal, 
             ordergenes, 
             "single_cell", 
             cluster_by = "col",
             facet_by = "Cluster",
             ceiling = 2,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))

### recluster C2+C4+C5
stromalC2C4C5 = normsoup_subStromal[,normsoup_subStromal$Cluster%in%c("Cluster2","Cluster4","Cluster5")]
vargenes = subset_genes(stromalC2C4C5, method = "CV", cutoff = 0.97, threshold = 0)
vargenes = setdiff(vargenes,pcgenes)
stromalC2C4C5 = dim_reduce(stromalC2C4C5, genelist = vargenes, 
                                 pre_reduce = "vPCA", nComp = 50, nVar = 0.65, 
                                 tSNE_perp = 70, scale = F)
stromalC2C4C5 = cluster_sc(stromalC2C4C5, 
                                 dimension = "2d", 
                                 method = "density", 
                                 num_clust = 6,
                                 xcol = "x", ycol = "y")
plot_tsne_metadata(stromalC2C4C5, color_by = "Cluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot_heatmap(stromalC2C4C5, 
             ordergenes, 
             "single_cell", 
             cluster_by = "col",
             facet_by = "Cluster",
             ceiling = 2,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))
saveRDS(stromalC2C4C5, "results/Stromal/norm_ex_sc_soupcnt_stromalC2C4C5.rds") 

normsoup_subStromal$Cluster[colnames(normsoup_subStromal)%in%colnames(stromalC2C4C5[,stromalC2C4C5$Cluster=="Cluster1"])] = "Cluster2"
normsoup_subStromal$Cluster[colnames(normsoup_subStromal)%in%colnames(stromalC2C4C5[,stromalC2C4C5$Cluster=="Cluster2"])] = "Cluster4"
normsoup_subStromal$Cluster[colnames(normsoup_subStromal)%in%colnames(stromalC2C4C5[,stromalC2C4C5$Cluster%in%c("Cluster3","Cluster4","Cluster5","Cluster6")])] = "Cluster5"
plot_tsne_metadata(normsoup_subStromal, color_by = "Cluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Stromal/tSNE_Stromal_newClusters.pdf", h=5.5, w=5)
saveRDS(normsoup_subStromal, "results/Stromal/norm_ex_sc_soupcnt_subStromal.rds") 

###
all_cluster_markers = scan("results/Stromal/markers/stromal_markers", what = character())
stromalheatmap = plot_heatmap(raw_subStromal, 
                              all_cluster_markers, 
                              "bulk", 
                              cluster_by = "row",
                              cluster_type = "kmeans",
                              k = 7,
                              show_k = T,
                              text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Stromal/heatmap_Stromal_markers_7k.pdf", h=8, w=4)
saveRDS(stromalheatmap, "results/Stromal/stromalheatmap.rds")
stromalheatmap = readRDS("results/Stromal/stromalheatmap.rds")
ordergenes = names(stromalheatmap[[2]]$cluster)

plot_heatmap(normsoup_subStromal, 
             ordergenes[ordergenes%in%rownames(normsoup_subStromal)], 
             "single_cell", 
             cluster_by = "col",
             facet_by = "Cluster",
             ceiling = 2,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Stromal/single_cell_heatmap_Stromal_markers_7k.pdf", h=8, w=4)


#########################################################################
### Subclustering Muscle cells
#########################################################################

normsoup_subMuscle<-subset_ex_sc(normsoup, variable = "celltype", select = c("Muscle")) #502 cells
plot_tsne_metadata(normsoup_subMuscle, color_by = "sample", title = "Muscle Cells")#C10+C19

pcgenes = scan("results/markers/markers_Principal_DEmarkers.tsv", what = character())
normnotfraction_Muscle = normnotfraction[,colnames(normsoup_subMuscle)]
vargenes = subset_genes(normnotfraction_Muscle, method = "CV", cutoff = 0.9, threshold = 1)
vargenes = setdiff(vargenes,pcgenes)

normsoup_subMuscle = dim_reduce(normsoup_subMuscle, genelist = vargenes, 
                                 pre_reduce = "vPCA", nComp = 50, nVar = 0.8, 
                                 tSNE_perp = 64, scale = F)
normsoup_subMuscle = cluster_sc(normsoup_subMuscle, 
                                 dimension = "2d", 
                                 method = "density", 
                                 num_clust = 7,
                                 xcol = "x", ycol = "y")
plot_tsne_metadata(normsoup_subMuscle, color_by = "Cluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Muscle/tSNE_Muscle_Cluster.pdf", h=5.5, w=5)
plot_tsne_metadata(normsoup_subMuscle, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Muscle/tSNE_Muscle_sample.pdf", h=5.5, w=5)

# merge clusters
normsoup_subMuscle$mergedCluster = "Cluster1"
normsoup_subMuscle$mergedCluster[normsoup_subMuscle$Cluster=="Cluster1"] = "Cluster2"
normsoup_subMuscle$mergedCluster[normsoup_subMuscle$Cluster=="Cluster6"] = "Cluster2"
normsoup_subMuscle$mergedCluster[normsoup_subMuscle$Cluster=="Cluster3"] = "Cluster2"
normsoup_subMuscle$mergedCluster[normsoup_subMuscle$Cluster=="Cluster7"] = "Cluster2"
plot_tsne_metadata(normsoup_subMuscle, color_by = "mergedCluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Muscle/tSNE_Muscle_mergedCluster.pdf", h=5.5, w=5)

saveRDS(normsoup_subMuscle, "results/Muscle/norm_ex_sc_soupcnt_subMuscle.rds") 
normsoup_subMuscle = readRDS(file = "results/Muscle/norm_ex_sc_soupcnt_subMuscle.rds")

findDEmarkers(raw,
              pd = pData(normsoup_subMuscle),
              DEgroup="mergedCluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Muscle/markers/",
              minCells = 0.01 )

### plot marker genes with norm values not fraction
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
normnotfraction_Muscle = normnotfraction[,colnames(normsoup_subMuscle)]
pData(normnotfraction_Muscle) = pData(normsoup_subMuscle)[colnames(normnotfraction_Muscle),]
plot_tsne_gene(normnotfraction_Muscle, "Atp6ap2", log_scale = T)
plot_tsne_gene(normnotfraction_Muscle, "Cox17", log_scale = T)

#AggBulk Muscle
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
raw_subMuscle = raw[,colnames(normsoup_subMuscle)]
pData(raw_subMuscle) = pData(normsoup_subMuscle)[colnames(raw_subMuscle),]
raw_subMuscle = calc_agg_bulk(raw_subMuscle, aggregate_by = c("Cluster"))
###

all_cluster_markers = scan ("results/Muscle/markers/muscle_markers", what = character())
muscleheatmap = plot_heatmap(raw_subMuscle, 
                             all_cluster_markers, 
                             "bulk", 
                             cluster_by = "row",
                             cluster_type = "kmeans",
                             k=2,
                             show_k = T,
                             text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Muscle/heatmap_muscle_markers.pdf", h=5, w=5)
saveRDS(muscleheatmap, file= "results/Muscle/muscleheatmap.rds")
geneorder = names(muscleheatmap[[2]]$cluster)
plot_heatmap(normnotfraction_Muscle, 
             geneorder, 
             "single_cell", 
             cluster_by = "col",
             facet_by = "mergedCluster",
             ceiling = 1.5,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Muscle/single_cell_HeatMapMusclesubCluster.pdf", h=5, w=10)

#########################################################################
### Subclustering Principal cells
#########################################################################
epidfraction_subPrincipal = subset_ex_sc(epidfraction, variable = "celltype", select = c("Principal"))
vargenes = subset_genes(epidfraction_subPrincipal, method = "CV", cutoff = 0.90, threshold = 0)
epidfraction_subPrincipal = dim_reduce(epidfraction_subPrincipal, genelist = vargenes, 
                                       pre_reduce = "vPCA", nComp = 30, nVar = 0.93, 
                                       tSNE_perp = 65, scale = F)
epidfraction_subPrincipal = cluster_sc(epidfraction_subPrincipal, 
                                       dimension = "2d", 
                                       method = "density", 
                                       num_clust = 15,
                                       xcol = "x", ycol = "y")
plot_tsne_metadata(epidfraction_subPrincipal, color_by = "Cluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot_tsne_metadata(epidfraction_subPrincipal, color_by = "UMI_sum", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot_tsne_metadata(epidfraction_subPrincipal, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))

plot_tsne_metadata(epidfraction_subPrincipal, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
saveRDS(epidfraction_subPrincipal,"results/Principal/epidfraction_subPrincipal.rds") 
raw = readRDS("results/rawcounts_epididymis.rds")
epidfraction_subPrincipal = readRDS("results/Principal/epidfraction_subPrincipal.rds")
epidnorm_subPrincipal = readRDS("results/Principal/epidnorm_subPrincipal.rds")

findDEmarkers(raw,
              pd = pData(epidfraction_subPrincipal),
              DEgroup="Cluster",
              sizefactor = "sizefactor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Principal/markers/Cluster_Fraction/",
              minCells = 0.01 )


plot_tsne_metadata(epidfraction_subPrincipal, 
                   color_by = "Cluster", 
                   size = 0.8, 
                   text_sizes = c(20, 10, 5, 10, 10, 5),
                   facet_by = "Cluster") 


#########################################################################
### Subclustering Basal cells
#########################################################################

normsoup_subBasal<-subset_ex_sc(normsoup, variable = "celltype", select = c("Basal")) #1412 cells
plot_tsne_metadata(normsoup_subBasal, color_by = "sample", title = "Basal Cells")
plot_tsne_metadata(normsoup_subBasal, color_by = "UMI_sum", title = "Basal Cells")

pcgenes = scan("results/markers/markers_Principal_DEmarkers.tsv", what = character())
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
pData(normnotfraction) = pData(normsoup)[colnames(normnotfraction),]
normnotfraction_Basal = normnotfraction[,colnames(normsoup_subBasal)]
vargenes = subset_genes(normnotfraction_Basal, method = "CV", cutoff = 0.95, threshold = 1)
vargenes = setdiff(vargenes,pcgenes)
normsoup_subBasal = dim_reduce(normsoup_subBasal, genelist = vargenes, 
                                   pre_reduce = "vPCA", nComp = 50, nVar = 0.9, 
                                   tSNE_perp = 90, scale = F)

normsoup_subBasal = cluster_sc(normsoup_subBasal, dimension = "2d", 
                              method = "density", num_clust = 11, 
                               xcol = "x", ycol = "y")

plot_tsne_metadata(normsoup_subBasal, color_by = "Cluster", title = "Basal", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Basal/tSNE_Basal_Cluster.eps", h=5.5, w=5)
plot_tsne_metadata(normsoup_subBasal, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Basal/tSNE_Basal_sample.pdf", h=5.5, w=5)

# merge a few clusters back
normsoup_subBasal$mergedCluster[normsoup_subBasal$Cluster%in%c("Cluster2","Cluster4","Cluster7")] = "Cluster1"
normsoup_subBasal$mergedCluster[normsoup_subBasal$Cluster%in%c("Cluster3","Cluster5","Cluster6","Cluster8", "Cluster9", "Cluster10")] = "Cluster2"
normsoup_subBasal$mergedCluster[normsoup_subBasal$Cluster%in%c("Cluster1")] = "Cluster3"
normsoup_subBasal$mergedCluster[normsoup_subBasal$Cluster%in%c("Cluster11")] = "Cluster4"
plot_tsne_metadata(normsoup_subBasal, color_by = "mergedCluster", title = "Basal", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Basal/tSNE_Basal_mergedCluster.eps", h=5.5, w=5)


saveRDS(normsoup_subBasal, "results/Basal/norm_ex_sc_soupcnt_subBasal.rds") 
normsoup_subBasal = readRDS(file = "results/Basal/norm_ex_sc_soupcnt_subBasal.rds")

findDEmarkers(raw,
              pd = pData(normsoup_subBasal),
              DEgroup="mergedCluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Basal/markers/",
              minCells = 0.01 )

### plot marker genes with norm values not fraction
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
normnotfraction_Basal = normnotfraction[,colnames(normsoup_subBasal)]
pData(normnotfraction_Basal) = pData(normsoup_subBasal)[colnames(normnotfraction_Basal),]
plot_tsne_metadata(normnotfraction_Basal, color_by = "mergedCluster", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot_tsne_gene(normnotfraction_Basal,"Cldn4") # Cluster 1
plot_tsne_gene(normnotfraction_Basal,"Ctsl") # Cluster 1
plot_tsne_gene(normnotfraction_Basal,"Cd9") # Cluster 1
plot_tsne_gene(normnotfraction_Basal,"Actb") 
plot_tsne_gene(normnotfraction_Basal,"Slc39a1") 
plot_tsne_gene(normnotfraction_Basal,"Lrrc8a")
plot_tsne_gene(normnotfraction_Basal,"Cd83")
plot_tsne_gene(normnotfraction_Basal,"Cd14")
plot_tsne_metadata(normsoup_subBasal, color_by = "UMI_sum", title = "Basal", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))

#AggBulk Basal
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
raw_subBasal = raw[,colnames(normsoup_subBasal)]
pData(raw_subBasal) = pData(normsoup_subBasal)[colnames(raw_subBasal),]
raw_subBasal = calc_agg_bulk(raw_subBasal, aggregate_by = c("mergedCluster"))

###
all_cluster_markers = scan("results/Basal/markers/all_cluster_markers.txt", what = character())
basalheatmap = plot_heatmap(raw_subBasal, all_cluster_markers, "bulk", 
                            cluster_by = "row",
                            cluster_type = "kmeans",
                            k=4,
                            show_k = T,
                            text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Basal/HeatMapBasalsubCluster.eps", h=5.5, w=5)
geneorder = names(basalheatmap[[2]]$cluster)
plot_heatmap(normnotfraction_Basal, 
             geneorder, 
             "single_cell", 
             cluster_by = "col",
             facet_by = "mergedCluster",
             ceiling = 1.5,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Basal/single_cell_HeatMapBasalsubCluster.pdf", h=5, w=10)

#########################################################################
### Subclustering Immune cells
#########################################################################

normsoup_subImmune = subset_ex_sc(normsoup, variable = "celltype", select = c("Immune")) #440 cells
plot_tsne_metadata(normsoup_subImmune, color_by = "sample", title = "Immune Cells")
pcgenes = scan("results/markers/markers_Principal_DEmarkers.tsv", what = character())
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
pData(normnotfraction) = pData(normsoup)[colnames(normnotfraction),]
normnotfraction_Immune = normnotfraction[,colnames(normsoup_subImmune)]
vargenes = subset_genes(normnotfraction_Immune, method = "CV", cutoff = 0.9, threshold = 1)
vargenes = setdiff(vargenes,pcgenes)
normsoup_subImmune = dim_reduce(normsoup_subImmune, genelist = vargenes, 
                               pre_reduce = "vPCA", nComp = 50, nVar = 0.85, 
                               tSNE_perp = 20, scale = F)
normsoup_subImmune = cluster_sc(normsoup_subImmune, 
                                    dimension = "2d", 
                                    method = "density", 
                                    num_clust = 7,
                                    xcol = "x", ycol = "y")

plot_tsne_metadata(normsoup_subImmune, color_by = "Cluster", title = "Immune", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Immune/tSNE_Immune_Cluster.pdf", h=5.5, w=5)
plot_tsne_metadata(normsoup_subImmune, color_by = "sample", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Immune/tSNE_Immune_sample.pdf", h=5.5, w=5)

# merge clusters
normsoup_subImmune$mergedCluster = normsoup_subImmune$Cluster
normsoup_subImmune$mergedCluster[normsoup_subImmune$Cluster=="Cluster4"] = "Cluster1" # endothelial
normsoup_subImmune$mergedCluster[normsoup_subImmune$Cluster=="Cluster7"] = "Cluster1" # endothelial
normsoup_subImmune$mergedCluster[normsoup_subImmune$Cluster=="Cluster6"] = "Cluster4"
plot_tsne_metadata(normsoup_subImmune, 
                   color_by = "mergedCluster", 
                   title = "Immune", 
                   size = 0.8, 
                   text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("results/Immune/tSNE_Immune_mergedCluster.pdf", h=4.5, w=4)
ggplot(pData(normsoup_subImmune), aes(x,y, colour=mergedCluster)) +
  geom_point(size=1) +
  theme_void() +
  theme(legend.position = "none")
ggsave("results/Immune/void_tSNE_Immune_mergedCluster.pdf", h=4, w=4)

saveRDS(normsoup_subImmune, "results/Immune/norm_ex_sc_soupcnt_subImmune.rds") 
normsoup_subImmune = readRDS(file = "results/Immune/norm_ex_sc_soupcnt_subImmune.rds")
table(normsoup_subImmune$mergedCluster,normsoup_subImmune$sample)

findDEmarkers(raw,
              pd = pData(normsoup_subImmune),
              DEgroup="Cluster",
              sizefactor = "size_factor",
              lib_size = "UMI_sum",
              batchID="sample",
              outdir= "results/Immune/markers/",
              minCells = 0.01 )

### plot marker genes with norm values not fraction
normnotfraction = readRDS(file = "ed/batch/soupXfixed/corrected_counts/norm_ex_sc_soupcnt.rds")
pData(normnotfraction) = pData(normsoup)[colnames(normnotfraction),]
normnotfraction_Immune = normnotfraction[,colnames(normsoup_subImmune)]
pData(normnotfraction_Immune) = pData(normsoup_subImmune)[colnames(normnotfraction_Immune),]
pdf("results/Immune/immune_cluster_marker_genes.pdf")
plot_tsne_gene(normnotfraction_Immune, "Cxcl10", log_scale = T, title = "Macrophages") # macrophages (cluster 5)
plot_tsne_gene(normnotfraction_Immune, "Ccr7", log_scale = T, title = "Dendritic cells") # dendritic cells (cluster 4)
plot_tsne_gene(normnotfraction_Immune, "Cd3g", log_scale = T, title = "T cells (not CD4 or CD8, maybe NK)") # T cells (not CD4 or CD8, maybe NK) (cluster 2)
plot_tsne_gene(normnotfraction_Immune, "Tm4sf1", log_scale = T, title = "Endothelial") # endothelial cells (cluster 1)
plot_tsne_gene(normnotfraction_Immune, "Trem2", log_scale = T, title = "Trem2+ Macrophages") # Trem2+ macrophages (cluster 1)
dev.off()

#AggBulk Immune
raw = readRDS(file = "ed/batch/soupXfixed/corrected_counts/merged_ex_sc_soupcnt.rds")
raw_subImmune = raw[,colnames(normsoup_subImmune)]
pData(raw_subImmune) = pData(normsoup_subImmune)[colnames(raw_subImmune),]
raw_subImmune = calc_agg_bulk(raw_subImmune, aggregate_by = c("Cluster"))

###
all_cluster_markers = scan ("results/Immune/markers/markers_immune", what = character())
immuneheatmap = plot_heatmap(raw_subImmune, 
                             all_cluster_markers, 
                             "bulk", 
                             cluster_by = "row",
                             cluster_type = "kmeans",
                             k=5,
                             show_k = T,
                             text_sizes = c(20, 10, 10, 10, 10, 5))
saveRDS(immuneheatmap, file = "results/Immune/immuneheatmap.rds")
#Getting the genes from the heatmap
genelist = immuneheatmap[[2]]$cluster
reordergenes = c(rev(names(immuneheatmap[[2]]$cluster)[immuneheatmap[[2]]$cluster==4]),
                 names(immuneheatmap[[2]]$cluster)[immuneheatmap[[2]]$cluster==3],
                 names(immuneheatmap[[2]]$cluster)[immuneheatmap[[2]]$cluster==1],
                 names(immuneheatmap[[2]]$cluster)[immuneheatmap[[2]]$cluster==2],
                 names(immuneheatmap[[2]]$cluster)[immuneheatmap[[2]]$cluster==5])
# plot heatmap again but with pre-defined order
plot_heatmap(raw_subImmune, 
             reordergenes, 
             type = "bulk", 
             cluster_by = "col",
             text_sizes = c(20, 10, 10, 10, 10, 5),
             gene_names = F,
             title = "Immune")
ggsave("results/Immune/heatmap_Immune_markers.eps", h=8, w=4)

plot_heatmap(raw_subImmune, 
             reordergenes, 
             "single_cell", 
             cluster_by = "col",
             facet_by = "Cluster",
             ceiling = 2,
             group_names = F,
             text_sizes = c(20, 10, 10, 10, 10, 5))
ggsave("results/Immune/single_cell_HeatMapImmunesubCluster.pdf", h=5, w=10)

#########################################################################
### Adding Mitochondrial UMI counts to dataset
#########################################################################
epidnorm = readRDS(file = "results/normcounts_epididymis.rds")
mt = read.table("results/all_chrM_gene_umis", stringsAsFactors = F, header = F)
rownames(mt) = mt$V2
colnames(mt) = c("gene_mtUMIs","cells")
pd = merge(pData(epidnorm),mt,by=0, all.x = T)
pd = pd[,grep("cells", colnames(pd), invert = T)]
rownames(pd) = pd[,1]
pd = pd[,-1]
pd[is.na(pd)] = 0
pData(epidnorm) = pd[colnames(epidnorm),]
epidnorm$log_mtUMIs = log1p(epidnorm$gene_mtUMIs)
plot_tsne_metadata(epidnorm, color_by = "gene_mtUMIs", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot_tsne_metadata(epidnorm, color_by = "log_mtUMIs", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot(density(epidnorm$log_mtUMIs))
epidnorm$mtUMIs_norm = epidnorm$gene_mtUMIs/(epidnorm$UMI_sum+epidnorm$gene_mtUMIs)
plot_tsne_metadata(epidnorm, color_by = "mtUMIs_norm", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
ggsave("ed/plots/mtUMIS_tSNE.png")
epidnorm$log_mtUMIs_norm = log1p(epidnorm$mtUMIs_norm)
plot_tsne_metadata(epidnorm, color_by = "log_mtUMIs_norm", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
plot_tsne_metadata(epidnorm, color_by = "UMI_sum", size = 0.8, text_sizes = c(20, 10, 5, 10, 10, 5))
quantile(epidnorm$mtUMIs_norm,seq(0,1,0.1))
plot(density(epidnorm$mtUMIs_norm))
ggplot(pData(epidnorm), aes(mtUMIs_norm)) +
  geom_histogram(binwidth = 0.01) +
  #facet_wrap(~Cluster) +
  geom_vline(xintercept = 0.3)
table(epidnorm$Cluster)
ggsave("ed/plots/mtUMIS_distribution.png")
saveRDS(epidnorm,"results/normcounts_epididymis.rds")
