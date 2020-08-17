###
library(visNetwork)
library(dplyr)
library(purrr)
library(ggraph)
library(tidygraph)
library(SignallingSingleCell)
library(googleVis)
library(clusterProfiler)
library(org.Mm.eg.db)
source("receptor_ligand/mod_id_rl.R")
set.seed(100)

###############################################################################################
##### Build Network Table
###############################################################################################

epidnorm = readRDS(file = "~/UmassDropbox/randoLab/SCepididymis/data/normcounts_epididymis.rds")

plot_tsne_metadata(epidnorm,color_by = "sub_C_celltype")

epidnorm <- calc_agg_bulk(epidnorm, 
                          aggregate_by = c("sub_C_celltype"), 
                          cutoff_frac = 0.2)

### Use known receptor-ligand database (mouse orthologs from Ramilowski et al., 2015)
load("data/mouse_Receptor_Ligand_Data.Rdata")
epidnorm <- mod_id_rl(epidnorm, Receptor_Ligand_Data)


###############################################################################################
##### Build iGraph Network
###############################################################################################

epidnorm_network <- calc_rl_connections(epidnorm, 
                                        nodes = "sub_C_celltype")
save(epidnorm_network, file = "data/epidnorm_network_sub_C_celltype.Rdata")

# filt_epidnorm_network = filter_rl_network(epidnorm_network, 
#                                            filter_type = "network", 
#                                            filter_by = "log10_Connection_product", 
#                                            cutoff = 2,
#                                            direction = ">")
# 
# 
# # Build the igraph object
# epidnorm_network_igraph <- build_rl_network(input = epidnorm_network,
#                                             #merge_all = T, # when you use a group_by
#                                             value = "log10_Connection_product",
#                                             prefix = "epidnorm_network")
# save(epidnorm_network_igraph, file = "receptor_ligand/epidnorm_igraph_sub_C_celltype.Rdata")
# 
# epidnorm_network_igraph_comparative <- build_rl_network(input = epidnorm_network, 
#                                                                     merge_all = F, 
#                                                                     comparitive = T,
#                                                                     group_by = "genotype",
#                                                                     from = "FX",
#                                                                     to = "WT",
#                                                                     value = "log10_Connection_product", 
#                                                                     prefix = "receptor_ligand/comparative_epidnorm_network")
# 
# #Analyze the main body
# epidnorm_network_C1_analyzed <- analyze_rl_network(epidnorm_network_igraph,
#                                                    subset = 1,
#                                                    subset_on = "connected",
#                                                    prefix = "epidnorm_network_C1_analyzed")
# save(epidnorm_network_C1_analyzed, file = "receptor_ligand/epidnorm_network_C1_analyzed_sub_C_celltype.Rdata")
# 

###############################################################################################
##### Cluster and Label Ligands
###############################################################################################

vals <- which(fData(epidnorm)[,"networks_ligands"])
ligs <- rownames(fData(epidnorm))[vals]
length(ligs)

kres_ligs <- plot_heatmap(epidnorm,
                          genes = ligs,
                          type = "bulk",
                          cluster_by = "row",
                          facet_by = "sub_C_celltype",
                          pdf_format = "tile",
                          scale_by = "row",
                          cluster_type = "kmeans",
                          text_sizes = rep(10,5),
                          k = 13,
                          show_k = T)
ggsave("receptor_ligand/heatmap_filt20_ligs_k13_sub_C_celltype.pdf", h=15,w=15)
saveRDS(kres_ligs, file = "receptor_ligand/kres_filt20_ligs_k13_sub_C_celltype.rds")
kres_ligs = readRDS("receptor_ligand/kres_filt20_ligs_k10_sub_C_celltype.rds")

plot_heatmap(epidnorm,
             genes = names(kres_ligs[[2]]$cluster),
             type = "bulk",
             cluster_by = F,
             facet_by = "sub_C_celltype",
             pdf_format = "tile",
             scale_by = "row",
             text_sizes = c(10,10,5,10,10),
             gene_names = T)
ggsave("receptor_ligand/named_heatmap_filt20_ligs_k13_sub_C_celltype.pdf", h=15,w=15)

lcluster = names(kres_ligs[[2]]$cluster)
ldata = fData(epidnorm)[lcluster,]
ldata$cluster = kres_ligs[[2]]$cluster
write.table(ldata, "receptor_ligand/heatmap_filt20_ligs_k10_sub_C_celltype.tsv", sep = "\t")

###############################################################################################
##### Cluster and Label Receptors
###############################################################################################

vals <- which(fData(epidnorm)[,"networks_Receptors"])
recs <- rownames(fData(epidnorm))[vals]
length(recs)

kres_recs <- plot_heatmap(epidnorm,
                          genes = recs,
                          type = "bulk",
                          cluster_by = "row",
                          facet_by = "sub_C_celltype",
                          pdf_format = "tile",
                          scale_by = "row",
                          cluster_type = "kmeans",
                          text_sizes = rep(10,5),
                          k = 14,
                          show_k = T)
ggsave("receptor_ligand/heatmap_filt20_recs_k13_sub_C_celltype.pdf", h=15,w=15)
saveRDS(kres_recs, file = "receptor_ligand/kres_filt20_recs_k13_sub_C_celltype.rds")
kres_recs = readRDS("receptor_ligand/kres_filt20_recs_k13_sub_C_celltype.rds")

plot_heatmap(epidnorm,
             genes = names(kres_recs[[2]]$cluster),
             type = "bulk",
             cluster_by = F,
             facet_by = "sub_C_celltype",
             pdf_format = "tile",
             scale_by = "row",
             text_sizes = c(10,10,5,10,10),
             gene_names = T)
ggsave("receptor_ligand/named_heatmap_filt20_recs_k13_sub_C_celltype.pdf", h=15,w=15)


rcluster = names(kres_recs[[2]]$cluster)
rdata = fData(epidnorm)[rcluster,]
rdata$cluster = kres_recs[[2]]$cluster
write.table(rdata, "receptor_ligand/heatmap_filt20_recs_k13_sub_C_celltype.tsv", sep = "\t")

###############################################################################################
##### Summary Network
###############################################################################################
innet = epidnorm_network$full_network
colnames(innet)[1] = "celltype_lig"
colnames(innet)[3] = "celltype_rec"
### add ligand cluster identities
lcluster = as.data.frame(kres_ligs[[2]]$cluster)
colnames(lcluster) = "cluster_lig"
lcluster = tibble::rownames_to_column(lcluster, "Ligand")
lcluster$cluster_lig_cell[lcluster$cluster_lig==1] = "Basal"
lcluster$cluster_lig_cell[lcluster$cluster_lig==2] = "Multiple"
lcluster$cluster_lig_cell[lcluster$cluster_lig==3] = "StromalIL6"
lcluster$cluster_lig_cell[lcluster$cluster_lig==4] = "AngryVas"
lcluster$cluster_lig_cell[lcluster$cluster_lig==5] = "ClearNarrow"
lcluster$cluster_lig_cell[lcluster$cluster_lig==6] = "Basal"
lcluster$cluster_lig_cell[lcluster$cluster_lig==7] = "Macs"
lcluster$cluster_lig_cell[lcluster$cluster_lig==8] = "Fibroblast"
lcluster$cluster_lig_cell[lcluster$cluster_lig==9] = "Macs"
lcluster$cluster_lig_cell[lcluster$cluster_lig==10] = "PC07_PC11_PC14"
lcluster$cluster_lig_cell[lcluster$cluster_lig==11] = "PC05_PC15"
lcluster$cluster_lig_cell[lcluster$cluster_lig==12] = "Muscle"
lcluster$cluster_lig_cell[lcluster$cluster_lig==13] = "Endothelial"
innet = merge(lcluster,innet,by="Ligand")
### add receptor cluster identities
rcluster = as.data.frame(kres_recs[[2]]$cluster)
colnames(rcluster) = "cluster_rec"
rcluster = tibble::rownames_to_column(rcluster, "Receptor")
rcluster$cluster_rec_cell[rcluster$cluster_rec==1] = "Basal"
rcluster$cluster_rec_cell[rcluster$cluster_rec==2] = "Clear"
rcluster$cluster_rec_cell[rcluster$cluster_rec==3] = "Macs"
rcluster$cluster_rec_cell[rcluster$cluster_rec==4] = "Tcell"
rcluster$cluster_rec_cell[rcluster$cluster_rec==5] = "Muscle"
rcluster$cluster_rec_cell[rcluster$cluster_rec==6] = "Endothelial"
rcluster$cluster_rec_cell[rcluster$cluster_rec==7] = "Stromal"
rcluster$cluster_rec_cell[rcluster$cluster_rec==8] = "Muscle"
rcluster$cluster_rec_cell[rcluster$cluster_rec==9] = "Basal"
rcluster$cluster_rec_cell[rcluster$cluster_rec==10] = "Dendritic"
rcluster$cluster_rec_cell[rcluster$cluster_rec==11] = "AngryVas"
rcluster$cluster_rec_cell[rcluster$cluster_rec==12] = "ClearNarrow"
rcluster$cluster_rec_cell[rcluster$cluster_rec==13] = "Macs"
rcluster$cluster_rec_cell[rcluster$cluster_rec==14] = "Principal"
innet = merge(rcluster,innet,by="Receptor")
write.table(innet, "receptor_ligand/full_network.tsv", sep = "\t")

##########################################
# Prepare summary input for Sankey plot
##########################################
innet = read.table("receptor_ligand/full_network.tsv", stringsAsFactors = F)
summary <- plyr::count(innet[,c("cluster_lig_cell","cluster_rec_cell")])
summary$cluster_lig_cell <- paste0("lig_", summary$cluster_lig_cell)
summary$cluster_rec_cell <- paste0("rec_", summary$cluster_rec_cell)
summary[,4] <- log10(summary[,3])
summary = summary %>% 
  mutate(connectotal = sum(freq)) %>% 
  group_by(cluster_lig_cell) %>% 
  mutate(ligtotal = sum(freq)) %>% 
  group_by(cluster_rec_cell) %>% 
  mutate(rectotal = sum(freq)) %>%
  rowwise() %>%
  mutate(pval = phyper(freq-1, rectotal, connectotal-rectotal, ligtotal, lower.tail = F)) %>%
  mutate(color = case_when(pval < 0.01 ~ "sig", pval >= 0.01 ~ "ns"))
summary$ligfrac = (summary$freq/summary$ligtotal*100)
summary$recfrac = (summary$freq/summary$rectotal*100)
# sank <- gvisSankey(summary, from = "cluster_lig_cell", to = "cluster_rec_cell", weight = "freq")
# plot(sank)
# write.table(summary, "receptor_ligand/summary_sankey_plot_heatmap_clusters.tsv", sep = "\t")
# summary = read.table("receptor_ligand/summary_sankey_plot_heatmap_clusters.tsv", stringsAsFactors = F)

### Ggplot alluvial style plot with custom colors 
library(ggplot2)
library(ggalluvial)
ggplot(summary,
       aes(y = freq, axis1 = cluster_lig_cell, axis2 = cluster_rec_cell)) +
  geom_alluvium(aes(fill=color), width = 1/12, alpha=0.9) +
  scale_fill_manual(values = c("#e0e0e0","#d73027")) + 
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_text(stat = "stratum", infer.label = TRUE) +
  scale_x_discrete(limits = c("cluster_lig_cell", "cluster_rec_cell"), expand = c(.05, .05)) +
  theme_void()+
  ylab(NULL)+xlab(NULL) +
  ggtitle("Epididymis connections")
ggsave("receptor_ligand/heatmap_clusters_sankey_ggplot.svg",h=10,w=15)