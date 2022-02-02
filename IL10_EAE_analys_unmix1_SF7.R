##BY: Patrick O'Connell
## 05/20/2021


#This analysis starts from FCS files and goes from there. W/ downsampling. 

#Packages
library(CATALYST)
library(tidyverse)
library(flowCore)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(diffcyt)

setwd("~/IL10_EAE_R_analys_unmix1_SF7")

##-------------------------Making metadata file----------------------------------##
#Load the content of the metadata.xlsx file into R (file name must match name of the fcs file you generate for each sample)
setwd("~/IL10_EAE_R_analys_unmix1_SF7")
metadata_filename <- "eae_metadata_altX.xlsx"
md <- read_excel(metadata_filename)

# Define condition variables as named in metadata
md$condition <- factor(md$condition, levels = c("EAE", "naive"))
md$IV_pos <- factor(md$IV_pos, levels = c("yes", "no"))
md$IL_10_pos <- factor(md$IL_10_pos, levels = c("yes", "no"))
md$genotype <- factor(md$genotype, levels = c("ERAP1_KO", "WT", "SLAMF7_KO"))
head(data.frame(md))
##-----------------------------------------------------------##

#Read FCS files in as flowset 
setwd("~/IL10_EAE_R_analys_unmix1_SF7/clean_il10_EAE_fcs")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

##------------------making panel file-----------------------##
setwd("~/IL10_EAE_R_analys_unmix1_SF7")
panel <- "neuroimmune_panel_IL10.xlsx"
panel <- read_excel(panel) 

#check
all(panel$fcs_colname %in% colnames(fcs_raw))
##------------------making panel file-----------------------##

#see how many cells per file
fsApply(fcs_raw, nrow)

#downsample to 100,000 cells per sample
# Define a downsampling ceiling
sampling.ceiling <- 100000
# Being reproducible is a plus
set.seed(666)

# BUILD A DOWNSAMPLED FLOWSET
fcs_raw.dsamp <- fsApply(fcs_raw, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling.ceiling, nrow(ff)))
  ff[idx,]  
})

#check it
fcs_raw.dsamp
fsApply(fcs_raw.dsamp, nrow)

#make Catalyst data object and arcsinh t-for all channels by 6000

#First 3 must be present and sample_id must match condition (i think)
factors <- list(factors = c("file_name", "condition", "sample_id", "genotype", "IV_pos", "IL_10_pos"))

daf_X <- prepData(
  fcs_raw.dsamp,
  features = NULL,
  md = md,
  md_cols = factors,
  panel = panel,
  transform = TRUE,
  cofactor = 6000,
  by_time = FALSE,
  FACS = TRUE
)

#remove erap1 samples
daf_X_sf7_1 <- filterSCE(daf_X, genotype != "ERAP1_KO")

#remove naive samples since we only have EAE for SF7 genotype
daf_X_sf7 <- filterSCE(daf_X_sf7_1, condition == "EAE")

levels(daf_X_sf7$condition)
levels(daf_X_sf7$genotype)


#Make MDS plot  
pbMDS(daf_X_sf7, by = "sample_id", color_by = "genotype")

#Heatmap of all markers per sample_id and condition w/ heiracheal clustering.
plotExprHeatmap(daf_alt, 
                bin_anno = TRUE, 
                scale = "first",
                row_anno = "condition"
                )

#Perform FlowSOM clustering
daf_X_sf7 <- cluster(daf_X_sf7, 
               features = "type", 
               xdim = 10, 
               ydim = 10, 
               maxK = 22, 
               seed = 1234, 
               verbose = TRUE
               ) 

#plot cluster heatmap
plotClusterHeatmap(daf_X_sf7,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta22", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
                   ) 

#plot marker expression by cluster
plotClusterExprs(daf_X_sf7, k = "meta22")  

# run UMAP                          
set.seed(777)                               
daf_X_sf7 <- runDR(daf_X_sf7, 
             dr = "UMAP", 
             cells = 10000,
             features = "type" 
             )

#save main analysis object
saveRDS(daf_X_sf7, file = "./daf_X_sf7.rds")

#UMAP visualizations
plotDR(daf_X_sf7, "UMAP", color_by = "meta22")

pz <- plotDR(daf_X_sf7, "UMAP", color_by = "meta22")               
pz + theme(axis.line = element_line(colour = NA),   #use this theme for all plots
    axis.ticks = element_line(colour = NA), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(colour = NA), 
    axis.text = element_text(colour = NA), 
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 14), 
    panel.background = element_rect(fill = NA)) 

## Facet per sample                                          
plotDR(daf, "UMAP", color_by = "meta20", facet = "sample_id")

## Facet per condition                                       
plotDR(daf_alt2, "UMAP", color_by = "meta20", facet = "IV_pos")

#plot FlowSOM codes to visualize similarity of clusters
plotCodes(daf_X, k = "meta22")

#Merge and annotate clusters
merging_tableX <- "merged_clusters_X.xlsx"
merging_table1 <- read_excel(merging_tableX)
head(data.frame(merging_table1))  

# convert to factor with merged clusters in desired order                
merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("Microglia", "BAMs", "MdCs", "pDCs", "DCs", "CCR2+ monocytes", "Neutrophils",             
                                                "CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Debris"))        

daf_X_sf7 <- mergeClusters(daf_X_sf7, k = "meta22",                                  
                     table = merging_table1, 
                     id = "merging1",
                     overwrite = TRUE
                     )   

#remove debris cluster
#make a new SCE object when you do this!!
daf_X_sf7 <- filterSCE(daf_X_sf7, cluster_id != "Debris", k = "merging1")

#change cluster colors
clust_color <- c("#996600", "#00CCAA", "#6698FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC")

daf_X_sf7 <- readRDS("./daf_X_SF7.rds")
#UMAP of merged clusters
plotDR(daf_X_sf7, "UMAP", color_by = "merging1", k_pal = clust_color) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
                                                                axis.ticks = element_line(colour = NA), 
                                                                panel.grid.major = element_line(linetype = "blank"), 
                                                                panel.grid.minor = element_line(linetype = "blank"), 
                                                                axis.title = element_text(colour = NA), 
                                                                axis.text = element_text(colour = NA), 
                                                                legend.text = element_text(size = 10), 
                                                                legend.title = element_text(size = 14), 
                                                                panel.background = element_rect(fill = NA)) 

#heatmap of merged cluster markers 
plotExprHeatmap(daf_X_sf7, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color, 
                features = c("CD45","CD11b","Tim3", "CD4", "CD8", "NK1.1","CD19","NKG2D","Ly6C","Ly6G", "MHCII","SiglecH", "CD11c","B220","CD49b","IgD", "CCR2","CD38", "CD206","AF", "IL-10","Lyve1", "CD90", "LAG3"),
                hm_pal = brewer.pal(9, "Greys")
)


#ridgeplots of markers on merged clusters
plotClusterExprs(daf_alt_clean_rerun, k = "merging1")

#save main analysis object
saveRDS(daf_X_sf7, file = "./daf_X_sf7.rds")

#-------------------------------compare groups now (w/o removing IV or IL-10 cell groups------------------------------##
#plot abundances as barplot
plotAbundances(daf_X_sf7, k = "merging1", by = "sample_id", k_pal = clust_color)


#plot abundances as boxplot
plotAbundances(daf_X_sf7, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
               )


##To analyze now we can generate a GLMM
ei <- metadata(daf_X_sf7)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

#how to access merged cluster levels
levels(cluster_ids(daf_X_sf7, k = "merging1"))

#cluster levels of FlowSOM clusters
levels(cluster_ids(daf_X_sf7))

da_res1 <- diffcyt(daf_X_sf7,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  

#examine output
names(da_res1)

cluster_stats <- rowData(da_res1$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res1$res)$p_adj < FDR_cutoff)

#plot expression of markers by cluster (boxplot)
plotPbExprs(daf_X_sf7,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "genotype",
            group_by = "genotype",
            geom = "both", 
            jitter = TRUE
                   )


#make dup of SCE object to perserve it
daf_X_sf7_test <- daf_X_sf7
rowData(daf_X_sf7_test)$marker_class <- "state"

ei <- metadata(daf_X_sf7_test)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(daf_X_sf7_test), cols_design = "genotype")

da_res99 <- diffcyt(daf_X_sf7_test,                                            
                    formula = da_formula99, contrast = contrast,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-limma",  
                    design = design99,
                    #markers_to_test = all_mark,  #works, except I can't specify which markers using this command
                    clustering_to_use = "merging1", verbose = FALSE)  

#examine output
cluster_stats99 <- rowData(da_res99$res) %>%
  as.data.frame(.)

table(rowData(da_res99$res)$p_adj < FDR_cutoff)

tbl_DS <- rowData(da_res99$res)

#plot heatmap of results
plotDiffHeatmap(daf_X_sf7_test, tbl_DS, fdr = 0.05, 
                sort_by = "lfc", col_anno = "genotype", all = TRUE)

#make specific for NK, CD8, CD4 T, B cells
k <- metadata(da_res99$res)$clustering_name
sub <- filterSCE(daf_IL10_only, cluster_id %in% c("NK cells", "CD8+ T cells", "CD4+ T cells", "B cells"), k = k)
plotDiffHeatmap(sub, tbl_DS, all = TRUE, normalize = FALSE)

#alt heatmap comparison for above
#plot heatmap of results
plotDiffHeatmap(daf_X_sf7_test, tbl_DS, fdr = 0.05, all = T,
                sort_by = "padj", col_anno = "genotype")

##----------------------Subset out WT samples and generate map of SF7 expression-----------------------##

#load SCE object back in
daf_X_sf7 <- readRDS("daf_X_sf7.rds")

#remove SLAMF7-KO samples
daf_X_sf7_wt <- filterSCE(daf_X_sf7, genotype == "WT")
levels(daf_X_sf7_wt$genotype)

#UMAP of SLAMF7 expression
plotDR(daf_X_sf7_wt, "UMAP", color_by = "SLAMF7", a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) 




##---------now make histogram plots of SLAMF7 expression for each cluster comparing WT and SLAMF7-KO------##
#manually extract df containing merging1 cluster ids matched to SLAMF7 expression

#add merging1 cluster to metadata
daf_X_sf7$merg_clust <- cluster_ids(daf_X_sf7, "merging1")
#make df
ggdf <- data.frame(colData(daf_X_sf7), t(assay(daf_X_sf7, "exprs")), check.names = FALSE)
#make plot
fg <- ggplot(ggdf, aes_string("SLAMF7", fill = "genotype")) + 
  geom_density(alpha = 0.7) + 
  facet_wrap(ggdf$merg_clust) +
  scale_fill_manual(values = c("#FF3333", "#636363"))

fg + theme(axis.line = element_line(size = 0.5, 
    linetype = "solid"), panel.background = element_rect(fill = NA)) 

##------------------------------------------------------------------------------------------------------##


##----------------------see what fraction of each cluster/sample are IV+-----------------------##

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_X_sf7, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IV_clust <- table(
  IV_pos = daf_X_sf7$IV_pos, 
  cluster = cluster_ids(daf_X_sf7, "merging1")) 

plot_clust_table <- prop.table(IV_clust, "cluster") %>%
  as.data.frame(.)

#make plot
rr <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rr + theme(axis.ticks = element_line(colour = "black", 
    size = 1), axis.title = element_text(size = 15), 
    axis.text = element_text(colour = "black", 
        hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
    axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
    )
##---------------------------------------------------------------------------------------------##

#I'm going to leave the IV+ cells in for now. we can see what are the major immune cell and IL-10+ cell
#contributions to the genotypes in our dataset and can remove later should we choose. 

##---------------------Now lets look at the IL-10+ compartment during EAE----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_X_sf7, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IL10_clust <- table(
  IL10_pos = daf_X_sf7$IL_10_pos, 
  cluster = cluster_ids(daf_X_sf7, "merging1")) 

plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
cc <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
cc + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
##--------------------------------------------------------------------------------------------##




##---------------------Now lets reconstruct the IL-10+ compartment during EAE----------------------##

#remove IL_10_neg samples
daf_IL10_eae <- filterSCE(daf_X_sf7, IL_10_pos == "yes")

levels(daf_IL10_eae$condition)
levels(daf_IL10_eae$IL_10_pos)
levels(daf_IL10_eae$sample_id)

#update cluster colors since not all clusters are here
clust_color <- c("#996600", "#00CCAA", "#6698FD", "#FFCB64", "#8800CC", "#00994D", 
                 "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC")

clust_color3 <- c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                  "#FF97A7", "#006633", "#FF65CA", "#DD97FC")
                  

#UMAP of IL-10+ by genotype
plotDR(daf_IL10_eae, "UMAP", color_by = "merging1", k_pal = clust_color3, facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=1)

#comapre cell numbers
n_cells(daf_IL10_eae)
n_cells(daf_X_sf7)
table(cluster_ids(daf_IL10_eae, k = "merging1"))

IL10_eae_by_IV <- table(
  IL10_pos = daf_IL10_eae$IL_10_pos, 
  IV_pos = daf_IL10_eae$IV_pos, 
  genotype = daf_IL10_eae$genotype) 
IL10_eae_by_IV

#pie chart of cluster contributions to IL-10 compartment in WT and SLAMF7-KO (WT first)
`%notin%` <- Negate(`%in%`) #make this opposite of %in% operator to help
to_remove2 <- c("Neutrophils")
{IL10_pie_eae <- table(
  genotype = daf_IL10_eae$genotype, 
  cluster = cluster_ids(daf_IL10_eae, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT") %>%
  dplyr::filter(., cluster %notin% to_remove2)}

#make plot
hf <- ggplot(data=IL10_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                               "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
hf

#now SLAMF7-KO
{IL10_pie_eae_sf7 <- table(
  genotype = daf_IL10_eae$genotype, 
  cluster = cluster_ids(daf_IL10_eae, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "SLAMF7_KO") %>%
  dplyr::filter(., cluster %notin% to_remove2)}

#make plot
fg <- ggplot(data=IL10_pie_eae_sf7, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                               "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
fg

##Now compare cluster IL-10+ cluster frequency b/w WT and sf7-KO
eie <- metadata(daf_IL10_eae)$experiment_info 
(da_formula12 <- createFormula(eie,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res12 <- diffcyt(daf_IL10_eae,                                            
                    formula = da_formula12, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "merging1", verbose = FALSE)    

#examine output
names(da_res12)

cluster_stats2 <- rowData(da_res12$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res12$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##



##--------------------- IL-10+ compartment during EAE (IV+ removed, 7/14/2021)----------------------##
                                #This data is published in the paper#
#remove IL_10_neg samples and IV+
daf_IL10_eae <- filterSCE(daf_X_sf7, IL_10_pos == "yes", IV_pos == "no")

levels(daf_IL10_eae$condition)
levels(daf_IL10_eae$IL_10_pos)
levels(daf_IL10_eae$sample_id)
levels(daf_IL10_eae$IV_pos)


#update cluster colors since not all clusters are here
clust_color3 <- c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                  "#FF97A7", "#006633", "#FF65CA", "#DD97FC")


#UMAP of IL-10+ by genotype
plotDR(daf_IL10_eae, "UMAP", color_by = "merging1", k_pal = clust_color3, facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=1)

#comapre cell numbers
n_cells(daf_IL10_eae)
n_cells(daf_X_sf7)
table(cluster_ids(daf_IL10_eae, k = "merging1"))

IL10_eae_by_IV <- table(
  IL10_pos = daf_IL10_eae$IL_10_pos, 
  IV_pos = daf_IL10_eae$IV_pos, 
  genotype = daf_IL10_eae$genotype) 
IL10_eae_by_IV

#pie chart of cluster contributions to IL-10 compartment in WT and SLAMF7-KO (WT first)
`%notin%` <- Negate(`%in%`) #make this opposite of %in% operator to help
to_remove2 <- c("Neutrophils")
{IL10_pie_eae <- table(
  genotype = daf_IL10_eae$genotype, 
  cluster = cluster_ids(daf_IL10_eae, "merging1")) %>%
    as.data.frame(.) %>%
    dplyr::filter(., genotype == "WT") %>%
    dplyr::filter(., cluster %notin% to_remove2)}

#make plot
hf <- ggplot(data=IL10_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                               "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
hf

#now SLAMF7-KO
{IL10_pie_eae_sf7 <- table(
  genotype = daf_IL10_eae$genotype, 
  cluster = cluster_ids(daf_IL10_eae, "merging1")) %>%
    as.data.frame(.) %>%
    dplyr::filter(., genotype == "SLAMF7_KO") %>%
    dplyr::filter(., cluster %notin% to_remove2)}

#make plot
fg <- ggplot(data=IL10_pie_eae_sf7, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                               "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
fg

##Now compare cluster IL-10+ cluster frequency b/w WT and sf7-KO
ei <- metadata(daf_IL10_eae)$experiment_info 
(da_formula12 <- createFormula(ei,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res12 <- diffcyt(daf_IL10_eae,                                            
                    formula = da_formula12, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "merging1", verbose = FALSE)    

#examine output
names(da_res12)

cluster_stats2 <- rowData(da_res12$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res12$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##


##--------------------- IL-10+ compartment during EAE (only IV-, non-resident; 7/14/2021)----------------------##

#remove IL_10_neg samples and IV-
daf_IL10_eae1 <- filterSCE(daf_X_sf7, IL_10_pos == "yes", IV_pos == "yes")

levels(daf_IL10_eae1$condition)
levels(daf_IL10_eae1$IL_10_pos)
levels(daf_IL10_eae1$sample_id)
levels(daf_IL10_eae1$IV_pos)


#update cluster colors since not all clusters are here
clust_color3 <- c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF9864", 
                  "#FF97A7", "#006633", "#FF65CA", "#DD97FC")


#UMAP of IL-10+ by genotype
plotDR(daf_IL10_eae1, "UMAP", color_by = "merging1", k_pal = clust_color3, facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=1)

#comapre cell numbers
n_cells(daf_IL10_eae1)
n_cells(daf_X_sf7)
table(cluster_ids(daf_IL10_eae1, k = "merging1"))


#pie chart of cluster contributions to IL-10 compartment in WT and SLAMF7-KO (WT first)
to_remove2 <- c("Neutrophils", "BAMs", "pDCs", "CCR2+ monocytes", "CD8+ T cells")
{IL10_pie_eae1 <- table(
  genotype = daf_IL10_eae1$genotype, 
  cluster = cluster_ids(daf_IL10_eae1, "merging1")) %>%
    as.data.frame(.) %>%
    dplyr::filter(., genotype == "WT") %>%
    dplyr::filter(., cluster %notin% to_remove2)}

#make plot
hf <- ggplot(data=IL10_pie_eae1, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#6698FD", "#8800CC", 
                               "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
hf

#now SLAMF7-KO
{IL10_pie_eae_sf71 <- table(
  genotype = daf_IL10_eae1$genotype, 
  cluster = cluster_ids(daf_IL10_eae1, "merging1")) %>%
    as.data.frame(.) %>%
    dplyr::filter(., genotype == "SLAMF7_KO") %>%
    dplyr::filter(., cluster %notin% to_remove2)}

#make plot
fg <- ggplot(data=IL10_pie_eae_sf71, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#6698FD", "#8800CC", 
                              "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
fg

##Now compare cluster IL-10+ cluster frequency b/w WT and sf7-KO
ei <- metadata(daf_IL10_eae1)$experiment_info 
(da_formula12 <- createFormula(ei,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res12 <- diffcyt(daf_IL10_eae1,                                            
                    formula = da_formula12, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "merging1", verbose = FALSE)    

#examine output
names(da_res12)

cluster_stats2 <- rowData(da_res12$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res12$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##
