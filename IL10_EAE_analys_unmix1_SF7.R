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

##---------------------Now subset B cells from complete exp and analyze----------------------##

#subset B cells
daf_B <- filterSCE(daf_X_sf7, cluster_id == "B cells", k = "merging1")

levels(daf_B$condition)
levels(daf_B$IL_10_pos)
levels(daf_B$sample_id)

#Make MDS plot  
pbMDS(daf_B, by = "sample_id", color_by = "genotype", size_by = TRUE)

#re-cluster just the B cells
daf_B <- cluster(daf_B, 
                 features = NULL, 
                 xdim = 10, 
                 ydim = 10, 
                 maxK = 8, 
                 seed = 7689, 
                 verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_B,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta8", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(493)                               
daf_B <- runDR(daf_B, 
               dr = "UMAP", 
               cells = 10000,
               features = NULL 
)

#UMAP B cells    
plotDR(daf_B, "UMAP", color_by = "meta8", facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

# run tSNE                          
set.seed(493)                               
daf_B <- runDR(daf_B, 
               dr = "TSNE", 
               cells = 10000,
               features = NULL 
)

#tSNE B cells    
plotDR(daf_B, "TSNE", color_by = "meta8") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)



#Merge and annotate B clusters
merging_table_B <- "merged_clusters_B.xlsx"
merging_table2 <- read_excel(merging_table_B)
head(data.frame(merging_table2))  

# convert to factor with merged clusters in desired order                
merging_table2new_cluster <- factor(merging_table2$new_cluster,         
                                     levels = c("Debris", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6"))        

daf_B <- mergeClusters(daf_B, k = "meta8",                                  
                       table = merging_table2, 
                       id = "B_clusters",
                       overwrite = TRUE
)   

#remove debris cluster
#make a new SCE object when you do this!!
daf_B_clean <- filterSCE(daf_B, cluster_id != "Debris", k = "B_clusters")

#save objects
saveRDS(daf_B, file = "./daf_B.rds")
saveRDS(daf_B_clean, file = "./daf_B_clean.rds")
daf_B_clean <- readRDS("./daf_B_clean.rds")
#change cluster colors
clust_color4 <- c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D")

#UMAP B cells
plotDR(daf_B_clean, "UMAP", color_by = "B_clusters") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged B cluster markers
plotExprHeatmap(daf_B_clean, 
                bin_anno = FALSE, 
                scale = "first",
                q = 0.1,
                k = "B_clusters",
                by = "cluster_id", 
                features = c("SLAMF7", "CD38", "NKG2D", "IL-10", "Ly6C", "AF", "CD19", "IgD", "MHCII", "SiglecH", "B220", "CCR2", "Tim3", "LAG3"),
                hm_pal = brewer.pal(9, "Greys")
)

#make UMAP colored by important markers on just B cells
plotDR(daf_B_clean, "UMAP", color_by = c("IL-10", "B220", "IgD", "SLAMF7","CCR2", "MHCII", "Tim3", "CD38", "Ly6C"), facet_by = NULL, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#UMAP B cells show IL-10_pos
plotDR(daf_B_clean, "UMAP", color_by = "IL_10_pos", a_pal = c("#FF9864", "#000000")) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA),                 #come back to this to make it pretty
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)





#------------------compare B cluster frequencies by genotype (EAE)-----------------------------#

#plot abundances as barplot
plotAbundances(daf_B_clean, k = "B_clusters", by = "sample_id")


#plot abundances as boxplot
plotAbundances(daf_B_clean, 
               k = "B_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##To analyze now we can generate a GLMM
eim <- metadata(daf_B_clean)$experiment_info 
(da_formula14 <- createFormula(eim,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res14 <- diffcyt(daf_B_clean,                                            
                    formula = da_formula14, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "B_clusters", verbose = FALSE)  

#examine output
cluster_stats4 <- rowData(da_res14$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res14$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##

##----------------------see what fraction of B cells from each cluster/sample are IV+-----------------------##

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_B_clean, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IV_clust_b <- table(
  IV_pos = daf_B_clean$IV_pos, 
  cluster = cluster_ids(daf_B_clean, "B_clusters")) 

plot_clust_table <- prop.table(IV_clust_b, "cluster") %>%
  as.data.frame(.)

#make plot
rg <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rg + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
##---------------------------------------------------------------------------------------------##

##---------------------Now lets see what fraction of each cluster are IL_10+----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_B_clean, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IL10_clust <- table(
  IL10_pos = daf_B_clean$IL_10_pos, 
  cluster = cluster_ids(daf_B_clean, "B_clusters")) 

plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
cn <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
cn + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
##--------------------------------------------------------------------------------------------##



##----------------------Compare IL-10+ B cells b/w genotypes by IV+ and cluster-----------------------##
#load main analysis object back in
daf_B_clean <- readRDS(file = "./daf_B_clean.rds")  


#Start w/ EAE
#remove IL-10 neg cells
daf_B_clean_IL10 <- filterSCE(daf_B_clean, IL_10_pos == "yes")

#pie chart of cluster contributions to IL-10+ NK cells in WT and SLAMF7-KO (wt first)
IL10_b_pie_eae <- table(
  genotype = daf_B_clean_IL10$genotype, 
  cluster = cluster_ids(daf_B_clean_IL10, "B_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hk <- ggplot(data=IL10_b_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D"))
hk

#pie chart of cluster contributions to IL-10+ B cells in WT and SLAMF7-KO (SLAMF7 now)
IL10_b_pie_eae2 <- table(
  genotype = daf_B_clean_IL10$genotype, 
  cluster = cluster_ids(daf_B_clean_IL10, "B_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "SLAMF7_KO")

#make plot
hl <- ggplot(data=IL10_b_pie_eae2, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D"))
hl

#plot abundances as boxplot
dff <- plotAbundances(daf_B_clean_IL10, 
               k = "B_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
) + ylim(0,NA) + scale_color_manual(values = c("#000000", "#ff1a1a")) + scale_fill_manual(values = c("#ffffff", "#ffffff"))
dff$facet$params$ncol <- 3   
dff

#heatmap of merged B cluster markers (IL-10+ only)
plotExprHeatmap(daf_B_clean_IL10, 
                bin_anno = FALSE, 
                scale = "first",
                q = 0.1,
                k = "B_clusters",
                by = "cluster_id", 
                features = c("SLAMF7", "CD38", "NKG2D", "IL-10", "Ly6C", "AF", "CD19", "IgD", "MHCII", "SiglecH", "B220", "CCR2", "Tim3", "LAG3"),
                hm_pal = brewer.pal(9, "Greys")
)


n_cells(daf_B_clean_IL10)

##Now compare cluster B IL-10+ cluster frequency b/w WT and SLAMF7-KO

eiz <- metadata(daf_B_clean_IL10)$experiment_info 
(da_formula17 <- createFormula(eiz,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res17 <- diffcyt(daf_B_clean_IL10,                                            
                    formula = da_formula17, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "B_clusters", verbose = FALSE)    

#examine output
cluster_stats7 <- rowData(da_res17$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res17$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------------##


                 
##---------------------Now subset T cells from complete exp and analyze----------------------##
#load main analysis object back in
daf_X_sf7 <- readRDS(file = "./daf_X_sf7.rds")

#subset T cells
daf_T <- filterSCE(daf_X_sf7, cluster_id %in% c("CD8+ T cells", "CD4+ T cells"), k = "merging1")

levels(daf_T$condition)
levels(daf_T$IL_10_pos)
levels(daf_T$sample_id)

#Make MDS plot  
pbMDS(daf_T, by = "sample_id", color_by = "genotype", size_by = TRUE)

#re-cluster just the T cells     
daf_T <- cluster(daf_T, 
                  features = NULL, 
                  xdim = 10, 
                  ydim = 10, 
                  maxK = 15, 
                  seed = 7689, 
                  verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_T,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta15", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(797)                               
daf_T <- runDR(daf_T, 
                dr = "UMAP", 
                cells = 10000,
                features = NULL 
)

#UMAP T cells
plotDR(daf_T, "UMAP", color_by = "meta15", facet_by = "IL_10_pos") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#Merge and annotate T clusters
merging_table_T <- "merged_clusters_T.xlsx"
merging_table3 <- read_excel(merging_table_T)
head(data.frame(merging_table3))  

# convert to factor with merged clusters in desired order                
merging_table3$new_cluster <- factor(merging_table3$new_cluster,         
                                    levels = c("Debris", "CD4_1", "CD4_2", "CD4_Treg_1", "CD4_Treg_2", 
                                               "CD8_1", "CD8_2", "CD8_3", "NKT cells", "DN_T cells", "MHC-II+/CD11c+ CD4+ T cells"))        

daf_T <- mergeClusters(daf_T, k = "meta15",                                  
                        table = merging_table3, 
                        id = "T_clusters",
                        overwrite = TRUE
)   

#remove debris cluster
#make a new SCE object when you do this!!
daf_T_clean <- filterSCE(daf_T, cluster_id != "Debris", k = "T_clusters")

#save objects
saveRDS(daf_T, file = "./daf_T.rds")
saveRDS(daf_T_clean, file = "./daf_T_clean.rds")

#load main analysis object back in
daf_T <- readRDS(file = "./daf_T.rds")

#change cluster colors
clust_color5 <- c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#BB32FE", 
                  "#CC0022", "#FF9999", "#CC00CC", "#636363")

#UMAP T cells 
plotDR(daf_T_clean, "UMAP", color_by = "T_clusters", k_pal = clust_color5, facet_by = "IL_10_pos") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#UMAP T cells selected markers
plotDR(daf_T_clean, "UMAP", color_by = c("Ly6C", "CD4", "CD8", "CD38", "CD11b", "NK1.1", "LAG3", "IL-10"), ncol = 2, k_pal = clust_color5, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged T cluster markers
plotExprHeatmap(daf_T_clean, 
                bin_anno = FALSE, 
                scale = "first",
                k = "T_clusters",
                q = 0.1,
                by = "cluster_id", 
                k_pal = clust_color5,
                features = c("CD38", "CD8", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD4", "AF", "LAG3", "Tim3", "CD49b", "NK1.1"),
                hm_pal = brewer.pal(9, "Greys")
)



#------------------compare T cluster frequencies by genotype (EAE)-----------------------------#

#plot abundances as barplot
plotAbundances(daf_T_clean, k = "T_clusters", by = "sample_id", k_pal = clust_color5)


#plot abundances as boxplot
plotAbundances(daf_T_clean, 
               k = "T_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color5,
               shape = NULL
)


##To analyze now we can generate a GLMM
eib <- metadata(daf_T_clean)$experiment_info 
(da_formula20 <- createFormula(eib,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 
contrast <- createContrast(c(0, 1))

da_res20 <- diffcyt(daf_T_clean,                                            
                    formula = da_formula20, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "T_clusters", verbose = FALSE)  

#examine output
cluster_stats40 <- rowData(da_res20$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res20$res)$p_adj < FDR_cutoff)

#now comparing markers on all clusters
#make dup of SCE object to preserve it
daf_T_clean_test <- daf_T_clean
rowData(daf_T_clean_test)$marker_class <- "state"

ei <- metadata(daf_T_clean_test)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(daf_T_clean_test), cols_design = "genotype")

da_res99 <- diffcyt(daf_T_clean_test,                                            
                    formula = da_formula99, contrast = contrast,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-limma",  
                    design = design99,
                    clustering_to_use = "T_clusters", verbose = FALSE)  

#examine output
cluster_stats99 <- rowData(da_res99$res) %>%
  as.data.frame(.)

table(rowData(da_res99$res)$p_adj < FDR_cutoff)

tbl_DS <- rowData(da_res99$res)

#plot heatmap of results
plotDiffHeatmap(daf_T_clean_test, tbl_DS, fdr = 0.05, 
                sort_by = "lfc", col_anno = "genotype", all = TRUE)

##-----------------------------------------------------------------------------------##

##----------------------see what fraction of T cells from each cluster/sample are IV+-----------------------##

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_T_clean, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IV_clust_T <- table(
  IV_pos = daf_T_clean$IV_pos, 
  cluster = cluster_ids(daf_T_clean, "T_clusters")) 

plot_clust_table <- prop.table(IV_clust_T, "cluster") %>%
  as.data.frame(.)

#make plot
rn <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rn + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)}
##---------------------------------------------------------------------------------------------##

##---------------------Now lets see what fraction of each cluster are IL_10+----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_T_clean, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IL10_clust <- table(
  IL10_pos = daf_T_clean$IL_10_pos, 
  cluster = cluster_ids(daf_T_clean, "T_clusters")) 

plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
cv <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
cv + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)}
##--------------------------------------------------------------------------------------------##


##----------------------Compare IL-10+ T cells b/w genotypes by IV+ and cluster-----------------------##

#now during EAE
#remove IL-10 neg cells
daf_T_clean_IL10_eae <- filterSCE(daf_T_clean, IL_10_pos == "yes")

#pie chart of cluster contributions to IL-10+ T cells in WT and SF7-KO (eae) (wt first)
IL10_T_pie_eae <- table(
  genotype = daf_T_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_T_clean_IL10_eae, "T_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hf <- ggplot(data=IL10_T_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#BB32FE", 
                               "#CC0022", "#FF9999", "#CC00CC", "#636363"))
hf

#pie chart of cluster contributions to IL-10+ T cells in WT and SF7-KO (eae) (SF7 now)
IL10_T_pie_eae2 <- table(
  genotype = daf_T_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_T_clean_IL10_eae, "T_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "SLAMF7_KO")

#make plot
hz <- ggplot(data=IL10_T_pie_eae2, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#BB32FE", 
                               "#CC0022", "#FF9999", "#CC00CC", "#636363"))
hz

#plot abundances as boxplot
plotAbundances(daf_T_clean_IL10_eae, 
               k = "T_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##Now compare cluster T IL-10+ cluster frequency b/w WT and SF7-KO    
eif <- metadata(daf_T_clean_IL10_eae)$experiment_info 
(da_formula39 <- createFormula(eif,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res39 <- diffcyt(daf_T_clean_IL10_eae,                                            
                    formula = da_formula39, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "T_clusters", verbose = FALSE, )  

#examine output
cluster_stats39 <- rowData(da_res39$res) %>%
  as.data.frame(.)

table(rowData(da_res39$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------------##




##------------------Now subset microglia-----------------------------------------------------##

#load main analysis object back in
daf_X_sf7 <- readRDS(file = "./daf_X_sf7.rds")

#subset microglia cells
daf_micro <- filterSCE(daf_X_sf7, cluster_id == "Microglia", k = "merging1")

levels(daf_micro$condition)
levels(daf_micro$IL_10_pos)
levels(daf_micro$sample_id)

#Make MDS plot  
pbMDS(daf_micro, by = "sample_id", color_by = "genotype", size_by = TRUE)

#re-cluster just the microglia     
daf_micro <- cluster(daf_micro, 
                 features = NULL, 
                 xdim = 10, 
                 ydim = 10, 
                 maxK = 6, 
                 seed = 7089, 
                 verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_micro,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta6", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(790)                               
daf_micro <- runDR(daf_micro, 
               dr = "UMAP", 
               cells = 10000,
               features = NULL 
)

#UMAP microglia
plotDR(daf_micro, "UMAP", color_by = "meta6", facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#save objects    
saveRDS(daf_micro, file = "./daf_micro.rds")

#load object back
daf_micro <- readRDS(file = "./daf_micro.rds")

#Merge and annotate microglia clusters
merging_table_micro <- "merged_clusters_micro.xlsx"
merging_table4 <- read_excel(merging_table_micro)
head(data.frame(merging_table4))  

# convert to factor with merged clusters in desired order                
merging_table4$new_cluster <- factor(merging_table4$new_cluster,         
                                     levels = c("Debris", "BAMs", "Microglia_1", "Microglia_2"))        

daf_micro <- mergeClusters(daf_micro, k = "meta6",                                  
                       table = merging_table4, 
                       id = "micro_clusters",
                       overwrite = TRUE
)   

#remove debris  and BAMs 
#make a new SCE object when you do this!!
daf_micro_clean <- filterSCE(daf_micro, cluster_id %in% c("Microglia_1", "Microglia_2", "Microglia_3"), k = "micro_clusters")

#save objects
saveRDS(daf_micro, file = "./daf_micro.rds")
saveRDS(daf_micro_clean, file = "./daf_micro_clean.rds")

#load back in
daf_micro_clean <- readRDS("./daf_micro_clean.rds")

#change cluster colors
clust_color6 <- c("#66E4FD", "#FF97A7")

#UMAP microglia
plotDR(daf_micro_clean, "UMAP", color_by = "micro_clusters", facet_by = "IL_10_pos", k_pal = clust_color6, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged micro cluster markers
plotExprHeatmap(daf_micro_clean, 
                bin_anno = FALSE, 
                scale = "first",
                k = "micro_clusters",
                by = "cluster_id", 
                k_pal = clust_color6,
                features = c("CD11b", "CD4", "Ly6C", "IL-10", "SiglecH", "CD11c", "B220", "CD38", "AF", "CD206", "Tim3", "Lyve1", "NKG2D", "CD90"),
                hm_pal = brewer.pal(9, "Greys")
)

#UMAP microglia (by important markers)
plotDR(daf_micro_clean, "UMAP", color_by = c("B220", "CD90", "NKG2D", "CD4", "CD38", "Ly6C", "SiglecH", "AF", "CD11c", "Tim3"), k_pal = clust_color6, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)




#------------------compare microglia cluster frequencies by genotype (EAE)-----------------------------#

#plot abundances as barplot
plotAbundances(daf_micro_clean, k = "micro_clusters", by = "sample_id", k_pal = clust_color6)


#plot abundances as boxplot
plotAbundances(daf_micro_clean, 
               k = "micro_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color6,
               shape = NULL
)


##To analyze now we can generate a GLMM     DIDNT DO. NO OBVIOUS DIFF. 
eiv <- metadata(daf_micro_clean_eae)$experiment_info 
(da_formula31 <- createFormula(eiv,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res31 <- diffcyt(daf_micro_clean_eae,                                            
                    formula = da_formula31, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "micro_clusters", verbose = FALSE)  

#examine output
cluster_stats31 <- rowData(da_res31$res) %>%
  as.data.frame(.)

table(rowData(da_res31$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##




##----------------------see what fraction of microglia from each cluster/sample are IV+ (eae)-----------------------##

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_micro_clean, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IV_clust_micro <- table(
  IV_pos = daf_micro_clean$IV_pos, 
  cluster = cluster_ids(daf_micro_clean, "micro_clusters")) 
  
  plot_clust_table <- prop.table(IV_clust_micro, "cluster") %>%
    as.data.frame(.)
  
  #make plot
  rp <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#00CC22", "#000000"))
  
  #make purrty
  rp + theme(axis.ticks = element_line(colour = "black", 
                                       size = 1), axis.title = element_text(size = 15), 
             axis.text = element_text(colour = "black", 
                                      hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
             axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
  )}
##---------------------------------------------------------------------------------------------##

##---------------------Now lets see what fraction of each cluster are IL_10+ (eae)----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_micro_clean, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IL10_clust <- table(
  IL10_pos = daf_micro_clean$IL_10_pos, 
  cluster = cluster_ids(daf_micro_clean, "micro_clusters")) 
  
  plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
    as.data.frame(.)
  
  #make plot
  ct <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#FF9864", "#000000"))
  
  #make purrty
  ct + theme(axis.ticks = element_line(colour = "black", 
                                       size = 1), axis.title = element_text(size = 15), 
             axis.text = element_text(colour = "black", 
                                      hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
             axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
  )}

##--------------------------------------------------------------------------------------------##


##----------------------Compare IL-10+ T cells b/w genotypes by IV+ and cluster (eae)-----------------------##
# during EAE
#remove IL-10 neg cells
daf_micro_clean_IL10_eae <- filterSCE(daf_micro_clean, IL_10_pos == "yes")

#pie chart of cluster contributions to IL-10+ microglia in WT and ERAP1-KO (eae) (wt first)
{IL10_micro_pie_eae <- table(
  genotype = daf_micro_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_micro_clean_IL10_eae, "micro_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hd <- ggplot(data=IL10_micro_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#66E4FD", "#FF97A7", "#9965FD"))
hd}

#pie chart of cluster contributions to IL-10+ T cells in WT and ERAP1-KO (eae) (SF7 now)
{IL10_micro_pie_eae2 <- table(
  genotype = daf_micro_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_micro_clean_IL10_eae, "micro_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "SLAMF7_KO")

#make plot
hz <- ggplot(data=IL10_micro_pie_eae2, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#66E4FD", "#FF97A7", "#9965FD"))
hz}

#plot abundances as boxplot
plotAbundances(daf_micro_clean_IL10_eae, 
               k = "micro_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##Now compare cluster microglia IL-10+ cluster frequency b/w WT and SF7-KO   ignore, all from same clust
eig <- metadata(daf_micro_clean_IL10_eae)$experiment_info 
(da_formula34 <- createFormula(eig,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res34 <- diffcyt(daf_micro_clean_IL10_eae,                                            
                    formula = da_formula34, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "micro_clusters", verbose = FALSE, transform = FALSE)  

#examine output
cluster_stats34 <- rowData(da_res34$res) %>%
  as.data.frame(.)

table(rowData(da_res34$res)$p_adj < FDR_cutoff)

##--------------------------------------------------------------------------------------------##



##--------------------- IL-10+ compartment during EAE (IV+ removed, 7/14/2021)----------------------##

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
