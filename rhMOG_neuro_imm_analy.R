#This analysis starts from FCS files and goes from there. W/ downsampling. 
# This is the neuroimmune profiling of SF7-KO and WT mice subjected to EAE w/ rhMOG35-55 peptide
# 06122022

#Packages
library(CATALYST)
library(tidyverse)
library(flowCore)
library(readxl)
library(ggplot2)
library(RColorBrewer)


##-------------------------Making metadata file----------------------------------##
#Load the content of the metadata.xlsx file into R (file name must match name of the fcs file you generate for each sample)
setwd("~/Documents/Research Projects/SLAMF7 neoruimmune project/SF7_rhMOG_EAE_neuroimmune_06032022")
metadata_filename <- "eae_metadata_alt.xlsx"
md <- read_excel(metadata_filename)

# Define condition variables as named in metadata
md$condition <- factor(md$condition, levels = c("EAE", "naive"))
md$genotype <- factor(md$genotype, levels = c("WT", "SLAMF7_KO"))
head(data.frame(md))
##-----------------------------------------------------------##

#Read FCS files in as flowset 
setwd("~/Documents/Research Projects/SLAMF7 neoruimmune project/SF7_rhMOG_EAE_neuroimmune_06032022/clean_living_fcs")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

##------------------making panel file-----------------------##
setwd("~/Documents/Research Projects/SLAMF7 neoruimmune project/SF7_rhMOG_EAE_neuroimmune_06032022")
panel <- "neuroimmune_panel_alt.xlsx"
panel <- read_excel(panel) 

#check
all(panel$fcs_colname %in% colnames(fcs_raw))
##------------------making panel file-----------------------##

#see how many cells per file
fsApply(fcs_raw, nrow)

#downsample to 60,000 cells per sample
# Define a downsampling ceiling
sampling.ceiling <- 60000
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


#First 3 must be present and sample_id must match condition
factors <- list(factors = c("file_name", "condition", "sample_id", "genotype"))

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


# view number of events per sample
table(daf_X$condition)
table(daf_X$sample_id)
table(daf_X$file_name)
table(daf_X$genotype)



#view metadata parameters
names(colData(daf_X))

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")  

#load RDS object
daf_X <- readRDS(file = "./daf_X.rds")

#plot global marker expression by condition
p <- plotExprs(daf_X, color_by = "condition")
p$facet$params$ncol <- 6                   
p 

#plot number of cells per condition
plotCounts(daf_X, color_by = "condition")
table(daf_X$condition) 

#Make MDS plot  
pbMDS(daf_X, by = "sample_id", color_by = "genotype", shape_by = "condition")


#Perform FlowSOM clustering
daf_X <- cluster(daf_X, 
               features = "type", 
               xdim = 10, 
               ydim = 10, 
               maxK = 22, 
               seed = 1234, 
               verbose = TRUE
               ) 

#plot cluster heatmap
plotClusterHeatmap(daf_X,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta22", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
                   ) 

# run UMAP                          
set.seed(777)                               
daf_X <- runDR(daf_X, 
             dr = "UMAP", 
             cells = 10000,
             features = "type" 
             )

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")

#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#UMAP visualizations
plotDR(daf_X, "UMAP", color_by = "meta22")

pz <- plotDR(daf_X, "UMAP", color_by = "meta22")               
pz + theme(axis.line = element_line(colour = NA),   #use this theme for all plots
    axis.ticks = element_line(colour = NA), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(colour = NA), 
    axis.text = element_text(colour = NA), 
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 14), 
    panel.background = element_rect(fill = NA)) 



#Merge and annotate clusters
merging_tableX <- "merged_clusters_X.xlsx"
merging_table1 <- read_excel(merging_tableX)
head(data.frame(merging_table1))  

# convert to factor with merged clusters in desired order                
merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("Microglia", "BAMs", "MdCs", "pDCs", "DCs", "CCR2+ monocytes", "Neutrophils",             
                                                "CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Debris"))        

daf_X <- mergeClusters(daf_X, k = "meta22",                                  
                     table = merging_table1, 
                     id = "merging1",
                     overwrite = TRUE
                     )   

#remove debris cluster
#make a new SCE object when you do this!!
daf_X_clean <- filterSCE(daf_X, cluster_id != "Debris", k = "merging1")

#change cluster colors
clust_color <- c("#6B5B95", "#ECDB54", "#E94B3C", "#6F9FD8", "#944743", "#DBB1CD", "#EC9787", "#00A591", "#6C4F3D", "#BC70A4", "#BFD641")

#UMAP of merged clusters  Print at 1600x1400
plotDR(daf_X_clean, "UMAP", color_by = "merging1", k_pal = clust_color) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
                                                                axis.ticks = element_line(colour = NA), 
                                                                panel.grid.major = element_line(linetype = "blank"), 
                                                                panel.grid.minor = element_line(linetype = "blank"), 
                                                                axis.title = element_text(colour = NA), 
                                                                axis.text = element_text(colour = NA), 
                                                                legend.text = element_text(size = 10), 
                                                                legend.title = element_text(size = 14), 
                                                                panel.background = element_rect(fill = NA)) +
                                                                geom_point(size=0.05)


#UMAP of merged clusters facet by genotype Print at 3000x1400
plotDR(daf_X_clean, "UMAP", color_by = "merging1", k_pal = clust_color, facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.05)

#UMAP of merged clusters facet by condition Print at 3000x1400
plotDR(daf_X_clean, "UMAP", color_by = "merging1", k_pal = clust_color, facet_by = "condition") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.05)

#heatmap of merged cluster markers
plotExprHeatmap(daf_X_clean, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color, 
                hm_pal = brewer.pal(9, "Greys")
)

#-------------------------------compare genotypes now (naive first)------------------------------##
library(diffcyt)

#remove all EAE samples
daf_X_naive <- filterSCE(daf_X_clean, condition == "naive")

#plot abundances as barplot
plotAbundances(daf_X_naive, k = "merging1", by = "sample_id", k_pal = clust_color)


#plot abundances as boxplot
obj66 <- plotAbundances(daf_X_naive, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
               )

obj66 + scale_color_manual(values = c("#000000", "#ff3333")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)


##To analyze now we can generate a GLMM
ei <- metadata(daf_X_naive)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res1 <- diffcyt(daf_X_naive,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  

#examine output
names(da_res1)

cluster_stats <- rowData(da_res1$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res1$res)$p_adj < FDR_cutoff)

#Show all results
topTable(da_res1, show_props = TRUE, format_vals = TRUE, digits = 2)

#plot expression of markers by cluster (boxplot)
plotPbExprs(daf_X_naive,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "genotype",
            group_by = "genotype",
            geom = "both", 
            jitter = TRUE
                   )


##now compare markers on all clusters
#now trying to manually set all markers to be "state" markers on the raw SCE object as a work around

#make dup of SCE object to perserve it
daf_X_sf7_test <- daf_X_naive
rowData(daf_X_sf7_test)$marker_class <- "state"

ei <- metadata(daf_X_sf7_test)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "genotype",      
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
                sort_by = "padj", col_anno = "genotype", all = TRUE)

##------------------------------------------------------------------------------------##


#-------------------------------compare genotypes now (EAE)------------------------------##
#load main analysis object back in
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
#remove all naive samples
daf_X_eae <- filterSCE(daf_X_clean, condition == "EAE")

#plot abundances as barplot
plotAbundances(daf_X_eae, k = "merging1", by = "sample_id", k_pal = clust_color)


#plot abundances as boxplot
obj77 <- plotAbundances(daf_X_eae, 
                        k = "merging1", 
                        by = "cluster_id",
                        group_by = "genotype",
                        shape = NULL,
)

obj77 + scale_color_manual(values = c("#000000", "#ff3333")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)


##To analyze now we can generate a GLMM
ei <- metadata(daf_X_eae)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res2 <- diffcyt(daf_X_eae,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  


cluster_stats2 <- rowData(da_res2$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res2$res)$p_adj < FDR_cutoff)

#Show all results
topTable(da_res2, show_props = TRUE, format_vals = TRUE, digits = 2)

#plot expression of markers by cluster (boxplot)
plotPbExprs(daf_X_eae,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "genotype",
            group_by = "genotype",
            geom = "both", 
            jitter = TRUE
)


##now compare markers on all clusters
#now trying to manually set all markers to be "state" markers on the raw SCE object as a work around

#make dup of SCE object to perserve it
daf_X_sf7_test2 <- daf_X_eae
rowData(daf_X_sf7_test2)$marker_class <- "state"

ei <- metadata(daf_X_sf7_test2)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(daf_X_sf7_test2), cols_design = "genotype")

da_res99 <- diffcyt(daf_X_sf7_test2,                                            
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
plotDiffHeatmap(daf_X_sf7_test2, tbl_DS, fdr = 0.05, lfc = 0.2,
                sort_by = "padj", col_anno = "genotype", all = TRUE)


##------------------------------------------------------------------------------------##
