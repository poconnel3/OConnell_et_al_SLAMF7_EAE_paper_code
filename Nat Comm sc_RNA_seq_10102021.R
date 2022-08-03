###Re-analysis of Nat Comm MS paper that used scRNA-seq. 
#Here I am splitting cells up by location, and from each location I am splitting cells into SF7+ and SF7-
#I will then compare DEGs and pathways b/w SF7+ and neg in just B cells (from all, HC only, and MS only)
#I will then run a cell-cell interaction prediciton package and compare the interactions and receptor/ligand pairs
#b/w SF7+ and SF7- cells (splititng up by location) to see cell-cell interactions enriched w/ SF7 expression

##10/10/2021
#Note: this is 10X 3' data


setwd("~/Documents/Research Projects/SLAMF7 neoruimmune project/scRNA_seq_interact_preditct")
library(Seurat)
library(tidyverse)
#Note: I analyzed this dataset previously and will just load in a seurat object w/ clean cells and annotations. 

#load Seurat object
csf2 <- readRDS("./csf2.rds")
head(csf2@meta.data)
#----------------Now we starting analyzing the annotated dataset----------------------#

#plot cells by cluster
DimPlot(csf2, label = FALSE, pt.size = 0.01)
#plot cells by location
DimPlot(csf2, group.by = "location", label = FALSE)
# How many cells are in each group
table(Idents(csf2))
# How many cells are in each case type?
table(csf2$condition)
# clusters by condition and location
table(Idents(csf2), csf2$condition)
table(Idents(csf2), csf2$location)
prop.table(table(Idents(csf2), csf2$condition), margin = 2)


#compare SF7 expression by location and condition
DefaultAssay(csf2) <- "RNA"
VlnPlot(csf2, features = c("SLAMF7"), slot = "counts", split.by = "location", assay = "RNA")
VlnPlot(csf2, features = c("SLAMF7"), slot = "counts", split.by = "condition", assay = "RNA")

#look at prohibin-2 expression (1/8/2022) Nothing intesting. 
FeaturePlot(csf2, features = c("rna_SLAMF7", "PHB2"), blend = TRUE, pt.size = 0.01)
FeatureScatter(csf2, feature1 = "rna_SLAMF7", feature2 = "PHB2", slot = "data")

#plot all cells by cluster, just CSF (nice colors)
nice_pal <- c("#79C360",    "#3F8EAA",    "#E52829",    "#9471B4",   "#A6CEE3",     "#ED8F47", 
              "#B89B74",    "#FDB762",   "#B15928", "#DDD399")
DimPlot(csf.only, label = FALSE, pt.size = 0.01, cols = nice_pal)





#---------Predict cell-cell interactions for all SF7+ cells and all SF7- cells---------#

#this package look like it might be nice since it integrates many approaches:https://saezlab.github.io/liana/articles/liana_tutorial.html 
#and this looks like a nice way to prioritize important interactions and visualize: https://github.com/CostaLab/CrossTalkeR
library(CrossTalkeR)
library(liana)  


#---lets start w/ CellPhoneDB---
#Since in addition to B cells, I want to know which cell interactions are enriched in SF7+ cells 
#I am following this: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md
#Taking CellPhoneDB output and using CellTalkeR to compare and visualize



#mark all cells as SLAMF7 pos or neg
sf7_pos <- subset(x = csf2, subset = SLAMF7 > 0.5) #positive cutoff is heuristic based on visual inspection of vionlin plots
pos_names <- rownames(sf7_pos@meta.data)
# Mutate a column in original  metadata
csf2$barcode <- rownames(csf2@meta.data)
csf2@meta.data <- csf2@meta.data %>% 
  mutate(SLAMF7 = ifelse((csf2$barcode %in% pos_names), "Pos",  "Neg"))
saveRDS(csf2, "./csf2.rds")

#CSF SF7+
csf.only <- subset(csf2, subset = location == "CSF")
DefaultAssay(csf.only) <- "RNA"
csf.SF7.pos <- subset(csf.only, subset = SLAMF7 == "Pos")
#save cell cluster names as metadata to help squidpy run
csf.SF7.pos[["clust.names"]] <- Idents(object = csf.SF7.pos)
csf.SF7.pos <- SetIdent(csf.SF7.pos, value = "clust.names")
# Run liana
liana.csf.SF7pos <- liana_wrap(csf.SF7.pos,
                         method = c("sca", "natmi", "logfc", "connectome",
                                    "cellchat", "squidpy", "cellphonedb"),
                         resource = c('OmniPath'),
                         # CellChat requires normalized data
                         cellchat.params = list(.normalize=TRUE))  
saveRDS(liana.csf.SF7pos, "./liana.csf.SF7pos.rds")
liana.csf.SF7pos %>% glimpse()

#CSF SF7-
csf.SF7.neg <- subset(csf.only, subset = SLAMF7 == "Neg")
#save cell cluster names as metadata to help squidpy run
csf.SF7.neg[["clust.names"]] <- Idents(object = csf.SF7.neg)
csf.SF7.neg <- SetIdent(csf.SF7.neg, value = "clust.names")
# Run liana
liana.csf.SF7neg <- liana_wrap(csf.SF7.neg,
                               method = c("sca", "natmi", "logfc", "connectome",
                                          "cellchat", "squidpy", "cellphonedb"),
                               resource = c('OmniPath'),
                               # CellChat requires normalized data
                               cellchat.params = list(.normalize=TRUE))  
saveRDS(liana.csf.SF7neg, "./liana.csf.SF7neg.rds")
liana.csf.SF7neg %>% glimpse()

#PBMC SF7+
pbmc.only <- subset(csf2, subset = location == "PBMC")
DefaultAssay(pbmc.only) <- "RNA"
pbmc.SF7.pos <- subset(pbmc.only, subset = SLAMF7 == "Pos")
#save cell cluster names as metadata to help squidpy run
pbmc.SF7.pos[["clust.names"]] <- Idents(object = pbmc.SF7.pos)
pbmc.SF7.pos <- SetIdent(pbmc.SF7.pos, value = "clust.names")
# Run liana
liana.pbmc.SF7pos <- liana_wrap(pbmc.SF7.pos,
                               method = c("sca", "natmi", "logfc", "connectome",
                                          "cellchat", "squidpy", "cellphonedb"),
                               resource = c('OmniPath'),
                               # CellChat requires normalized data
                               cellchat.params = list(.normalize=TRUE))  
saveRDS(liana.pbmc.SF7pos, "./liana.pbmc.SF7pos.rds")
liana.pbmc.SF7pos %>% glimpse()

#PBMC SF7-
pbmc.SF7.neg <- subset(pbmc.only, subset = SLAMF7 == "Neg")
#save cell cluster names as metadata to help squidpy run
pbmc.SF7.neg[["clust.names"]] <- Idents(object = pbmc.SF7.neg)
pbmc.SF7.neg <- SetIdent(pbmc.SF7.neg, value = "clust.names")
# Run liana
liana.pbmc.SF7neg <- liana_wrap(pbmc.SF7.neg,
                                method = c("sca", "natmi", "logfc", "connectome",
                                           "cellchat", "squidpy", "cellphonedb"),
                                resource = c('OmniPath'),
                                # CellChat requires normalized data
                                cellchat.params = list(.normalize=TRUE))  
saveRDS(liana.pbmc.SF7neg, "./liana.pbmc.SF7neg.rds")
liana.pbmc.SF7neg %>% glimpse()


#analyze and visualize CellPhoneDb results w/ CellTalkeR
liana.csf.SF7pos <- readRDS("./liana.csf.SF7pos.rds")
liana.csf.SF7neg <- readRDS("./liana.csf.SF7neg.rds")
glimpse(liana.csf.SF7neg)

#Now w/ dev version
devtools::install_github("https://github.com/CostaLab/CrossTalkeR@v1.3.0", ref = "devel", build_vignettes = TRUE)

#reformat column names to match new format
liana.csf.SF7pos.rename <- liana.csf.SF7pos[[7]]
liana.csf.SF7pos.rename <- rename(liana.csf.SF7pos.rename, gene_A = ligand)
liana.csf.SF7pos.rename <- rename(liana.csf.SF7pos.rename, gene_B = receptor)
liana.csf.SF7pos.rename <- rename(liana.csf.SF7pos.rename, MeanLR = lr.mean)
#add columns for type_gene_A (ligand) and type_gene_B (receptor). start by making all rows have the same value
liana.csf.SF7pos.rename <- liana.csf.SF7pos.rename %>% mutate(type_gene_A = "Ligand") %>% mutate(type_gene_B = "Receptor")


liana.csf.SF7neg.rename <- liana.csf.SF7neg[[7]]
liana.csf.SF7neg.rename <- rename(liana.csf.SF7neg.rename, gene_A = ligand)
liana.csf.SF7neg.rename <- rename(liana.csf.SF7neg.rename, gene_B = receptor)
liana.csf.SF7neg.rename <- rename(liana.csf.SF7neg.rename, MeanLR = lr.mean)
#add columns for type_gene_A (ligand) and type_gene_B (receptor). start by making all rows have the same value
liana.csf.SF7neg.rename <- liana.csf.SF7neg.rename %>% mutate(type_gene_A = "Ligand") %>% mutate(type_gene_B = "Receptor")


write_csv(liana.csf.SF7pos.rename, "./liana.csf.SF7pos.rename.csv")
write_csv(liana.csf.SF7neg.rename, "./liana.csf.SF7neg.rename.csv")

paths <- c('SF7neg' = './liana.csf.SF7neg.rename.csv',
           'SF7pos' = './liana.csf.SF7pos.rename.csv')

# Generating the report and the object
csf.ctr.rep <- generate_report(lrpaths = paths, 
                               genes=NULL,
                               out_path = "~/Documents/Research Projects/SLAMF7 neoruimmune project/scRNA_seq_interact_preditct/",
                               threshold = 50, # threshold of prune edges 0=keep all
                               out_file='2csf.ctr.rep.html') 
saveRDS(csf.ctr.rep, "./csf.ctr.rep.rds")

#Here is a nice resource on understanding pagerank scores: https://www.andreaperlato.com/graphpost/page-rank-in-network-analysis/
#basically, a higher pagerank score means that node is more important in a network. A higher score
#means that node has more in-links from other important nodes. 


##---------------1/6/2022----------------------##
#making nice tables of CSF GGI networks
#here are the colors
colours <- csf.ctr.rep@colors

csf.ctr.rep <- readRDS("./csf.ctr.rep.rds")
pgrnk <- as.data.frame(csf.ctr.rep@rankings$SF7pos_x_SF7neg$Pagerank) %>%
  rownames_to_column(var = "Cell_subsets")

ggplot(pgrnk, aes(x = reorder(Cell_subsets,csf.ctr.rep@rankings$SF7pos_x_SF7neg$Pagerank) , y = csf.ctr.rep@rankings$SF7pos_x_SF7neg$Pagerank, fill = Cell_subsets)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(size = 0.9), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        panel.background = element_rect(fill = NA)) + labs( x="Cell subsets", y="Pagerank score") +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  scale_fill_manual(values = c("#A6CEE3",    "#3F8EAA",    "#79C360",    "#B89B74",    "#E52829",    "#FDB762",    "#ED8F47", 
                               "#9471B4",    "#DDD399",    "#B15928" ))



#Now nice plot of Fisher exact test results
library(tidyverse)
fisher <- as.data.frame(csf.ctr.rep@stats$SF7pos_x_SF7neg) %>%
  dplyr::filter(csf.ctr.rep@stats$SF7pos_x_SF7neg$p <= 0.05)

ggplot(fisher, aes(x = reorder(columns_name, lodds) , y = lodds, fill = p)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(size = 0.9), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        panel.background = element_rect(fill = NA)) + labs( x="Interactions", y="Log odds ratio") +
  geom_hline(yintercept = 0, linetype="dashed", color="black") + 
  scale_fill_gradient2(high='#e6f0ff', mid='#3385ff', low='#005ce6', space='Lab')
  



#--now repeat above but for PMBCs--#
liana.pbmc.SF7pos <- readRDS("./liana.pbmc.SF7pos.rds")
liana.pbmc.SF7neg <- readRDS("./liana.pbmc.SF7neg.rds")

#reformat column names to match new format
liana.pbmc.SF7pos.rename <- liana.pbmc.SF7pos[[7]]
liana.pbmc.SF7pos.rename <- rename(liana.pbmc.SF7pos.rename, gene_A = ligand)
liana.pbmc.SF7pos.rename <- rename(liana.pbmc.SF7pos.rename, gene_B = receptor)
liana.pbmc.SF7pos.rename <- rename(liana.pbmc.SF7pos.rename, MeanLR = lr.mean)
#add columns for type_gene_A (ligand) and type_gene_B (receptor). start by making all rows have the same value
liana.pbmc.SF7pos.rename <- liana.pbmc.SF7pos.rename %>% mutate(type_gene_A = "Ligand") %>% mutate(type_gene_B = "Receptor")


liana.pbmc.SF7neg.rename <- liana.pbmc.SF7neg[[7]]
liana.pbmc.SF7neg.rename <- rename(liana.pbmc.SF7neg.rename, gene_A = ligand)
liana.pbmc.SF7neg.rename <- rename(liana.pbmc.SF7neg.rename, gene_B = receptor)
liana.pbmc.SF7neg.rename <- rename(liana.pbmc.SF7neg.rename, MeanLR = lr.mean)
#add columns for type_gene_A (ligand) and type_gene_B (receptor). start by making all rows have the same value
liana.pbmc.SF7neg.rename <- liana.pbmc.SF7neg.rename %>% mutate(type_gene_A = "Ligand") %>% mutate(type_gene_B = "Receptor")


write_csv(liana.pbmc.SF7pos.rename, "./liana.pbmc.SF7pos.rename.csv")
write_csv(liana.pbmc.SF7neg.rename, "./liana.pbmc.SF7neg.rename.csv")

paths2 <- c('SF7neg' = './liana.pbmc.SF7neg.rename.csv',
           'SF7pos' = './liana.pbmc.SF7pos.rename.csv')

# Generating the report and the object
pbmc.ctr.rep <- generate_report(lrpaths = paths2, 
                               genes=NULL,
                               out_path = "~/Documents/Research Projects/SLAMF7 neoruimmune project/scRNA_seq_interact_preditct/pbmc_ctr_out/",
                               threshold = 0, # threshold of prune edges 0=keep all
                               out_file='pbmc.ctr.rep.html')
saveRDS(pbmc.ctr.rep, "./pbmc.ctr.rep.rds")


