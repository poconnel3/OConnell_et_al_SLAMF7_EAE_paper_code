###Re-analysis of Nat Comm MS paper that used scRNA-seq. 
## Paper: https://www.nature.com/articles/s41467-019-14118-w 
#Here I am splitting cells up by location, and from each location I am splitting cells into SF7+ and SF7-
#I will then run a cell-cell interaction prediction package and compare the interactions and receptor/ligand pairs
#b/w SF7+ and SF7- cells (splitting up by location) to see cell-cell interactions enriched w/ SF7 expression

##10/10/2021
#Note: this is 10X 3' data


library(Seurat)
library(tidyverse)

# Load the all datasets
CSF.data22 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC1csf") #%>%
data22 <- CreateSeuratObject(counts = CSF.data22, project = "HC1 CSF", min.cells = 3, min.features = 200)

CSF.data23 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC1pbmc") #%>%
data23 <- CreateSeuratObject(counts = CSF.data23, project = "HC1 PBMC", min.cells = 3, min.features = 200)

CSF.data24 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC2csf") #%>%
data24 <- CreateSeuratObject(counts = CSF.data24, project = "HC2 CSF", min.cells = 3, min.features = 200)

CSF.data25 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC2pbmc") #%>%
data25 <- CreateSeuratObject(counts = CSF.data25, project = "HC2 PBMC", min.cells = 3, min.features = 200)

CSF.data26 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC3csf") #%>%
data26 <- CreateSeuratObject(counts = CSF.data26, project = "HC3 CSF", min.cells = 3, min.features = 200)

CSF.data27 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC3pbmc") #%>%
data27 <- CreateSeuratObject(counts = CSF.data27, project = "HC3 PBMC", min.cells = 3, min.features = 200)

CSF.data28 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC4csf") #%>%
data28 <- CreateSeuratObject(counts = CSF.data28, project = "HC4 CSF", min.cells = 3, min.features = 200)

CSF.data29 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC4pbmc") #%>%
data29 <- CreateSeuratObject(counts = CSF.data29, project = "HC4 PBMC", min.cells = 3, min.features = 200)

CSF.data30 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC5csf") #%>%
data30 <- CreateSeuratObject(counts = CSF.data30, project = "HC5 CSF", min.cells = 3, min.features = 200)

CSF.data31 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC5pbmc") #%>%
data31 <- CreateSeuratObject(counts = CSF.data31, project = "HC5 PBMC", min.cells = 3, min.features = 200)

CSF.data32 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/HC6csf") #%>%
data32 <- CreateSeuratObject(counts = CSF.data32, project = "HC6 CSF", min.cells = 3, min.features = 200)

CSF.data33 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS1csf") #%>%
data33 <- CreateSeuratObject(counts = CSF.data33, project = "MS1 CSF", min.cells = 3, min.features = 200)

CSF.data34 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS1pbmc") #%>%
data34 <- CreateSeuratObject(counts = CSF.data34, project = "MS1 PBMC", min.cells = 3, min.features = 200)

CSF.data35 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS2csf") #%>% 
data35 <- CreateSeuratObject(counts = CSF.data35, project = "MS2 CSF", min.cells = 3, min.features = 200)

CSF.data36 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS2pbmc") #%>%
data36 <- CreateSeuratObject(counts = CSF.data36, project = "MS2 PBMC", min.cells = 3, min.features = 200)

CSF.data37 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS3csf") #%>%
data37 <- CreateSeuratObject(counts = CSF.data37, project = "MS3 CSF", min.cells = 3, min.features = 200)

CSF.data38 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS3pbmc") #%>%
data38 <- CreateSeuratObject(counts = CSF.data38, project = "MS3 PBMC", min.cells = 3, min.features = 200)

CSF.data39 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS4csf") #%>%
data39 <- CreateSeuratObject(counts = CSF.data39, project = "MS4 CSF", min.cells = 3, min.features = 200)

CSF.data40 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS4pbmc") #%>%
data40 <- CreateSeuratObject(counts = CSF.data40, project = "MS4 PBMC", min.cells = 3, min.features = 200)

CSF.data41 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS5csf") #%>%
data41 <- CreateSeuratObject(counts = CSF.data41, project = "MS5 CSF", min.cells = 3, min.features = 200)

CSF.data42 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS5pbmc") #%>%
data42 <- CreateSeuratObject(counts = CSF.data42, project = "MS5 PBMC", min.cells = 3, min.features = 200)

CSF.data43 <- Read10X(data.dir ="~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/Nat_Comm_sc_RNA_seq_reanalysis/GSE138266_RAW/MS6csf") #%>%
data43 <- CreateSeuratObject(counts = CSF.data43, project = "MS6 CSF", min.cells = 3, min.features = 200)

#now need to merge all together (do not integrate)
CSF_all <- merge(data22, y = c(data23, data24, data25, data26, data27, data28, data29, data30, data31, data32, data33, data34, data35, data36, data37, data38, data39, data40, data41, data42, data43), 
                 add.cell.ids = c("HC1 CSF", "HC1 PBMC", "HC2 CSF", "HC2 PBMC", "HC3 CSF", "HC3 PBMC", "HC4 CSF", "HC4 PBMC", "HC5 CSF", "HC5 PBMC", "HC6 CSF", "MS1 CSF", "MS1 PBMC", "MS2 CSF", "MS2 PBMC", "MS3 CSF", "MS3 PBMC", "MS4 CSF", "MS4 PBMC", "MS5 CSF", "MS5 PBMC", "MS6 CSF"), project = "NatComm")
CSF_all 

saveRDS(CSF_all, file = "./csf_all.rds")

#remove all non-Seurat object data files 
rm(list=ls(pattern="CSF"))

#remove all single-Seurat objects
rm(list=ls(pattern="data"))

CSF_all <- readRDS(file = "./csf_all.rds")

head(colnames(CSF_all))

unique(sapply(X = strsplit(colnames(CSF_all), split = "_"), FUN = "[", 1))

table(CSF_all$orig.ident)

#calculate percent mitochondrial genes and add it to new column in metadata
CSF_all[["percent.mt"]] <- PercentageFeatureSet(CSF_all, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(CSF_all@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(CSF_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(CSF_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CSF_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(CSF_all, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2
plot3

#Only take data that has cells that have more then 200 and less then 6000 transcripts and less the 5% mito genes
csf_all <- subset(CSF_all, subset = nFeature_RNA > 200 & nFeature_RNA < 5800 & percent.mt < 5)

# run sctransform (normalizes, scales, finds variable features, and removes unwanted variability) (really slow w/ this many cells)
csf_all <- SCTransform(csf_all, vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS(csf_all, file = "./csf_all.rds")

# visualize and cluster
set.seed(4646)
csf_all <- RunPCA(csf_all, verbose = FALSE)
csf_all <- RunUMAP(csf_all, dims = 1:30, verbose = FALSE)

csf_all <- FindNeighbors(csf_all, dims = 1:30, verbose = FALSE)
csf_all <- FindClusters(csf_all, verbose = FALSE)
DimPlot(csf_all, label = TRUE) + NoLegend()

# canonical marker genes as violin plots.
VlnPlot(csf_all, features = c("CD8A", "NCAM1", "CD4", "CD19", "CD3D", "CD3E", "SLAMF7", "TRAC", "CD8B", "FOXP3", "TRDC", "GNLY", "FCGR3A", "SELL", "CD79A", "IGHD", "LYZ", "S100A8", "CD14", "TCF4", "GNG11"), 
        pt.size = 0.2, ncol = 5)

# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(csf_all, features = c("CD8A", "NCAM1", "CD4", "CD19", "CD3D", "CD3E", "SLAMF7", "TRAC", "CD8B", "FOXP3", "TRDC", "GNLY", "FCGR3A", "SELL", "CD79A", "IGHD", "LYZ", "S100A8", "CD14", "TCF4", "GNG11"), pt.size = 0.2, 
            ncol = 5)

table(Idents(csf_all))

###need to see top genes in all clusters and remove ones w/ low gene number 
##and ones that have only mitochondrial and ribisomal transcripts 
VlnPlot(csf_all, features = c("nFeature_RNA"), x.lab.rot = TRUE)
VlnPlot(csf_all, features = c("percent.mt"), x.lab.rot = TRUE)

##remove cluster 23, 24, low quality cells. 
csf_all <- subset(csf_all, idents = c("23", "24"), invert = TRUE)

#Pre-process again
csf1 <- SCTransform(csf_all, vars.to.regress = "percent.mt", verbose = FALSE)

csf1 <- RunPCA(csf1, verbose = FALSE)
csf1 <- RunUMAP(csf1, dims = 1:25, verbose = FALSE)

csf1 <- FindNeighbors(csf1, dims = 1:25, verbose = FALSE)
csf1 <- FindClusters(csf1, verbose = FALSE)
DimPlot(csf1, label = TRUE) + NoLegend()

saveRDS(csf1, file = "./csf1.rds")

#Now combine like clusters and manually annotate same as Nat Comm paper
DimPlot(csf1, label = TRUE) + NoLegend()

#Find top De genes in each cluster (placed out of order)
#Identify cluster markers
csf.markers <- FindAllMarkers(csf1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
csf.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(csf.markers %>% group_by(cluster) %>% top_n(10, avg_logFC), 'csf.markers.tsv', sep='\t')
#Heat map of top DE genes per cluster
top10 <- csf.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(csf1, features = top10$gene) + NoLegend()
saveRDS(csf.markers, file = "./csf.marker.rds")

csf1.markers <- readRDS(file = "./csf.marker.rds")

#Feature plot of cannonical markers from oroginal paper on my map
FeaturePlot(csf1, features = c("CD8B", "CD4", "CD3E", "IL7R", "CCR7", "FOXP3", "TRDC", "GNLY", "FCGR3A", "SELL", "TRAC", "CD27", "CD38", "TCF4"), pt.size = 0.2, 
            ncol = 4)

#remove some more junk (annotated below and re-SCT and plot)
csf1 <- subset(csf1, idents = c("3", "21", "29"), invert = TRUE)

#name clusters
new.cluster.ids <- c("0", "CD8 T cells", "2", "Apoptotic", "4", "Granulocytes", "6", "7", "NK cells 1", "CD8 T cells", "10", "B cells 1", "Monocytes 2", "13", "NK cells 2", 
                     "15", "Granulocytes", "CD8 T cells", "18", "mDC2", "B cells 2", "Platelets", "Tregs", "23", "Monocytes 1", "Plasma cells", "NK cells 1", "pDCs", "mDC1", "Junk")

#Pre-process again
csf2 <- SCTransform(csf1, vars.to.regress = "percent.mt", verbose = FALSE)

csf2 <- RunPCA(csf2, verbose = FALSE)
csf2 <- RunUMAP(csf2, dims = 1:20, verbose = FALSE)

csf2 <- FindNeighbors(csf2, dims = 1:20, verbose = FALSE)
csf2 <- FindClusters(csf2, verbose = FALSE)
DimPlot(csf2, label = TRUE) #+ NoLegend()

saveRDS(csf2, file = "./csf2.rds")
#-----------------------------------------------------------#


#-----Restart from old pre-processed data------#
#Note: I analyzed this data set previously and will just load in a seurat object w/ clean cells and annotations. 

#load Seurat object
csf2 <- readRDS("./csf2.rds")
head(csf2@meta.data)
#----------------Now we starting analyzing the annotated data set----------------------#

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

##-------------------Subset out B cells
all_B <- subset(csf2, idents = c("B cells"))
#compare SF7 expression by location
VlnPlot(all_B, features = c("SLAMF7"), slot = "counts", assay = "RNA", split.by = "location")
VlnPlot(all_B, features = c("SLAMF7"), slot = "counts", assay = "RNA", split.by = "condition")
table(Idents(all_B), all_B$condition)

#mark all B cells as SLAMF7 pos or neg
# Subset on the expression level of a gene/feature
sf7_pos <- subset(x = all_B, subset = SLAMF7 > 0.5) #positive cutoff is heuristic based on visual inspection of violin plots
pos_names <- rownames(sf7_pos@meta.data)
# Mutate a column in original  metadata
all_B$barcode <- rownames(all_B@meta.data)
all_B@meta.data <- all_B@meta.data %>% 
  mutate(SLAMF7 = ifelse((all_B$barcode %in% pos_names), "Pos",  "Neg"))

#see how many SF7+ B cells b/w conditions (1/22/2022)
table(all_B$SLAMF7, all_B$condition)
table(all_B$SLAMF7, all_B$location, all_B$condition)
prop.table(table(all_B$SLAMF7, all_B$condition), margin = 2)




#---------Predict cell-cell interactions for all SF7+ cells and all SF7- cells---------#
library(CrossTalkeR)
library(liana)  #this package is a PITA to load. I can tell using this will be fun....
# actually runs nicely, just takes a while. 

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


#analyze and visualize CellPhoneDb results w/ CellTalkeR
liana.csf.SF7pos <- readRDS("./liana.csf.SF7pos.rds")
liana.csf.SF7neg <- readRDS("./liana.csf.SF7neg.rds")

#Must install dev version of CrossTalkeR
devtools::install_github("https://github.com/CostaLab/CrossTalkeR@v1.3.0", ref = "devel", build_vignettes = TRUE)

#reformat column names to match new format
liana.csf.SF7pos.rename <- liana.csf.SF7pos[[7]] #This selects only CellPhoneDB output
liana.csf.SF7pos.rename <- rename(liana.csf.SF7pos.rename, gene_A = ligand)
liana.csf.SF7pos.rename <- rename(liana.csf.SF7pos.rename, gene_B = receptor)
liana.csf.SF7pos.rename <- rename(liana.csf.SF7pos.rename, MeanLR = lr.mean)
#add columns for type_gene_A (ligand) and type_gene_B (receptor). start by making all rows have the same value
liana.csf.SF7pos.rename <- liana.csf.SF7pos.rename %>% mutate(type_gene_A = "Ligand") %>% mutate(type_gene_B = "Receptor")


liana.csf.SF7neg.rename <- liana.csf.SF7neg[[7]] #This selects only CellPhoneDB output
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

#Here is a nice resource on underatanding pagerank scores: https://www.andreaperlato.com/graphpost/page-rank-in-network-analysis/
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
  
