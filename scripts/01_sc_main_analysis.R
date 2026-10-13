rm(list = ls())
options(stringsAsFactors = F)
gc()
#  
save.image("recent.RData")


af[["RNA"]] <- as(af[["RNA"]], "Assay")#5 to 4

af[["RNA"]] <- as(object = af[["RNA"]], Class = "Assay5")#4 to 5

table(af$Cell_Anno)
cwd_color<-c("royalblue",  "indianred")
region_color<-c('#FB9D32','#DC51AD','#9E6FF3','#BBE379','#2b2b2b','#27B0D9')

cell_color<-c("#D10D0F","#00A6D0", "#FFB600","#42A24E",'#0F7240','#507BA8',"#354C85",
           "#83C5BD", '#5d2b84',"#BBB8D7"   )
hip_color<-c("#ffabc3", "#d9da00", "#827cb9", "#00bea8", "#cc8b00","#a9ceff")
           
DG_color<-c( "#d9da00", "#827cb9" ,'#f3758a','#7cec37' )



library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(ggpubr)
library(Biobase)
library(Matrix)
library(RcppML) 
library(magick)
library(GeneNMF)#R>=4.3
library(scRNAtoolVis)
library(Nebulosa)
library(reshape2)
library(ggnetwork)
library(ggimage)
library(mascarade)
library(UCell)
library(ggrastr)
library(qs)



setwd("./")
getwd()



pdf(file = "01.vlnplot.pdf",width = 8,height = 5)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA"), ncol = 4)
dev.off()

plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#
pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot2+plot3+plot_layout(ncol = 3)      #plot_layout
dev.off()

plot3
# Remove mitochondrial and ribosomal genes
af <- SCTransform(af, assay = "Xenium", verbose = FALSE)

af <- Seurat::RunPCA(af)


pct <- af[["pca"]]@stdev / sum(af[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(af)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

af <- Seurat::RunPCA(af,npcs = 16)

af <- RunHarmony(af, group.by.vars = "Type")

namecon="CWD"     
nametreat="IWD"


#
conaf=c("CWD289","CWD291","CWD361")
#
af$oldType=af$Type  
af$Type=factor(ifelse(af$Type %in% conaf,namecon,nametreat))


Idents(af)=af$Type

table(af$oldType)


stage_vector <- rep(NA, length(af.ST$group))


stage_vector[af.ST$group %in% c("CWD1", "CWD2",
                                'CWD3')] <- 'CWD'
stage_vector[af.ST$group %in% c("IWD1", "IWD2",
                                'IWD3')] <- 'IWD'
#
af.ST$Type <- stage_vector

af <- FindNeighbors(af, dims = 1:16,reduction = "harmony")

af <- FindClusters(af, resolution = 0.6)

# 
head(Idents(af), 5)

# 
head(af@meta.data)
table(af@meta.data$seurat_clusters)
# 

af.markers <- FindAllMarkers(af.t_LEC, only.pos = F, min.pct = 0.25, 
                             logfc.threshold = 0.1)
write.csv(af.markers,file = "tumor_LEC_cluster_markers.csv")
af.markers_cluster <- af.markers[af.markers$cluster == 8, ]

myMarkers=c("COL4A1",'TIMP2','APP','CCN1','PECAM1','SYNM','SPTBN1')
pdf(file = "cluster.pdf",width =4.5,height = 4)
jjVolcano(diffData = af.markers, aesCol = c('#42AFFC','tomato'),#topGeneN=3,
          myMarkers=myMarkers,flip = T,legend.position = c(0.9,0.2),#polar=T,
          tile.col = colors_LEC)##cluster
dev.off()


af <- RunUMAP(af, dims = 1:16, reduction = "pca",min.dist = 0.7)

# 
pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = F, label.size = 2.5,pt.size = 1,
       group.by = "Type",
        cols= colors8)+#ggsci::scale_color_igv()+
  theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

af$seurat_clusters<-Idents(af)


af <- FindClusters(af, resolution = 0.3)
table(af@active.ident)

ann.ids <- c("Neurons_EX",  #
             "remove",  #1
             "remove",  #2
             "Oligo", #3
             "Neurons_IN", #4
             "Astro",#5
             "Neurons_EX",#6
             "Neurons_EX",#7
             "remove",#8
             "Neurons_EX",#9
             "remove",#10
             "Neurons_IN",#11
             "Neurons_IN",#12
             "OPC",#13
             "Micro",#14
             "Ependemal",#15
             "Neurons_EX",#16
             "Neurons_EX",#17
             "Neurons_EX",#18
             "remove",#19
             "remove",#20
             "Oligo",#21
             "Neurons_EX",#22
             "remove",#23
             "remove",#24
             "remove",#25
             "remove",#26
             "remove"#25
)

Idents(af)

afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens

af$Celltype=Idents(af)
########
pdf(file = "anno.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = F, label.size = 2.5,pt.size = 0.0001,
        group.by = "Celltype",
        cols= colors1)+#ggsci::scale_color_igv()+
  theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
# 
emb <- Embeddings(af,reduction = 'umap') %>% #   
  as.data.frame() %>% mutate(celltype=af@meta.data$Celltype)
maskTable <- generateMask(dims=emb[,1:2],
                          clusters=emb$celltype,#######
                          smoothSigma = 0.015,  
                          minDensity = 0.1,      
                          gridSize = 200,      
                          type = 'partition')

# 
pdf(file = "cluster_UMAP.pdf",width =7,height = 5.5)
DimPlot(af, label = F, pt.size = 1)&  geom_path(data=maskTable[cluster=='NK cell'],####
                                                aes(x=umap_1, y=umap_2,group=group), 
                                                linetype=6,size=1)&
  scale_color_manual(values =alpha(c("#eb6841","#edc951","#4f372d","#00a0b0",
                                     "#85218B", "#849271", "#CF9546", "#C67052", "#7A989A","#208A42","#C1AE8D"),0.5))
dev.off()


#

source(file = "./R/ks_scAverExp.R")
source(file = "./R/ks_scDotplot.R")

featureSets <- list(
                    "Neurons_EX" = c('Rbfox3',"Slc17a7", "Slc17a6"),
                    "Neurons_IN" = c("Slc32a1", "Gad1"),
                    "Astro" = c("Slc1a3",'Atp1b2'),
                    "Oligo" = c("Mag", "Mog"),#("Mbp", "Hapln2"
                    "OPC" = c("Pdgfra","Cspg4"),
                    "Micro" = c("Itgam",'Fcrls'),
                    "FBs" = c("Dcn",'Col1a1'),
                    "EC" = c("Pecam1"),
                    "Epen" = c("Dnah11"),
                    "Peri" = c("Pdgfrb")
              
                    )

Exp_scRNA <- ks_scAverExp(obj = sce,#seurat or ArchR proj
                          features = featureSets,
                          CellTypes = unique(sce$celltype))

scRNA_dotplot = ks_scDotplot(obj = sce,#seurat or ArchR proj
                             features = featureSets,
                             CellTypes = unique(sce$celltype),
                             avgPctMat= Exp_scRNA,
                             pal = cell_color)


af.markers <- FindAllMarkers(afd, only.pos = T, min.pct = 0.25, logfc.threshold = 0.75)

# get top 10 genes
topaf.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


# plot
pdf(file = "09-cell_marker.hetmap.pdf",width =7,height = 7)
AverageHeatmap(object=afd,
               markerGene = topaf.markers$gene,slot = "data",assays='RNA',
               clusterAnnoName =F,annoCol=T,myanCol= colors2[1:2],# htRange=c(-0.5, 0, 2),
               htCol =c('#435B95','white','#C3400C')
               )






afgenes=c("Rbfox3",#neu
          "Slc17a7", "Slc17a6",#ex
         "Slc32a1", "Gad1",#in
          "Slc1a3",'Atp1b2',#AST
          "Pdgfra","Cspg4",#OPC
          "Mog",'Mag',#OLI
         "Itgam", "C1qa",'Fcrls',#MICRO
          "Foxj1", "Pdgfrb",#EPN
           'Col1a1','Pecam1'          )


XX=c('Dlg4','Snap25','Vamp2','Syp','Grin1','Pkm','Ckb','Vdac2','Mfn2','Cxcl12','Ighm','Bdnf' )

XX = c("Sdhc",'Mpc2','Dnm1l','Atp5f1','Atp2b3','Nlgn2','Kcnc2','Pde1b','Ptpn5')
XX = c("Gpr124",'Tek','Ocln','Cdh5','Pecam1','Itgb1','Vdac2','Mpc2','Pygb','Ndufs1')

XX = c('Nr3c1','Tsc22d3','Foxo1','Pdk4','Slc2a1','Pcx','Pkm')
XX = c('Zbtb16','Il1rap','Cx3cr1','Il10ra','Tgfbr1','Tgfbr2')
XX=c('Foxo1','Slc16a1','Pkm','Pcx','Slc2a1','Akt1','Akt2')
XX=c('Slc16a3','Slc16a1')


# 
Idents(af)<-af$Type
pdf(file = "spatial_marker_dotplot.pdf",width =7,height = 3)
DotPlot(af, features = XX)+theme_bw()+#

theme(panel.grid = element_blank(),
           
axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+#
  
scale_color_gradientn(values = seq(0,1,0.2),
                      colours = c('#64A2CC',"white",'#E07B56',"#731929"))+#
 # c('#64A2CC',"white",'#E07B56',"#731929")
labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+coord_flip()
dev.off()
  ###################

################################
Idents(af)<-af$Celltype





#################
#
featureSets <- list(
  "Neurons_EX" = c('Rbfox3',"Slc17a7", "Slc17a6"),
  "Neurons_IN" = c("Slc32a1", "Gad1"),
  "Astro" = c("Slc1a3",'Atp1b2'),
  "Oligo" = c("Mag", "Mog"),#("Mbp", "Hapln2"
  "OPC" = c("Pdgfra","Cspg4"),
  "Micro" = c("Itgam",'Fcrls'),
  "FBs" = c("Dcn",'Col1a1'),
  "EC" = c("Pecam1"),
  "Epen" = c("Dnah11"),
  "Peri" = c("Pdgfrb")
  
)
#afgenes=read.table("stroke_geneFA.txt",header = F,sep = "\t")[,1]
pdf(file = "spatial_Rbfox3.pdf",width =6,height = 5)
plot_density(af, features = c("Foxo1"),
             pal = 'viridis', raster = T, size = 0.8) &
  theme_blank()&#
  theme(legend.frame = element_rect(colour = "black"),#
        legend.ticks = element_line(colour = "black", linewidth  = 0),#
        legend.key.width = unit(0.3, "cm"),#
        legend.key.height = unit(0.8, "cm"),#
        legend.title = element_text(color = 'black', face = "bold", size=8))#
dev.off()

table(af.ST$Cell_Anno)


##### 
source("./R/Singlecellratio_plotstat.R")


af<-af.region

af<- subset(af.ST, Type %in% c("CWD"))

af$celltype<-af$Cell_Anno
af$orig.ident<-af$Hip#######

Singlecellratio_plotstat(af, group_by = "Hip",###
                         meta.include = c("Hip","orig.ident"),#
                         color_by = 'cell.type'###
)+scale_fill_manual(values=cell_color)#

ggsave("ration.pdf",width = 5,height = 4.5)


source(file = "./R/ks_cellration_plot.R")
library(jjAnno)
af.LEC$LEC


ks_scFraction(group = af$celltype,
              celltype = af$Hip,pal = hip_color,#
              colorLab = cell_color , Total = T)#
ggsave("niche_ration.pdf",width = 6.5,height = 4.5)


af$orig.ident<-afe$group###
af$celltype<-afe$Cell_Anno  ####
my_comparisons<-list(c("IWD","CWD"))
Singlecellratio_plotstat(afe, group_by = "Type",
                         custom_fill_colors = c("royalblue",  "indianred"),
                         meta.include = c("Type","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Type',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =4)
ggsave("niche_2_3 cellratio_diff.pdf",width = 10,height = 5.5)






