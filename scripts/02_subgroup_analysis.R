

#group diff
source(file = "./R/vnplot.R")
gene_sig <- 'Pomc'
#cell diff

comparisons <- list()
comp=combn(names(table(af$Type)),2)
af$Type
#names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(afe,gene_signature = gene_sig, file_name = "Slc2a1", test_sign = comparisons,pta=0,
     cols=colors6,label="p.signif",group="Type",widplot=4,heiplot=4,ak=1)

af$fat_acid
######
source("./R/ks_VlnExp.R")


my_comparisons3 <- list(c("IWD", "CWD"))

ks_VlnExp(object = afe, group="Type",group_order=c("IWD",'CWD'),
          features='Pomc',cols=c(  "indianred","royalblue"),
          comparisons=c(list(c("IWD", "CWD"))
          ))


ggsave("Lepr .pdf",width = 4.5,height = 6.5)


table(afd$Hip)
rm(af)
afd[["RNA"]] <- as(afd[["RNA"]], "Assay")
table(af.ST$Cell_Anno)
af.ST<- subset(af.ST, Cell_Anno %in% c("Neurons_EX"))


#DefaultAssay(afd) <- "RNA"

afd <- NormalizeData(afd, normalization.method = "LogNormalize", scale.factor = 10000)
afd <- FindVariableFeatures(afd, selection.method = "vst", nfeatures = 2000)
scale.genes<-VariableFeatures(afd)
afd <-ScaleData(afd,features = scale.genes)
afd <- Seurat::RunPCA(afd, features = VariableFeatures(object = afd))
#afd <- RunHarmony(afd, group.by.vars = "Type")#### No need for PCA anymore

ElbowPlot(afd,reduction = "pca")
afd=FindNeighbors(afd, dims = 1:12, reduction = "pca")
afd <- RunUMAP(afd, reduction = "pca", dims = 1:12)

afd <- FindClusters(afd, resolution = 0.3)

table(afd@active.ident)
afd$sub.cluster<-Idents(afd)


afgenes=c('Pdgfrb','Cspg4','Ccdc113','Dnah11','Pecam1','Cldn5','Dcn','Col1a1')

DotPlot(afd, features = afgenes)+theme_bw()+# Remove background, rotate image
  
  theme(panel.grid = element_blank(),
        
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+#
  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#64A2CC',"white",'#E07B56',"#731929"))+#
  # c('#64A2CC',"white",'#E07B56',"#731929")
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))#+coord_flip()


ann.ids <- c("EC",  #0Injured Proximal tubular cell  not a cluster, it is my own index, cluster+1
             "FBs",  #1Injured Proximal tubular cell
             "Peri",  #2
             "Epen", #3
             "FBs", #4
             "Peri",#5
             'Epen',#6
             'remove',#7
             'remove'#8
             
)
Idents(afd)

afidens=mapvalues(Idents(afd), from = levels(Idents(afd)), to = ann.ids)
Idents(afd)=afidens
afd$POSTN=Idents(afd)
pdf(file = "subcluster_anno_UMAP.pdf",width =7,height = 5.5)
DimPlot(afd, reduction = "umap", label = F, label.size = 2.5,
       group.by = "Hip",
        pt.size = 1,cols= colors5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
table(Idents(afd))

#################
Idents(af, cells = colnames(afd)) <- Idents(afd)


table(Idents(af))

scRNA<- subset(scRNA, subset = Celltype != "Peri")

#Idents(af.1)<-"Proximal tubule cell"

af.ST
af<- subset(af1, idents = selected_idents)

af1$afd_subcluster<-Idents(af1)

Idents(af.ST, cells = colnames(af.stroma)) <- Idents(af.stroma)

table(af.ST$doamin2)
