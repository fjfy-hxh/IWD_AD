
write.csv(CellChatDB.mouse$interaction, "cellchat_mouse_interactions.csv", row.names = FALSE)
write.csv(CellChatDB.mouse$geneInfo, "cellchat_mouse_geneinfo.csv", row.names = FALSE)

object.list = qread("cellchat_object.qs")

af.chat = readRDS("af.chat.rds")
af.chat<-af
af.chat$cellType <- Idents(af.chat)

af.chat$Celltype<-af.chat$Cell_Anno
aflist=SplitObject(af.chat,split.by = "Type") 
table(af.chat$Celltype)


######
af.chat$Celltype <- gsub("Neurons_EX", "Neurons EX", af.chat$Celltype)
af.chat$Celltype <- gsub("Neurons_IN", "Neurons IN", af.chat$Celltype)
##################################
cellchat.Control = createCellChat(object = aflist$CWD, 
                                  group.by = "Celltype")

#########################
#new_levels <- c('B cells','CD4+ T cell','CD8+ T cell','Endo','Epi','FBs','MCs','NK cells',
#                        'PVL')

# 
#cellchat.Control@idents <- factor(cellchat.Control@idents, levels = new_levels)

# 
levels(cellchat.Control@idents) #

groupSize <- as.numeric(table(cellchat.Control@idents)) #
#
CellChatDB = CellChatDB.mouse #
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = CellChatDB

cellchat.Control <- subsetData(cellchat.Control)
cellchat.Control <- identifyOverExpressedGenes(cellchat.Control)
cellchat.Control <- identifyOverExpressedInteractions(cellchat.Control)

cellchat.Control <- computeCommunProb(object=cellchat.Control,raw.use = TRUE)

cellchat.Control <- computeCommunProbPathway(cellchat.Control)

cellchat.Control <- aggregateNet(cellchat.Control)
###################################################################################
cellchat.cABMR = createCellChat(object = aflist$IWD, group.by = "Celltype")

##
#new_levels <- c("B cells",  "Mono/mac", "DCs", "Mast", "Neu", 
#                "T cells", "NK cells", "Endo", "Epi", "Fibro")
# 
# 
levels(cellchat.cABMR@idents)

groupSize <- as.numeric(table(cellchat.cABMR@idents)) #
#
CellChatDB = CellChatDB.mouse #
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = CellChatDB
cellchat.cABMR@DB = CellChatDB.use

cellchat.cABMR <- subsetData(cellchat.cABMR)
cellchat.cABMR <- identifyOverExpressedGenes(cellchat.cABMR)
cellchat.cABMR <- identifyOverExpressedInteractions(cellchat.cABMR)

cellchat.cABMR <- computeCommunProb(object=cellchat.cABMR,raw.use = TRUE)
#
cellchat.cABMR <- computeCommunProbPathway(cellchat.cABMR)
#
cellchat.cABMR <- aggregateNet(cellchat.cABMR)



#############
object.list <- list('CWD' = cellchat.Control, 'IWD' = cellchat.cABMR)
object.list1 <- list(
  `FABP1- Tumor` = cellchat_object$`FABP1- Tumor`,
  `FABP1+ Tumor` = cellchat_object$`FABP1+ Tumor`
)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

qsave(object.list, file = "cellchat_object.qs")
#########################################
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = c("#6bb9d2","#d55640"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = c("#6bb9d2","#d55640"), measure = "weight")
pdf(file = "13.cellchat_compareInteractions.pdf",width =7,height = 3)
gg1 + gg2
dev.off()

###########################
pdf(file = "13.netVisual_circle_CON_count.pdf",width =6,height = 6)
netVisual_circle(cellchat.Control@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions in NAT")
dev.off()

pdf(file = "13.netVisual_circle4_CON_weight.pdf",width =6,height = 6)
netVisual_circle(cellchat.Control@net$weight, color.use=colors3,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights in NAT")
dev.off()

pdf(file = "13.netVisual_circle4_TUMOR_count.pdf",width =6,height = 6)
netVisual_circle(cellchat.cABMR@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions in TUMOR")
dev.off()
pdf(file = "13.netVisual_circle4_TUMOR_weight.pdf",width =6,height = 6)
netVisual_circle(cellchat.cABMR@net$weight, color.use=colors3,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights in TUMOR")
dev.off()



mat <- cellchat.Control@net$weight
pdf(file = "1.pdf",width =10,height = 10)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,  weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

mat1 <- cellchat.cABMR@net$weight
pdf(file = "2.pdf",width =10,height = 10)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat1)) {
  mat2 <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
  mat2[i, ] <- mat1[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat1),
                   title.name = rownames(mat1)[i])
}
dev.off()


pdf(file = "13.netVisual_diffInteraction.pdf",width =6,height = 6)
par(mfrow = c(1,1), xpd=TRUE)

netVisual_diffInteraction(cellchat, weight.scale = T,
                          vertex.size.max = 20,
                          edge.width.max = 20,arrow.width = 2,
                          color.use =c( '#f598b3', '#282568','#73c9e4',  '#156434','#9D0709') ,
                          color.edge = c("#F3AE63","#73558B"), comparison = c(2, 1),
                        #  targets.use = 5,
                       #   sources.use = NULL
                        #  title = "primary VS metastasis"
                        )
dev.off()

############################
cellchat_col <- c(
  "B cells" = "#4991C1",
  "DCs" = "#D51F26",
  'Endo'="#BDA7CB",
  'Epi'="#435B95",
  'Fibro'="#C6307C",
  'Mast'='#8D75AF',
  "MC" = "#C0937E",
  "Mono/mac" = "#7AA15E",
  'Neu'='#EE934E',
  'NK cells'="#E45D61",
  'T cells'="#67A59B"
)

pdf(file = "13.netVisual_heatmap1.pdf",width =12,height = 6)
a1=netVisual_heatmap(cellchat, measure=c("weight"),
                     color.use = c( "#FFB600",  "#83C5BD" ,"#5d2b84", "#354C85",
                                    "#507BA8","#D10D0F", "#00A6D0",
                                    "#42A24E", "#0F7240" ,
                                   "#BBB8D7"),
                     color.heatmap = c('#64A2CC',"#860808")
)
a2=netVisual_heatmap(cellchat,
                    # color.use = cellchat_col,
                     color.heatmap = c('#64A2CC',"#860808"))
a1+a2
dev.off()

pdf(file = "13.rankNet.pdf",width =5,height = 4.5)
rankNet(cellchat, mode = "comparison", stacked = T, 
        sources.use = 10, #
        targets.use = 2, #
        color.use = c("indianred","royalblue"),
        do.stat = TRUE, comparison = c(2,1))
dev.off()
#########
levels(cellchat@idents$joint) #
pdf(file = "13.netVisual_bubble.pdf",width =9,height = 6)
netVisual_bubble(cellchat, 
                 #sources.use = 2, # sender cells
                 targets.use = 2, # receiver cells
                 comparison = c(2, 1), angle.x = 45)
dev.off()
############################################
############################################################################################
source(file = "./R/ks_CELLCHAT_groupCells_bubble.R")

net <- list()
table(afc$Type)

for (i in 1:length(object.list)) {
  
  net1 <- subsetCommunication(object.list[[i]])
  net[[i]] <- net1
  names(net)[i] <- group_name[i]
}

# 
find_interaction_by_ligand <- function(data, ligand_name) {
  # 
  if (ligand_name %in% data$ligand) {
    # 
    result <- data$interaction_name[data$ligand == ligand_name]
    return(result)
  } else {
    # 
    return(paste(ligand_name, "not found in data$ligand"))
  }
}
# 

find_interaction_by_ligand(my_iCAF_data, "APP")
View(my_iCAF_data)
my_iCAF_data$interaction_name
my_iCAF_data$interaction_name
c('#f598b3', '#282568','#73c9e4',  '#156434','#9D0709')


object.list <- list('cellchat' = cellchat.Control, 'a WD' = cellchat.Control)
cellchat.net$interaction_name
#######
ks_CC_groupCells_bubble(object_list=object.list,#####
                        group_name = c("IWD",'CWD'), #
                        celltypes=c("Astro","EC"),#
                        LR_int =c('ANGPTL4_CDH5','VEGFA_VEGFR2'
                                  ),#c("VEGFA_VEGFR1R2",'VEGFA_VEGFR1','APP - CD74')
                        cell_cols = c('#cce4b2',"#FCCC3A", "#282568"),##
                        cols=c('#64A2CC','white','#E07B56',"#731929" ),
                        group_cols = c("#6bb9d2","#d55640"))

object.list$`HR+ HER2- BC`@meta
##################################################################################
source(file = "./R/cellchat_chord_plot.R")

cellchat
cellchat.Control<-object.list$CWD
net<-subsetCommunication(cellchat.Control)
View(cellchat.net)
cellchat.net<- subsetCommunication(cellchat.Control,
                                        #sources.use = "Neurons EX",  #
                                        targets.use = "EC"
                                   ) 

cellchat.net$source <-gsub("_", " ", cellchat.net$source)  # 
cellchat.net$target <-gsub("_", " ", cellchat.net$target)
unique(cellchat.net$interaction_name)

pdf(file = "ligand_receptor_chord.pdf",width =12,height = 10)
plot_cellchat_LR(cellchat_net = cellchat.net,#####
                 celltype_inter = c("EC","Astro"),#
                 celltype_color = c('#1F1C53'='EC',
                                    '#FCCC3A'="Astro"),
                 ligand_col="#FB8072",# ligand color
                 receptor_col="#1988B0",# receptor color
                 text_size=0.8,
                 Group="",### used for plotting in the past
                 legend_in_plot=T)
dev.off()

cellchat.filtered



