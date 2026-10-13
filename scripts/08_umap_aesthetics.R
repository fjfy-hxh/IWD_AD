#http://www.360doc.com/content/24/0430/20/82995563_1122020662.shtml


#############################
source("./R/umap_plot_functions.R")


###Small arrow aesthetic plot####
#plotData$UMAP_1
p <- plot_umap_with_unchull(
  seurat_obj = af,
  cluster_col = "Celltype",
  colors = colors3,
  point_size = 0.03#0.2-0.3
)

print(p)

#pdf(file = "Spatial major cluster UMAP.pdf",width =9,height = 6)
pdf(file = "LEC C1_C3 UMAP.pdf",width =6,height = 4)
p
dev.off()

###Hazy UMAP
p1 <- plotUMAP(
  seurat_obj = af.epi,
  reduction = "umap",
  cluster_col = "Type",
  color_palette = colors3,
  axis_line_size = 0.3,      # axis line width
  legend_point_size = 3,     # legend point size
  legend_text_size = 5,      # legend text size
  legend_ncol = 1 ,          # number of legend columns
  legend_alpha = 0.7 ,
  point_size = 0.2
)
p1
ggsave(filename ="Episub UMAP .png",height = 3,width = 4,plot = p1,dpi = 300,bg="white")










# Highlight only key cell populations


emb <- Embeddings(af,reduction = 'umap') %>% # extract coordinate data    
  as.data.frame() %>% mutate(celltype=af@meta.data$Celltype)
maskTable <- generateMask(dims=emb[,1:2],
                          clusters=emb$celltype,#######do not change this
                          smoothSigma = 0.015,  
                          minDensity = 0.1,      
                          gridSize = 200,      
                          type = 'partition')

colors<-c("grey","#FB8072")
plot6<-DimPlot(af ,reduction= 'umap', group.by = "Type",
               label =TRUE,pt.size = 1)#+ scale_color_manual(values=colors)+NoLegend()
  +geom_path(data=maskTable[cluster=='plaque'],######select a cell
                                                  aes(x=umap_1, y=umap_2,group=group), 
                                                  linetype=6,size=1)
 
  
plot6


af$Type

library(plot1cell)

human_data<-af
# ============1. Circlize plot to visualize cell clustering and meta data======
#####A clean environment may be needed here; loading other packages may cause errors
####Prepare data for ploting

circ_data <- prepare_circlize_data(human_data, scale = 0.8)

length(names(table(Idents(human_data))))

cluster_colors<-colors7[1:8]

(length(names(table(human_data$Type))))
rep_colors<-two_group_col3
#(length(names(table(human_data$oldType))))
#Typecolor<-colors2[1:27]
###plot and save figures
png(filename =  'circlize_plot.png', width = 8, height = 6,units = 'in', res = 500)
plot_circlize(circ_data, do.label = T, pt.size = 0.2,           
              col.use = cluster_colors ,bg.color = 'white',        
              kde2d.n = 500, repel = T, label.cex = 1.2)

add_track(circ_data, group = "Type",colors = rep_colors, track_num = 2) # add legend
#add_track(circ_data, group = "oldType", colors = Typecolor, track_num = 3)
legend("topright",   
       legend = unique(human_data$Type),   
       col = rep_colors,       pch = 15,   
       cex=0.5,       pt.cex=2,    
       box.lwd = 0,  
       bg=NULL,       y.intersp = 1,    
       box.lty=0)
 legend("topleft",   
        legend = unique(human_data$Celltype),   
        col = cluster_colors,       pch = 15,   
        cex=0.5,       pt.cex=2,    
        box.lwd = 0,  
        bg=NULL,       y.intersp = 1,    
        box.lty=0)


dev.off()

