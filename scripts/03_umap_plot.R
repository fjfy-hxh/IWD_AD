#http://www.360doc.com/content/24/0430/20/82995563_1122020662.shtml


#############################
source("./R/umap_with_unchull.R")



colors_domain<-c('#00afbb','#fe9a30')


#plotData$UMAP_1
p <- plot_umap_with_unchull(
  seurat_obj = af,
  cluster_col = "Celltype",
  colors = colors3,
  point_size = 0.03#0.2-0.3
)

print(p)


pdf(file = "UMAP.pdf",width =6,height = 4)
p
dev.off()

###
p1 <- plotUMAP(
  seurat_obj = af,
  reduction = "umap",
  cluster_col = "Type",
  color_palette = cwd_color,
  axis_line_size = 0.4,      # 
  legend_point_size = 3,     #
  legend_text_size = 5,      # 
  legend_ncol = 1 ,          # 
  legend_alpha = 0.7 
)
p1
ggsave(filename ="scRNA_Type UMAP.png",height = 3,width = 4,plot = p1,dpi = 300,bg="white")






