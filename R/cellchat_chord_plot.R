plot_cellchat_LR <- function(cellchat_net,#net.df <- subsetCommunication(cellchat_obj)
                             celltype_inter, # celltypes to display
                             celltype_color, # color setting for the displayed celltypes
                             ligand_col="#FB8072",# ligand color
                             receptor_col="#1988B0",# receptor color
                             text_size,# text size
                             legend_in_plot=T,# whether to show the legend in the plot; if F, save it to a separate pdf
                             Group) { # grouping, mainly used to distinguish groups
  
  #load packages
  
  library(tidyverse)
  library(circlize)
  library(ggsci)
  library(igraph)
  library(gtools)
  library(ComplexHeatmap)
  

  
  net.df <- cellchat_net # get the pairs of cell types or ligand receptors in the interactions
  net.df$cell_type_pair <- paste0(net.df$source, '_', net.df$target) # build the pairs of interaction cell types
  
  node_v1 <- vector() # the first node
  node_v2 <- vector() # the second node
  weight_v <- vector() # the edge weight
  
  for (i in 1:nrow(net.df)) {
    cci <- strsplit(net.df$interaction_name_2[i], split = ' ')[[1]] # split the interaction name
    ct.pair <- strsplit(net.df$cell_type_pair[i], split = '_')[[1]] # get the cell type pair
    prefix1 <- ct.pair[1]
    prefix2 <- ct.pair[2]
    # prefix1 <- substr(ct1, start = 1, stop = 3)
    # prefix2 <- substr(ct2, start = 1, stop = 3)
    node1 <- paste0(prefix1, '_', cci[1])
    id2 <- cci[length(cci)]
    weight <- net.df$prob[i]
    
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])))
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1))
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2)
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) {
    return()
  }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  # write.csv(x = as.data.frame(as.matrix(adj)), file = paste0(Group, '_mat.csv'))
  
  graph_adj <- as.data.frame(as.matrix(adj))
  
  
  Genes <- names(V(g))
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    return(x[1])
  })
  category <- data.frame(Genes = Genes, ID = ID)
  # write.csv(x = category, file = paste0(Group, '_cat.csv'))
  
  graph_module = category
  
  # LR grouping
  graph_module$LR <- ifelse(graph_module$Genes %in% g_df$node_v1,
                            "ligand", "receptor")
  
  
  #plot
  if(is.null(celltype_inter)==F){
    
    graph_module <- graph_module[graph_module$ID %in% celltype_inter,]
    str <- paste0("^(", paste(celltype_inter, collapse = "|"), ")")
    str_cells <- grepl(str, rownames(graph_adj))
    str_cells1 <- rownames(graph_adj)[str_cells]
    str_cells2 <- colnames(graph_adj)[str_cells]
    
    graph_adj <- graph_adj[str_cells1,]
    graph_adj <- graph_adj[,str_cells2]
    
  }
  
  
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  
  #set colors
  LR_color <- data.frame(LR=c("ligand","receptor"),
                         color1=c(ligand_col,receptor_col))
  graph_module <- left_join(graph_module,LR_color, by="LR")
  
  
  
  celltype_color <- as.data.frame(celltype_color)
  celltype_color$color2 <- rownames(celltype_color)
  colnames(celltype_color) <- c("ID",'color2')
  
  graph_module <- left_join(graph_module,celltype_color, by="ID")
  
  
  
  
  raw_edges <-  as.data.frame(cbind(get.edgelist(g), E(g)$weight)) %>%
    mutate(
      V1 = gsub('\\.', '-', V1),
      V2 = gsub('\\.', '-', V2),
      V3 = as.numeric(V3),
      V4 = 1
    )
  edges <- raw_edges %>%
    arrange(V3)
  
  
  #nodes
  nodes <-  unique(c(edges$V1, edges$V2))
  
  sectors <- sort(unique(c(raw_edges$V1, raw_edges$V2)))
  
  
  #https://vuetifyjs.com/en/styles/colors/#material-colors
  col_fun = colorRamp2(range(edges$V3), c("#FFFDE7", "#013220"))
  
  
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(-0.15,0.2))
  circos.initialize(sectors, xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA)
  
  # we go back to the first track and customize sector labels
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      node_text_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color2) %>%
        as.character()
      
      node_LR_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color1) %>%
        as.character()
      
      
      circos.rect(
        xlim[1],
        0,
        xlim[2],
        1,
        col = node_LR_color,
        border = NA
      )
      
      circos.text(
        mean(xlim),
        2,
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col= node_text_color,
        cex=text_size
      )
    },
    bg.border = NA
  ) 
  
  for (i in seq_len(nrow(edges))) {
    link <- edges[i,]
    circos.link(link[[1]],
                c(0, 1),
                link[[2]],
                c(0, 1),
                col = col_fun(link[[3]]),
                border = NA)
  }
  
  
  title(paste0("  ",Group))
  
  #plot legend
  
  
  if(legend_in_plot==T){
    lgd <- Legend(title ="Score", col_fun = col_fun,direction="horizontal",border="black")
    grid.draw(lgd)
    
    
    legend(1.1,0.5,
           pch=15,
           legend=c("Ligand","Receptor"),
           bty="n",
           col =c(ligand_col,receptor_col),
           cex=1,pt.cex=3,
           border="black") 
    
    legend(1.1,0,
           pch=20,
           legend=celltype_color$ID,
           bty="n",
           col = celltype_color$color2,
           cex=1,pt.cex=3,
           border="black") 
    
    
    
    
  }else{
    
    
    lgd <- Legend(title ="Score", col_fun = col_fun,direction="horizontal",border="black")
    pdf(
      paste(paste0(Group, "_legend.pdf"), sep = ""),
      width = 10,
      height = 10
    )
    grid.draw(lgd)
    dev.off()
    
    
    lgd_LR = Legend(labels = c("Ligand","Receptor"), 
                    legend_gp = gpar(fill = c(ligand_col,receptor_col)))
    
    pdf(
      paste(paste0(Group, "_legend_LR.pdf"), sep = ""),
      width = 10,
      height = 10
    )
    grid.draw(lgd_LR)
    dev.off()
    
    
    lgd_cells = Legend(labels =celltype_color$ID, 
                       legend_gp = gpar(fill = celltype_color$color2))
    
    pdf(
      paste(paste0(Group, "_legend_cells.pdf"), sep = ""),
      width = 10,
      height = 10
    )
    grid.draw(lgd_cells)
    dev.off()
    
    
  }
  
  
}
