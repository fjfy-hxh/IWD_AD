ks_CC_groupCells_bubble <- function(object_list,# cellchat objects stored as a list, e.g. group A, B and C each run cellchat and get a result; object_list = list(A, B, C)
                                    group_name,# group names as a vector, must correspond to object_list
                                    celltypes, # celltypes whose interactions are to be displayed
                                    LR_int, # ligand-receptor pairs to display
                                    cols=NULL,# density color
                                    cell_cols,# cell color
                                    group_cols# group color
){ 
  
  #p
  require(CellChat)
  require(ggplot2)
  require(ggh4x)
  
  
  
  
  if(is.list(object_list)==TRUE){
    
    
    print("go on")
    
  }else{
    
    print("Please check it")
  }
  
  
  net <- list()
  
  
  for (i in 1:length(object_list)) {
    
    net1 <- subsetCommunication(object_list[[i]])
    net[[i]] <- net1
    names(net)[i] <- group_name[i]
  }
  
  
  
  for (i in 1:length(net)) {
    
    net[[i]] <-  net[[i]] %>% filter(interaction_name %in%  LR_int) %>%
      filter(source %in% celltypes) %>% filter(target %in% celltypes)
    
    net[[i]]$ST <- paste0(net[[i]]$source,"_",net[[i]]$target,"~",net[[i]]$interaction_name)
    
    names(net)[i] <- group_name[i]
    
  }
  
  
  combinations <- expand.grid(a = celltypes, b = celltypes)
  combinations <- paste(combinations$a, combinations$b, sep = "_")
  com <- expand.grid(a = combinations, b = LR_int)
  com <- paste(com$a, com$b, sep = "~")
  
  
  spe_pa <- list()
  
  for (i in 1:length(net)) {
    
    spe_pa[[i]] <- setdiff(com, unique(net[[i]]$ST))
    
  }
  
  
  #0
  is_empty <- lapply(spe_pa, function(x) length(x) == 0)
  
  #fix group data
  df.net <- list()
  
  for (i in 1:length(net)) {
    
    if(is_empty[[i]]==T){
      
      df.net1 <- as.data.frame(matrix(NA, nrow = length(spe_pa[[i]]), ncol = ncol(net[[i]]), dimnames = list(NULL, colnames(net[[i]])))) 
    }else{
      
      df.net1 <- as.data.frame(matrix(NA, nrow = length(spe_pa[[i]]), ncol = ncol(net[[i]]), dimnames = list(NULL, colnames(net[[i]])))) 
      df.net1$ST <- spe_pa[[i]]
      str_a <- stringr::str_split(df.net1$ST, "~", simplify = T)
      df.net1$interaction_name <- str_a[,2]
      df.net1$source <- stringr::str_split(str_a[,1], "_", simplify = T)[,1]
      df.net1$target <- stringr::str_split(str_a[,1], "_", simplify = T)[,2]
      df.net1$prob <- 0
      df.net1$pval <- 1
    }
    
    
    
    df.net[[i]] <- df.net1
  }
  
  
  
  
  for (i in 1:length(net)) {
    
    if(is_empty[[i]]==T){
      
      net[[i]] <- rbind(net[[i]], df.net[[i]])
      net[[i]] <- net[[i]][-nrow(net[[i]]),]
      
    }else{
      
      net[[i]] <- rbind(net[[i]], df.net[[i]])
      
    }
    
  }
  
  
  
  for (i in 1:length(net)) {
    
    net[[i]]$group <- group_name[i]
    
  }
  
  
  #data final
  data_df <- do.call(rbind, net)
  
  if(is.null(cols)){
    cols = RColorBrewer::brewer.pal(n = 10, name = "Spectral")
    
  }else{
    
    cols=cols
  }
  
  
  
  data_df$target <- factor(data_df$target, levels = celltypes)
  data_df$source <- factor(data_df$source, levels = celltypes)
  data_df$interaction_name <- factor(data_df$interaction_name, levels = LR_int)
  
  
  data_df$pval <- -log10(data_df$pval+0.001)
  
  
  p =  ggplot(data_df, aes(target,interaction_name))+
    geom_point(aes(size= pval, color=prob))+
    geom_point(aes(size= pval), shape=21,stroke=1)+
    facet_grid2(group~source,switch = "x",
                strip  = strip_nested(background_x = elem_list_rect(fill = alpha(cell_cols,0.5),color = NA),
                                      background_y = elem_list_rect(fill = alpha(group_cols,0.5),color = NA)))+
    theme_bw()+
    scale_x_discrete(position = "top")+
    theme(axis.text = element_text(color='black', size=9),
          axis.text.x.top = element_text(angle = 90, hjust = 0),
          axis.title.x.top = element_text(vjust=5),
          plot.margin = margin(12,12,12,12),
          strip.text = element_text(size = 10,color='black',face = "bold"),
          strip.background = element_rect(fill = NA),
          legend.frame = element_rect(colour = "black"),
          legend.key.width = unit(0.3, "cm"),# modify the legend width
          legend.key.height = unit(0.5, "cm"),# modify the legend height
          legend.title = element_text(color = 'black', face = "bold", size=9))+
    scale_color_gradientn(colours = cols)+
    guides(size=guide_legend(title="-log10(pval)"), # set the legend size title
           colour=guide_colorbar(title="commun.\nprob"))+# set the legend color bar title
    labs(x="Source cells",
         y="")
  
  
  
  return(p)
  
}

