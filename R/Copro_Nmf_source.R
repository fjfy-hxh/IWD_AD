# ReactomePA package must be installed on first use
# BiocManager::install("ReactomePA")

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA) # Core: load the ReactomePA package

run_NMF_Reactome_analysis <- function(signature_matrix, 
                                      top_n = 100, 
                                      color_palette = "Spectral", 
                                      out_dir = "NMF_Reactome_Outputs") {
  
  plot_dir <- file.path(out_dir, "Loading_Plots")
  reactome_dir <- file.path(out_dir, "Reactome_Results")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reactome_dir, recursive = TRUE, showWarnings = FALSE)
  
  n_factors <- nrow(signature_matrix)
  message(paste("Starting Reactome pathway enrichment, total", n_factors, "factors..."))
  
  for (i in 1:n_factors) {
    message(paste("--> Processing Factor", i))
    
    # ---- 1. Data sorting and truncation ----
    vec <- signature_matrix[i, ]
    sorted_vec <- sort(vec, decreasing = TRUE)
    
    dat_plot <- data.frame(
      gene = names(sorted_vec),
      gene_loading = as.numeric(sorted_vec),
      rank = 1:length(sorted_vec),
      stringsAsFactors = FALSE
    )
    dat_plot$gene <- trimws(toupper(dat_plot$gene))
    
    if (nrow(dat_plot) > top_n) {
      dat_plot <- dat_plot[1:top_n, ]
    }
    
    # ---- 2. Draw and save the loading curve ----
    p_loading <- ggplot(dat_plot, aes(x = rank, y = gene_loading, color = gene_loading)) +
      geom_line(color = "grey90", size = 0.8) +
      geom_point(size = 1.2) +
      scale_color_gradientn(colors = rev(brewer.pal(11, color_palette))) +
      geom_text_repel(data = dat_plot[1:5, ], aes(label = gene), size = 2.5, box.padding = 0.2, point.padding = 0.2, max.overlaps = Inf, show.legend = FALSE) + 
      labs(title = paste("Factor", i), x = "Rank", y = "Loading") +
      theme_bw() +
      theme(plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt"), legend.position = "right", legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.5, "cm"), legend.text = element_text(size = 6), legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 7), panel.grid = element_blank()) +
      coord_cartesian(clip = "off")
    
    ggsave(filename = file.path(plot_dir, paste0("Factor_", i, "_Gene_rank.pdf")), plot = p_loading, width = 4.5, height = 3.5)
    
    # ---- 3. Reactome enrichment analysis ----
    # Genes must first be converted to ENTREZID before being passed to Reactome
    gene_ids <- tryCatch({
      bitr(dat_plot$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
    }, error = function(e) { return(NULL) })
    
    if (is.null(gene_ids) || nrow(gene_ids) == 0) {
      message(paste("   [Warning] Factor", i, "gene conversion to ENTREZID failed, skipping."))
      next
    }
    
    # First run: strict filtering mode
    eract <- tryCatch({
      enrichPathway(gene          = gene_ids$ENTREZID,
                    organism      = "human", # For mouse, change to "mouse"
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2)
    }, error = function(e) { return(NULL) })
    
    # Second run: fallback (relax the filtering criteria)
    if (is.null(eract) || nrow(eract) == 0) {
      eract <- tryCatch({
        enrichPathway(gene          = gene_ids$ENTREZID,
                      organism      = "human", 
                      pAdjustMethod = "none",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 1,
                      minGSSize     = 1)
      }, error = function(e) { return(NULL) })
      
      if (!is.null(eract) && nrow(eract) > 0) {
        message(paste0("   [Notice] Factor ", i, " strict filtering returned no results, switched to relaxed Reactome mode."))
      }
    }
    
    # ---- 4. Save and plot ----
    if (!is.null(eract) && nrow(eract) > 0) {
      # Convert ENTREZID back to SYMBOL so the genes are readable in the CSV
      eract <- setReadable(eract, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      # Save the CSV table
      write.csv(as.data.frame(eract), 
                file = file.path(reactome_dir, paste0("Factor_", i, "_Reactome_enrichment.csv")), 
                row.names = FALSE)
      
      # Draw and save the dotplot bubble chart
      p_react <- dotplot(eract, showCategory = 10, title = paste("Factor", i, "Reactome Enrichment")) +
        theme_bw() + 
        scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          panel.grid.major = element_line(color = "grey95"),
          panel.grid.minor = element_blank()
        )
      
      ggsave(filename = file.path(reactome_dir, paste0("Factor_", i, "_Reactome_dotplot.pdf")), 
             plot = p_react, width = 6.0, height = 4.5) # Slightly wider, since Reactome names are also quite long
    } else {
      message(paste("   [Warning] Factor", i, "has no enriched pathways at all in the Reactome database."))
    }
  }
  message("Reactome pipeline completed!")
}
################
ks_plot_nmf_factor <- function(nmf_res, 
                                   k_blocks = 10, 
                                   low_color = "#E0F3F8",  # 👈 Low-similarity color (default light blue-white)
                                   high_color = "#4A0080", # 👈 High-similarity color (default deep pure purple)
                                   use_zero_center = FALSE) { # Whether to center symmetrically at 0 (red-white-blue transition)
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. Extract the matrix and clean up the dimensions
  sig_mat <- nmf_res$signatures
  if (nrow(sig_mat) > ncol(sig_mat)) { sig_mat <- t(sig_mat) }
  rownames(sig_mat) <- paste0("factor", 1:nrow(sig_mat))
  
  # 2. Compute the correlation matrix
  cor_matrix <- cor(t(sig_mat), method = "pearson")
  diag(cor_matrix) <- 1
  hc <- hclust(dist(1 - cor_matrix), method = "ward.D2")
  
  min_val <- min(cor_matrix)
  
  # 3. Core change: build a fully continuous color mapping
  if (use_zero_center) {
    # If zero-centering is enabled, use a smooth continuous blue -> white -> red transition
    col_fun = colorRamp2(c(-1, 0, 1), c("#2166AC", "#FFFFFF", "#B2182B"))
  } else {
    # Linear interpolation between the two endpoints, fully continuous colors with no blocking
    col_fun = colorRamp2(c(min_val, 1), c(low_color, high_color))
  }
  
  # 4. Draw the heatmap
  p <- Heatmap(
    cor_matrix, 
    name = "Similarity", 
    col = col_fun,
    cluster_rows = hc,
    cluster_columns = hc,
    row_split = k_blocks,
    column_split = k_blocks,
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    rect_gp = gpar(col = "#E5E5E5", lwd = 0.5), 
    border = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontface = "bold", size = 10),
    column_names_gp = gpar(fontface = "bold", size = 10),
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      legend_height = unit(4, "cm")
    )
  )
  
  return(draw(p))
}
##################
###########
#NMF top marker heatmap
ks_plot_factor_markers <- function(nmf_res, 
                                   top_n = 5, 
                                   color_palette = c("#E0F2F1", "#26A69A", "#0A2240", "#050811")) { # 👈 Default is still the original refined palette
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. Extract the raw matrix
  sig_mat <- nmf_res$signatures
  if (is.null(rownames(sig_mat))) {
    rownames(sig_mat) <- paste0("F", 1:nrow(sig_mat))
  }
  
  # 2. Extract the top genes for each Factor
  top_genes_list <- lapply(1:nrow(sig_mat), function(i) {
    names(sort(sig_mat[i, ], decreasing = TRUE)[1:top_n])
  })
  selected_genes <- unique(unlist(top_genes_list))
  sub_mat <- sig_mat[, selected_genes, drop = FALSE]
  
  # 3. Row normalization (Min-Max Scaling 0-1)
  sub_mat_scaled <- t(apply(sub_mat, 1, function(x) {
    if(max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }))
  
  # 4. [Core upgrade]: automatically split 0 to 1 into equal breaks based on the number of colors passed in, producing a fully continuous gradient
  breaks <- seq(0, 1, length.out = length(color_palette))
  col_fun = colorRamp2(breaks, color_palette)
  
  # 5. Start drawing the heatmap
  p <- Heatmap(
    sub_mat_scaled, 
    name = "Gene contribution to factor (scaled)", 
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = "Factor",
    row_title_gp = gpar(fontface = "bold", size = 12),
    row_names_side = "left",
    row_names_gp = gpar(fontface = "bold", size = 10),
    column_names_rot = 90,
    column_names_gp = gpar(fontface = "italic", size = 9),
    rect_gp = gpar(col = "white", lwd = 0.2), # Thin light lines separating each cell
    border = TRUE,
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "leftcenter",
      legend_width = unit(5, "cm")
    )
  )
  
  return(draw(p, heatmap_legend_side = "top"))
}
##############
##NMF factor-distance function
ks_plot_factor_trajectories <- function(seurat_obj,
                                          factor_names = paste0("factor", 1:10),
                                          distance_var = "vertical_distance",
                                          custom_colors = "turbo") {
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # 1. Safely extract the specified factors and distance metadata from the Seurat object
  vars_to_fetch <- c(factor_names, distance_var)
  if (!all(vars_to_fetch %in% colnames(seurat_obj@meta.data)) && 
      !all(vars_to_fetch %in% rownames(seurat_obj))) {
    # Compatible with packages that store factors in metadata and others that store them in an Assay
    st_data <- FetchData(seurat_obj, vars = vars_to_fetch)
  } else {
    st_data <- FetchData(seurat_obj, vars = vars_to_fetch)
  }
  
  # 2. Convert to long format and normalize each Factor's weight to 0-1 independently within its group
  st_data_long <- st_data %>%
    pivot_longer(
      cols = all_of(factor_names), 
      names_to = "Factor", 
      values_to = "Weight"
    ) %>%
    group_by(Factor) %>%
    mutate(Weight_Scaled = (Weight - min(Weight, na.rm = TRUE)) / 
             (max(Weight, na.rm = TRUE) - min(Weight, na.rm = TRUE))) %>%
    ungroup()
  
  # Lock the factor ordering of Factor
  st_data_long$Factor <- factor(st_data_long$Factor, levels = factor_names)
  
  # Rename the distance column to make ggplot mapping easier
  colnames(st_data_long)[colnames(st_data_long) == distance_var] <- "Distance"
  
  # 3. Start building the base plot in a classic academic style
  p <- ggplot(st_data_long, aes(x = Distance, y = Weight_Scaled, color = Factor, fill = Factor)) +
    # Fit a smooth curve using a GAM model
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 1.2, alpha = 0.05) +
    # Use the high-transparency classic academic theme standard in high-impact papers
    theme_classic() +
    labs(
      x = "Vertical Distance",
      y = "Relative Factor Value (Scaled 0-1)",
      title = "Combined Trajectories of NMF Factors along Vertical Distance",
      color = "NMF Factors"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black", size = 10),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    ) +
    # Perfectly trim the legend: forcibly drop the redundant fill legend, keeping only the colored line legend
    guides(fill = "none") 
  
  # 4. Core: custom color control logic
  # Option A: a built-in viridis palette name string
  if (length(custom_colors) == 1 && 
      custom_colors %in% c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")) {
    p <- p + scale_color_viridis_d(option = custom_colors)
  } else {
    # Option B: a custom color vector was passed in (e.g. 10 custom hex color codes)
    if (length(custom_colors) < length(factor_names)) {
      warning("Fewer custom colors than factors; some lines may lose their color, the default palette has been applied automatically.")
      p <- p + scale_color_viridis_d(option = "turbo")
    } else {
      p <- p + scale_color_manual(values = custom_colors)
    }
  }
  
  return(p)
}
###############
#Compute approximate factor values
ks_plot_density_summary <- function(obj, 
                                    feature = "factor7", 
                                    group_by = "Type", 
                                    custom_colors = NULL,
                                    y_limits = NULL,     # New: truncate the Y-axis range, e.g. c(0, 5)
                                    x_limits = NULL) {   # New: truncate the X-axis range, e.g. c(0, 10)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  
  # 1. Auto-compatible: supports passing a Seurat object directly or a meta matrix
  if (inherits(obj, "Seurat")) {
    meta_df <- obj@meta.data
  } else if (is.data.frame(obj)) {
    meta_df <- obj
  } else {
    stop("Invalid input object: must be a Seurat object or a metadata data.frame")
  }
  
  # Check whether the column names exist
  if (!feature %in% colnames(meta_df)) stop(paste("Feature column not found in data:", feature))
  if (!group_by %in% colnames(meta_df)) stop(paste("Grouping column not found in data:", group_by))
  
  # 2. Summary statistics (fully aligned to the specified feature)
  summary_table <- meta_df %>%
    dplyr::group_by(across(all_of(group_by))) %>%
    dplyr::summarise(
      mean_score   = mean(.data[[feature]], na.rm = TRUE),
      median_score = median(.data[[feature]], na.rm = TRUE),
      sd_score     = sd(.data[[feature]], na.rm = TRUE),
      cell_count   = dplyr::n(),
      .groups = "drop"
    )
  
  # Print the tidy summary statistics to the console
  cat("\n📊 ---", feature, "summary statistics across", group_by, "groups ---\n")
  print(as.data.frame(summary_table))
  cat("----------------------------------------------------\n\n")
  
  # 3. Draw an academic-style density plot
  p <- ggplot(meta_df, aes(x = .data[[feature]], fill = .data[[group_by]])) +
    geom_density(alpha = 0.4, color = "#333333", linewidth = 0.5) + # Note: the size argument is now recommended to be linewidth in newer ggplot versions
    labs(
      x = paste(feature, "Score"),
      y = "Density",
      title = paste("Density Distribution of", feature),
      fill = group_by
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black", size = 10),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  
  # 4. Custom color control
  if (!is.null(custom_colors)) {
    p <- p + scale_fill_manual(values = custom_colors)
  } else {
    p <- p + scale_fill_brewer(palette = "Set2")
  }
  
  # 5. Axis truncation (Zoom in - core change)
  # coord_cartesian safely zooms the view without dropping data points
  if (!is.null(y_limits) || !is.null(x_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, xlim = x_limits)
  }
  
  return(p)
}
################
##Draw the factor spatial plot
ks_plot_spatial_feature <- function(seurat_obj,
                                    feature = "factor8",
                                    group_id = "A9",
                                    group_col = "group",
                                    x_col = "x_centroid",
                                    y_col = "y_centroid",
                                    reduction_name = "spatial",
                                    pt_size = 0.1,
                                    limits = c(400, 1400),
                                    xlim = NULL,
                                    ylim = NULL,
                                    show_axes = FALSE, # 👈 [New parameter]: whether to display axes and grid lines for positioning
                                    colors = c('#6090C1','#ACD2E5','#FEF9B7','#FEE395','#F2724D','#D7312D')) {
  library(Seurat)
  library(ggplot2)
  library(scales)
  
  # 1. Safely filter the sample/group
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Grouping column not found in metadata:", group_col))
  }
  
  plot_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_col]] %in% group_id, ]))
  
  # 2. Automatically inject centroid coordinates into the specified reduction
  plot_obj <- add_spatial_centroid_reduction(
    seurat_obj     = plot_obj, 
    x_col          = x_col, 
    y_col          = y_col, 
    reduction_name = reduction_name
  )
  
  # 3. Base FeaturePlot rendering
  p <- FeaturePlot(
    object = plot_obj, 
    features = feature, 
    reduction = reduction_name,
    pt.size = pt_size
  )
  
  # 4. Overlay the advanced gradient and truncation mechanism
  p <- p + scale_color_gradientn(
    colours = colorRampPalette(colors)(100), 
    limits  = limits,           
    oob     = scales::squish, 
    name    = ""
  )
  
  # 5. [Core change]: dynamically control the layer styling based on show_axes
  if (show_axes) {
    # 🕵♂️ Auxiliary positioning mode: show axis ticks, labels and a light grid
    p <- p + theme_bw() + 
      theme(
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        
        # Show axes and grid
        axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black", size = 9),
        panel.grid.major = element_line(color = "gray90", size = 0.3), # Auxiliary reference lines
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Spatial X Centroid", y = "Spatial Y Centroid")
  } else {
    # 📄 Final publication mode: completely strip off the axes and cluttered grid
    p <- p + theme_void() + 
      theme(
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2) 
      )
  }
  
  # 6. Dynamically control local window cropping and zooming
  p <- p + coord_fixed(ratio = 1, xlim = xlim, ylim = ylim)
  
  return(p)
}
##############
#Run NMF multiple times to compute AIC
ks_select_nnmf_k <- function(data, 
                             location, 
                             batch_id, 
                             k_candidates = 2:10) {
  library(ggplot2)
  library(dplyr)
  
  # 1. Basic dimension extraction (sum the row and column totals, immune to the risk of genes and spots being transposed)
  n_dim1 <- nrow(data)
  n_dim2 <- ncol(data)
  
  # Initialize the storage table
  aic_tb <- data.frame(
    K = integer(),
    GKL_loss = numeric(),
    param_num = integer(),
    AIC = numeric()
  )
  
  # 2. Loop over each K, fit NNMF and compute AIC
  for (k in k_candidates) {
    message("Fitting K = ", k)
    
    # Fit the nnmf model
    fit <- nnmf(
      data = data,
      noSignatures = k,
      location = location,
      batch = batch_id
    )
    
    # Extract the GKL reconstruction loss D(V||WH)
    gkl <- fit$error
    
    # Compute the total number of parameters
    nprm <- k * (n_dim1 + n_dim2)
    
    # Compute the AIC value
    aic <- 2 * gkl + 2 * nprm
    
    # Store into the table
    aic_tb <- rbind(aic_tb, data.frame(
      K = k,
      GKL_loss = gkl,
      param_num = nprm,
      AIC = aic
    ))
    
    # Release memory promptly to prevent out-of-memory during high-rank fits
    gc() 
  }
  
  # 3. Select the optimal K (minimum AIC)
  opt_K <- aic_tb$K[which.min(aic_tb$AIC)]
  min_aic <- min(aic_tb$AIC)
  
  cat("\n🏆 Done! Optimal noSignatures (K) with the lowest AIC =", opt_K, "\n\n")
  
  # 4. Use theme_classic to draw a publication-grade AIC curve
  p <- ggplot(aic_tb, aes(x = K, y = AIC)) +
    # Draw the base line and data points
    geom_line(color = "#333333", size = 1) +
    geom_point(color = "#333333", size = 3) +
    
    # Highlight the optimal K point in red (using the refined inverted triangle shape = 17)
    geom_point(data = filter(aic_tb, K == opt_K), aes(x = K, y = AIC), 
               color = "#D7312D", shape = 17, size = 5) +
    
    # Dynamically add a bold academic text label for the optimal K
    annotate("text", x = opt_K, y = min_aic, 
             label = paste0("Optimum K = ", opt_K), 
             vjust = -1.8, color = "#D7312D", fontface = "bold", size = 4.5) +
    
    # Force the X-axis ticks to strictly match your candidate K range
    scale_x_continuous(breaks = k_candidates) +
    
    # Classic academic theme styling
    theme_classic() +
    labs(
      x = "noSignatures (K)",
      y = "AIC Value",
      title = "AIC Curve for Selecting Optimal NMF Signatures"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black", size = 10)
    )
  
  # Display the rendered figure directly in the plot window
  print(p)
  
  # Return the result list
  return(list(
    aic_table = aic_tb,
    optimal_K = opt_K,
    plot = p
  ))
}
