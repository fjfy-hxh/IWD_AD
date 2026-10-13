
library(NNMF)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(magrittr)
library(writexl)
library(gprofiler2)
library(pheatmap)
library(CoPro)
source(file = "./R/Copro_Nmf_source.R")

###
af.nmf_copro<-af.ST
# 
data_raw <- t(as.matrix( af.nmf_copro@assays$Spatial$data))
# 
location_raw <- Embeddings( af.nmf_copro, reduction = "spatial")

#
keep_cells <- rowSums(data_raw) > 0
data <- data_raw[keep_cells, ]
location <- location_raw[keep_cells, ]
keep_genes <- colSums(data) > 0
data <- data[, keep_genes]
# 
batch_id = groupondist(location, size = 10000)
min_max_norm <- function(coord_vec) {
  # Formula: (x - min(x)) / (max(x) - min(x))
  return((coord_vec - min(coord_vec)) / (max(coord_vec) - min(coord_vec)))
}
location[, 1] <- min_max_norm(location[, 1])
location[, 2] <- min_max_norm(location[, 2])
str(data)
# 
res = nnmf(data = data, noSignatures = 10, location = location, batch = batch_id)
#qsave(v, file = "res.qs")
#res<-qread("res.qs")
res$weights # The weights for each signature
res$signatures # The signatures

ks_plot_nmf_factor(res,  k_blocks = 10, 
                   low_color = "#EAECEE", high_color = "#4A0080")
nrow(res$signatures)
# 
ks_plot_factor_markers(res, top_n = 5,
                       color_palette = c("#E0F2F1", "#26A69A", "#0A2240"))

top_gene = as.data.frame(t(topfeatures(signatures = res$signatures,
                       feature_names = names(keep_genes), ntop = 100)))[-1, ]


rownames(res$weights)=rownames(data)
colnames(res$weights) <- paste0("factor", 1:ncol(res$weights))
rownames(res$signatures) <- paste0("factor", 1:ncol(res$weights))
colnames(res$signatures) <- colnames(data)
weights_mat <- res$weights[colnames(af.ST), ]

top_name = as.data.frame(t(res$signatures))
class(top_name)

# 
af.ST <- AddMetaData(
  object = af.ST,
  metadata = weights_mat
)
class(dat$metaData)
#################
#
#
nmf_selection <- ks_select_nnmf_k(
  data         = data, 
  location     = location, 
  batch_id     = batch_id, 
  k_candidates = 3:12
)

# 
head(nmf_selection$aic_table)
# 
best_k <- nmf_selection$optimal_K

res1 <- nnmf(
  data = data,
  noSignatures = best_k,
  location = location,
  batch = batch_id
)


colnames(af.ST@meta.data)

nmf_factor=c('factor1','factor2','factor3','factor4','factor5','factor6',
             'factor7','factor8','factor9','factor10')
library(MySeuratWrappers)
######
af.nmfplot<- subset(af.ST, niche %in% c('niche_4'))
Idents(af.nmfplot)<-af.nmfplot$Type
x <- Idents(af.nmfplot)
levels(x) <- c("A", "M", "Y")
Idents(af.nmfplot) <- x
table(af.nmfplot$Type)
VlnPlot(af.nmfplot, features = nmf_factor,
        stacked=T,pt.size=0,
        cols = colors_Type,
        direction = "horizontal", #
        x.lab = '', y.lab = '')+#
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())#




# 
run_NMF_Reactome_analysis(
  signature_matrix = res$signatures, 
  top_n = 200, # 
  color_palette = "Spectral", 
  out_dir = "NMF_Reactome_Plots"
)
dim(res$signatures)



af.niche1<- subset(af.ST, niche %in% c('niche_1'))

###
ks_plot_factor_trajectories(af.ST, 
                            factor_names = paste0("factor", 1:10),
                            distance_var = "vertical_distance",
                            custom_colors = colors8)
table(af.epi1$radius_range)
##
ks_plot_density_summary(af.ST, feature = "Fib_1", group_by = "Type",
                        custom_colors=colors_Type#,x_limits = c(0,1000)
                        )
###
ks_plot_spatial_feature(
  seurat_obj = af.ST,
  feature    = "min_distance",
  group_id   = "A11",
  # xlim       = c(500, 2500), #
  # ylim       = c(500, 2500), # 
  pt_size    = 0.4 ,          #
  limits     = c(0, 300),
  show_axes  = F   # 
)


ggsave(filename = "Y10 factor9.png",  width  = 5,  height   = 5,  units = "in",dpi= 300)


