#########Ucell########


#
source('./R/ucell_pathway_scoring.R')
library(UCell)
library(stringr)
library(viridis)
library(clusterProfiler)
library(msigdbr)
#




#####################

senescence_gs <- read.gmt("REACTOME_CELLULAR_SENESCENCE.v2026.1.Hs.gmt")
senescence_genes <- senescence_gs$gene

cell_color <- c(
  "A" = "#F3AE63",
  "M" = "#73558B",
  "Y"   = "#93CBAE"
)
# Specify comparison groups
my_comps <- list(
  c("A", "M"),
  c("A", "Y"),
  c('M','Y')
)

# Run
af.ST <- RunSinglePathwayUCell(
  seu_object     = afd,
  pathway_genes  = senescence_genes,
  pathway_name   = "REACTOME_CELLULAR_SENESCENCE",
  cell_color     = cell_color,
  group_col      = "Type",
  cell_levels    = c("A", "M", "Y"),
  my_comparisons = my_comps
)


table(af.ST$dist)
