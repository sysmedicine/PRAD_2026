# Supplementary tables ----
data_dir <- "data"
results <- file.path("result", "tables")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
load(file.path(data_dir, "DEGs/DEG_geo_0.01.RData"))
load(file.path(data_dir, "DEGs/DEG_tcga_0.01.RData"))
geo_deg <- rbind(GEOup_dif_exp_01,GEOdown_dif_exp_01)
geo_deg <- geo_deg[,c('padj','up_down')] 
geo_deg$Genename <- rownames(geo_deg)
geo_deg$Cohortname <- 'GEO'
tcga_deg <- rbind(TCGAup_dif_exp_01,TCGAdown_dif_exp_01)
tcga_deg <- tcga_deg[,c('padj','up_down')] 
tcga_deg$Genename <- rownames(tcga_deg)
tcga_deg$Cohortname <- 'TCGA'
combined_deg <- rbind(geo_deg,tcga_deg)
combined_deg <- combined_deg[,c("Genename","Cohortname","up_down", "padj")]
write.csv(combined_deg,
          file = file.path(results, "combined_degs.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")
load(file.path(data_dir, "modules/list_geo_prad_num.RData"))
head(list_geo_prad_num)
geo_modules <- dplyr::bind_rows(list_geo_prad_num, .id = "List_ID")
geo_modules <- geo_modules[,c('gene','GEO_NUM','MM','Connectivity')]
write.csv(geo_modules,
          file = file.path(results, "geo_module_membership.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")
load(file.path(data_dir, "modules/list_tcga_prad_num.RData"))
head(list_tcga_prad_num)
tcga_modules <- dplyr::bind_rows(list_tcga_prad_num, .id = "List_ID")
tcga_modules <- tcga_modules[,c('gene','TCGA_NUM','MM','Connectivity')]
write.csv(tcga_modules,
          file = file.path(results, "tcga_module_membership.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")
message("Results saved to ", results)

# DEG overlap Venn diagram ----
plots <- "plot"
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
VennDiagram::venn.diagram(
  x = list(A = rownames(TCGAup_dif_exp_01), D = rownames(GEOdown_dif_exp_01),
           B = rownames(GEOup_dif_exp_01), C = rownames(TCGAdown_dif_exp_01)),
  filename = file.path(plots, "deg_overlap_venn_color.png"),
  disable.logging = TRUE,
  height = 4500, width = 4500, resolution = 3000, imagetype = "png",
  col = "transparent", fill = c("#FAE48E", "#FF8993", "#A0CDDE", "#9EFA8E"),
  alpha = 0.5, cex = 0.45, cat.cex = 0
)

message("Plots saved to ", normalizePath(plots))

# Archived overlap plots ----
if (FALSE) {
prad_legacy_root <- getwd()
{library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)
library(forcats)
require(matrixStats)
require(data.table)
require(igraph)
require(Hmisc)
require(compiler)
require(reshape2)
require(gridExtra)
require(ggpmisc)  
require(ggrepel)
library(tibble)
library("pheatmap")
library (VennDiagram)
}
setwd(file.path(prad_legacy_root, 'data/legacy_project/part3'))
# Jaccard overlap heatmaps ----
load('rdata/list_geo_prad_num.RData')
load('rdata/list_tcga_prad_num.RData')
load('rdata/11.2_jaccard_dd.RData')
load('rdata/11.1_jaccard_upup.RData')
  disease="AD"
  degree="severe"
  fdr_th=0.05
  int_super_pathway<-"others"#----change
  
  load('rdata/heatmap_up_tcgageo.RData')
  up_matirx <- similarity_matrix_up
  range(similarity_matrix_up)
  my.breaks <- c(seq(0, 0.1, by=0.01))
  my.colors <- c(colorRampPalette(colors = c( "#FEF9E7","#E67E22"))(length(my.breaks)))
  
  label_matrix<-matrix("",dim(dfupup)[1],dim(dfupup)[2])
  index_sig<-which(dfupup < 0.05)
  label_matrix[index_sig]<-"*"
  
  #up
  png("plot/heatmap_up_jaccard.png",units="in", width=10, height=5,res=900)
  pheatmap(up_matirx,
           cluster_cols = T,
           cluster_rows = T,
           clustering_method = "ward.D2",
           color =  my.colors, 
           breaks = my.breaks,
           legend_breaks = seq(0,0.1,by=0.05),
           display_numbers = label_matrix,
           fontsize_number = 20,
           fontsize = 14,
           angle_col = 315,
           treeheight_row =40,
           treeheight_col =40,
           cellwidth = 40,
           cellheight = 30)
  #main = "Overlpa of modules in TCGA GEO183019 (up_regulated)")
  #gaps_col = 2,
  #labels_col = c("CMA","Placebo"),
  #labels_row = rownames(up_matirx))
  dev.off()
#down
load('rdata/heatmap_down_tcgageo.RData')
down_matrix <- similarity_matrix_down
range(down_matrix)
my.breaks <- c(seq(0, 0.1, by=0.01))
my.colors <- c(colorRampPalette(colors = c("#EBF5FB","#2874A6"))(length(my.breaks)))
label_matrix<-matrix("",dim(dfdown)[1],dim(dfdown)[2])
index_sig<-which(dfdown < 0.05)
label_matrix[index_sig]<-"*"
png("plot/heatmap_down_jaccard.png",units="in", width=10, height=5,res=900)
pheatmap(down_matrix,
         cluster_cols = T,
         cluster_rows = T,
         clustering_method = "ward.D2",
         color =  my.colors, 
         breaks = my.breaks,
         legend_breaks = seq(0,0.1,by=0.05),
         display_numbers = label_matrix,
         fontsize_number = 20,
         fontsize = 14,
         angle_col = 315,
         treeheight_row =40,
         treeheight_col =40,
         cellwidth = 40,
         cellheight = 25)
#main = "Overlpa of modules in TCGA GEO183019 (up_regulated)")
#gaps_col = 2,
#labels_col = c("CMA","Placebo"),
#labels_row = rownames(up_matirx))
dev.off()
# DEG overlap Venn diagram ----
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/dif_exp_gene_up_down05.RData'))
res01 <- filter(res0_05 , padj < 0.01)
res01up <- filter(res01 , up_down == 'up')
res01down <- filter(res01 , up_down == 'down')
load(file.path(prad_legacy_root, 'data/legacy_project/part1/output/2 DEG/dif_exp_gene_up_down001.RData'))
geo_up <- rownames(res01up)
geo_down <- rownames(res01down)
tcga_up <- rownames(up_dif_exp_01)
tcga_down <- rownames(down_dif_exp_01)
venn.diagram(x= list(A = tcga_up,D = geo_down,B = geo_up,C = tcga_down), 
             filename = "plot/venn/My3.png", 
             height = 4500, 
             width = 4500,
             resolution =3000, 
             imagetype="png", 
             col ="transparent", 
             fill =c("#FAE48E","#FF8993","#A0CDDE","#9EFA8E"),
             alpha = 0.5, 
             # label.col = c("orange","white","darkorchid4","white","white", "white","white",
             #               "white","darkblue", "white", "white","white","white","darkgreen", 
             #               "white"), 
             cex = 0.45,
             # fontfamily = "serif", 
             # fontface = "bold",
             # cat.col =c("darkblue", "darkgreen", "orange","darkorchid4"), 
             cat.cex = 0,
             # cat.pos = 0, 
             # cat.dist = 0.07,
             # cat.fontfamily = "serif", 
             # rotation.degree = 270
             )
}
