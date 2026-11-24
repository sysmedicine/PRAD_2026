# WGCNA module overlap and BH correction ----
data_dir <- "data"
results <- file.path("result", "overlap")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
workbook <- file.path(data_dir, "module_overlap_BH.xlsx")
pairs_per_group <- c(Upregulated = 18L, Downregulated = 24L)
pair_columns <- c("tcga_module", "geo_module", "tcga_size", "geo_size",
            "intersection", "archived_jaccard", "archived_p", "archived_fdr", "archived_star")
read_pairs <- function(group) {
  range <- paste0("A8:I", 7L + pairs_per_group[[group]])
  x <- as.data.frame(readxl::read_excel(workbook, sheet = group, range = range, col_names = pair_columns))
  numeric_columns <- pair_columns[3:8]
  x[numeric_columns] <- lapply(x[numeric_columns], as.numeric)
  x$group <- group
  stopifnot(nrow(x) == pairs_per_group[[group]], !anyNA(x[pair_columns[3:5]]))
  x
}
module_pairs <- do.call(rbind, lapply(names(pairs_per_group), read_pairs))
background_gene_count <- 12788L
module_pairs$jaccard <- with(module_pairs, intersection / (tcga_size + geo_size - intersection))
module_pairs$p <- with(module_pairs, phyper(intersection - 1L, tcga_size, background_gene_count - tcga_size,
                             geo_size, lower.tail = FALSE))
module_pairs$fdr <- p.adjust(module_pairs$p, method = "BH")
module_pairs$significant <- with(module_pairs, p < 0.05 & fdr < 0.05)
stopifnot(identical(module_pairs$significant, !is.na(module_pairs$archived_star) & module_pairs$archived_star == "*"))
write.csv(module_pairs, file.path(results, "wgcna_module_overlap.csv"), row.names = FALSE)
print(with(module_pairs, table(group, significant)))
message("Results saved to ", results)
# Jaccard overlap heatmaps ----
plots <- "plot"
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
heatmap_colors <- list(
  Upregulated = c("#FEF9E7", "#E67E22"),
  Downregulated = c("#EBF5FB", "#2874A6")
)
for (direction in names(heatmap_colors)) {
  panel <- module_pairs[module_pairs$group == direction, ]
  tcga_modules <- unique(panel$tcga_module)
  geo_modules <- unique(panel$geo_module)
  heatmap_values <- matrix(NA_real_, length(tcga_modules), length(geo_modules),
                           dimnames = list(tcga_modules, geo_modules))
  cells <- cbind(match(panel$tcga_module, tcga_modules), match(panel$geo_module, geo_modules))
  stopifnot(!anyDuplicated(as.data.frame(cells)))
  heatmap_values[cells] <- panel$jaccard
  stopifnot(all(is.finite(heatmap_values)))
  stars <- matrix("", nrow(heatmap_values), ncol(heatmap_values), dimnames = dimnames(heatmap_values))
  stars[cells] <- ifelse(panel$significant, "*", "")
  pheatmap::pheatmap(
    heatmap_values,
    color = grDevices::colorRampPalette(heatmap_colors[[direction]])(100),
    breaks = seq(0, 0.1, length.out = 101),
    legend_breaks = c(0, 0.05, 0.1),
    clustering_method = "ward.D2",
    display_numbers = stars,
    fontsize_number = 14,
    angle_col = "315",
    main = paste(direction, "modules | Jaccard index"),
    filename = file.path(plots, paste0("wgcna_jaccard_", tolower(direction), ".pdf")),
    width = 8,
    height = 4
  )
}
message("Plots saved to ", normalizePath(plots))
# WGCNA soft-threshold selection ----
data_dir <- "data"
plots <- "plot"
sft_files <- c(TCGA = "tcga_sft_count.RData", GEO = "geo_sft_count.RData")
sft_objects <- c(TCGA = "sft2", GEO = "sft")
sft_data <- lapply(names(sft_files), function(cohort) {
  input <- new.env(parent = emptyenv())
  load(file.path(data_dir, sft_files[[cohort]]), envir = input)
  fit <- input[[sft_objects[[cohort]]]]$fitIndices
  fit <- fit[match(1:20, fit$Power), c("Power", "SFT.R.sq", "slope", "mean.k.")]
  stopifnot(all(is.finite(as.matrix(fit))))
  fit$signed_r2 <- -sign(fit$slope) * fit$SFT.R.sq
  fit
})
names(sft_data) <- names(sft_files)
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
soft_threshold_pdf <- file.path(plots, "wgcna_soft_threshold.pdf")
grDevices::pdf(soft_threshold_pdf, width = 8, height = 7, bg = "white")
tryCatch({
  par(mfrow = c(2, 2), mar = c(4.1, 4.5, 2.5, 1))
  for (cohort in names(sft_data)) {
    fit <- sft_data[[cohort]]
    plot(fit$Power, fit$signed_r2, type = "n",
         ylim = range(c(fit$signed_r2, 0.8)),
         xlab = "Soft-thresholding power", ylab = expression("Signed " * R^2),
         main = paste(cohort, "cohort"))
    abline(h = 0.8, col = "red")
    text(fit$Power, fit$signed_r2, labels = fit$Power, col = "red", cex = 0.8)
    plot(fit$Power, fit$mean.k., type = "n",
         xlab = "Soft-thresholding power", ylab = "Mean connectivity")
    text(fit$Power, fit$mean.k., labels = fit$Power, col = "red", cex = 0.8)
  }
}, finally = grDevices::dev.off())
message("Soft-threshold plot saved to ", normalizePath(soft_threshold_pdf))
# Archived TCGA network construction ----
if (FALSE) {
prad_legacy_root <- getwd()
rm(list=ls())
wd <- file.path(prad_legacy_root, 'data/legacy_project/part3')
setwd(wd)
{
  
library(dplyr)
library('DESeq2')
library('stringr')
library(DESeq2)
library(magrittr)
library(WGCNA)
library(ggplot2)
load(file.path(prad_legacy_root, 'data/legacy_project/part1/1 >1 express gene.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/part2 network/bg.RData'))
rm(result_count_normal,result_tpm_normal,result_tpm_prad)
options(stringsAsFactors = FALSE)
}
#rowname=sample colnames = gene
{
df <- result_count_prad[bg,]
vst <- vst(df,blind = F,fitType = "parametric")
datExpr0 <- t(vst)
rm(df,result_count_prad,vst)
}
{
load(file.path(prad_legacy_root, 'data/legacy_project/part3/-ed/m17/rdata/geo_sft_count.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/part3/-ed/m17/rdata/tcga_sft_count.RData'))
powers = c(c(1:20))
sft = pickSoftThreshold(datExpr0,
                        powerVector = powers, 
                        verbose = 5, 
                        blockSize = 21000, 
                        networkType = "signed hybrid")
#signed12 hybird 7 unsigned 4
png("result/sft/sft_geo1.png",units="in", width=5, height=5,res=600)
plot(sft$fitIndices[,1][1:20], -sign(sft$fitIndices[,3])[1:20]*sft$fitIndices[,2][1:20],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1][1:20], -sign(sft$fitIndices[,3])[1:20]*sft$fitIndices[,2][1:20],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
dev.off() 
png("result/sft/sft_geo2.png",units="in", width=5, height=5,res=600)
plot(sft$fitIndices[,1][1:20], sft$fitIndices[,5][1:20],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1][1:20], sft$fitIndices[,5][1:20], labels=powers, cex=cex1,col="red")
dev.off() 
png("result/sft/sft_tcga1.png",units="in", width=5, height=5,res=600)
plot(sft2$fitIndices[,1][1:20], -sign(sft2$fitIndices[,3])[1:20]*sft2$fitIndices[,2][1:20],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft2$fitIndices[,1][1:20], -sign(sft2$fitIndices[,3])[1:20]*sft2$fitIndices[,2][1:20],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
dev.off() 
png("result/sft/sft_tcga2.png",units="in", width=5, height=5,res=600)
plot(sft2$fitIndices[,1][1:20], sft2$fitIndices[,5][1:20],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft2$fitIndices[,1][1:20], sft2$fitIndices[,5][1:20], labels=powers, cex=cex1,col="red")
dev.off() }
# Adjacency and topological overlap ----
{
# https://www.rdocumentation.org/packages/WGCNA/versions/1.71/topics/adjacency
adjacency = adjacency(datExpr0, power = 12, type = "signed")
# https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/TOMsimilarity
TOM_signed_full = TOMsimilarity(adjacency,
                                     TOMType = "signed")
dissTOM = 1-TOM_signed_full
}
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), 
# "median" (= WPGMC) or "centroid" (= UPGMC).
                                    geneTree = hclust(as.dist(dissTOM), method = "mcquitty");
png("result/sft/gentree—tcga.png",units="in", width=9, height=5,res=600)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
dev.off() 
{
                                                                minModuleSize = 80;
dynamicMods = cutreeDynamic(
                            dendro = geneTree,
                                                                      cutHeight = 0.999,
                            minClusterSize = minModuleSize,
                                                                       method = "hybrid", # hybrid ,tree
                            distM = dissTOM,
                                                                    deepSplit = 2, #1-4
                                                                pamRespectsDendro = T,
                            pamStage = TRUE,
                            useMedoids = FALSE,
                            maxDistToLabel = NULL,
                            # maxPamDist = cutHeight,
                            respectSmallClusters = TRUE,
                            verbose = 2, indent = 0)
table(dynamicMods)
}
dynamicColors = labels2colors(dynamicMods)
{
# table(dynamicColors)
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
#                     dendroLabels = FALSE, 
#                     hang = 0.03,
#                     addGuide = TRUE, 
#                     guideHang = 0.05, 
#                     main = "Gene dendrogram and module colors")
}
{MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
shoppingMETree = hclust(as.dist(MEDiss), method = "average")}#mcquitty,average
{
# plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# MEDissThres = 0.3
# abline(h=MEDissThres, col = "red")
                                                                    MEDissThres = 0.15
merge = mergeCloseModules(datExpr0, dynamicColors, 
                          cutHeight = MEDissThres, verbose = 1)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
# png("result/sft/gentree—tcga.png",units="in", width=9, height=5,res=600)
#  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
#                      dendroLabels = FALSE, hang = 0.03, addGuide = T, guideHang = 0.05)
#  dev.off() 
moduleColors = mergedColors
} # gai
table(moduleColors)
length(table(moduleColors))
{colorOrder = c("grey", standardColors(50));
moduleLabels = data.frame(row.names = unique(moduleColors), 
                          cluster_number = paste0("Cluster-", 1:length(unique(unique(moduleColors)))),
                          data_number = paste0("TCGA-", 1:length(unique(unique(moduleColors)))))
moduleLabels["grey",] = "Unassigned"
MEs = mergedMEs;
wgcna_tcga_results <-  list()
wgcna_tcga_results$full_MEs = MEs
wgcna_tcga_results$full_module_unmerged = dynamicColors
wgcna_tcga_results$full_module_merged = moduleColors
wgcna_tcga_results$full_module_labels = moduleLabels
wgcna_tcga_results$full_gene_tree = geneTree
wgcna_tcga_results$genes = colnames(datExpr0)}
save(wgcna_tcga_results,file = 'rdata/wgcan_tcga_result.RData')
{
MEnames <- colnames(wgcna_tcga_results$full_MEs)
class <- wgcna_tcga_results$full_module_labels
d2 <- rownames(class)
bindclass <- cbind(class,d2)
v <- as.data.frame(wgcna_tcga_results$full_module_merged)
colnames(v) <- 'color'
v1 <- left_join(v, bindclass, by=c('color'='d2'))
wgcna_tcga_results$cluster_n <- v1$cluster_number
wgcna_tcga_results$data_n <- v1$data_number
wgcna_tcga_results$bindclass <- bindclass
e <- as.data.frame(cbind(wgcna_tcga_results$genes,wgcna_tcga_results$full_module_merged))
colnames(e) <- c('gene','color')
list_tcga_prad <- list()
for (i in d2) {
  #name <- paste('tcga_me',i,sep="")
  a <- e %>%
    dplyr::filter(color == i)
  #list1 <- list(list1,paste('tcga_me',i,sep=""))
  list_tcga_prad[[paste("ME",i,sep = "" )]] <- a
}
list_tcga_prad[["MEgrey"]] <- NULL
# save(list_tcga_prad,file = 'rdata/list_tcga_prad.RData')
class <- dplyr::filter(class, !cluster_number=='Unassigned')
names(list_tcga_prad) <- class$data_number
list_tcga_prad_num <- list_tcga_prad
length(list_tcga_prad_num)
}
save(list_tcga_prad_num,file = 'rdata/list_tcga_prad_num.RData')
# Module membership ----
{nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
modNames = substring(names(wgcna_tcga_results$full_MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, wgcna_tcga_results$full_MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
namemms <- left_join(as.data.frame(modNames), bindclass, by=c('modNames'='d2'))
names(geneModuleMembership) = namemms$data_number
names(MMPvalue) = namemms$data_number
dim(geneModuleMembership)
wgcna_tcga_results$MM <- geneModuleMembership
wgcna_tcga_results$MMpvalue <- MMPvalue
}
save(wgcna_tcga_results,file = 'rdata/wgcna_tcga_results.RData')
load('rdata/wgcna_tcga_results.RData')
{
for(i in names(list_tcga_prad_num)){
  
      MM <- tibble::rownames_to_column(wgcna_tcga_results[["MM"]] ,  var = 'gene')
      mm <- left_join(list_tcga_prad_num[[i]],MM,by=('gene'='gene'))
      mm <- mm[,c('gene','color',i)]
      colnames(mm) <- c('gene','color','MM')
      list_tcga_prad_num[[i]] <- mm
}
 
}
save(list_tcga_prad_num,file = 'rdata/list_tcga_prad_num.RData')
# Intramodular connectivity ----
# connectivity_tcga <- intramodularConnectivity(adjacency,wgcna_tcga_results$data_n,scaleByMax = FALSE)
# 
# NC_WGCNA_tcga <- fundamentalNetworkConcepts(adjacency, GS = NULL)
# save(NC_WGCNA_tcga, file = 'rdata/NC_WGCNA_tcga.RData')
load(file.path(prad_legacy_root, 'data/legacy_project/part3/-ed/m12/rdata/NC_WGCNA_tcga.RData'))#signed
for (i in names(list_tcga_prad_num)) {
  Connectivity <- NULL
  ScaledConnectivity <- NULL 
  ClusterCoef <- NULL
  MAR <- NULL
  TCGA_NUM <- NULL
  for (j in list_tcga_prad_num[[i]]$gene) {
    con <- NC_WGCNA_tcga[["Connectivity"]][[j]]
    Connectivity <- c(Connectivity,con)
    sc <- NC_WGCNA_tcga[["ScaledConnectivity"]][[j]]
    ScaledConnectivity <- c(ScaledConnectivity,sc)
    cl <- NC_WGCNA_tcga[["ClusterCoef"]][[j]]
    ClusterCoef <- c(ClusterCoef,cl)
    ma <- NC_WGCNA_tcga[["MAR"]][[j]]
    MAR <- c(MAR,ma)
    tc <- i
    TCGA_NUM <- c(TCGA_NUM,i)
  }
  bind <- cbind(list_tcga_prad_num[[i]],Connectivity,ScaledConnectivity,ClusterCoef,MAR,TCGA_NUM)
  list_tcga_prad_num[[i]] <- bind
}
View(list_tcga_prad_num)
save(list_tcga_prad_num,file = 'rdata/list_tcga_prad_num.RData')
load('rdata/list_tcga_prad_num.RData')
# rm(adjacency,bind,bindclass,)
# rm(nGenes,nSamples,modNames,geneModuleMembership,MMPvalue,namemms,a,TOM_signed_full,
#    v,v1,moduleLabels,METree,MEs,mergedMEs,merge,MEList,MEDiss,e,dissTOM,class)
}
# Archived GEO network construction ----
if (FALSE) {
prad_legacy_root <- getwd()
rm(list=ls())
wd <- file.path(prad_legacy_root, 'data/legacy_project/part3')
setwd(wd)
{library(dplyr)
library('stringr')
library(DESeq2)
library(magrittr)
library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE)
load(file.path(prad_legacy_root, 'data/legacy_project/WGCNA/rdata/0 GEO original.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/part2 network/bg.RData'))
prad <- as.matrix(prad0)
df <- prad[bg,]
df1 <- apply(df,c(1,2),round)
vst <- vst(df1,blind = F,fitType = "parametric")
datExpr0 <- t(vst)
}
# sampleTree = hclust(dist(datExpr0), method = "average");
rm(df,df1,prad,vst,powers,prad0)
# Adjacency and topological overlap ----
{
  # https://www.rdocumentation.org/packages/WGCNA/versions/1.71/topics/adjacency
  adjacency = adjacency(datExpr0, power = 16, type = "signed")
  # https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/TOMsimilarity
  TOM_signed_full = TOMsimilarity(adjacency,
                                  TOMType = "signed")
  dissTOM = 1-TOM_signed_full
  
}
                                    geneTree = hclust(as.dist(dissTOM), method = "mcquitty");
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
{
                                                                    minModuleSize = 80;
  dynamicMods = cutreeDynamic(
    dendro = geneTree,
                                                                      cutHeight = 0.999,
    minClusterSize = minModuleSize,
                                                                        method = "hybrid", # hybrid ,tree
    distM = dissTOM,
                                                                        deepSplit = 2, #1-4
  
                                                                pamRespectsDendro = T,
    pamStage = TRUE,
    useMedoids = FALSE,
    maxDistToLabel = NULL,
    respectSmallClusters = TRUE,
    # maxPamDist = cutHeight,
    verbose = 2, indent = 0)
  table(dynamicMods)
}
dynamicColors = labels2colors(dynamicMods)
{
  # table(dynamicColors)
  # plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
  #                     dendroLabels = FALSE, 
  #                     hang = 0.03,
  #                     addGuide = TRUE, 
  #                     guideHang = 0.05, 
  #                     main = "Gene dendrogram and module colors")
}
{MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")}
{
  # plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  # MEDissThres = 0.3
  # abline(h=MEDissThres, col = "red")
  
                                                                     MEDissThres = 0.15
  merge = mergeCloseModules(datExpr0, dynamicColors, 
                            cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors;
 
   mergedMEs = merge$newMEs;
  
  
  png("result/sft/gentree—geo.png",units="in", width=9, height=5,res=600)
   plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
                       dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  dev.off() 
  
  
    moduleColors = mergedColors
} 
table(moduleColors)
length(table(moduleColors))
{
colorOrder = c("grey", standardColors(50));
moduleLabels = data.frame(row.names = unique(moduleColors), 
                          cluster_number = paste0("Cluster-", 1:length(unique(unique(moduleColors)))),
                          data_number = paste0("GEO-", 1:length(unique(unique(moduleColors)))))
moduleLabels["grey",] = "Unassigned"
MEs = mergedMEs;
wgcna_geo_results <-  list()
wgcna_geo_results$full_MEs = MEs
wgcna_geo_results$full_module_unmerged = dynamicColors
wgcna_geo_results$full_module_merged = moduleColors
wgcna_geo_results$full_module_labels = moduleLabels
wgcna_geo_results$full_gene_tree = geneTree
wgcna_geo_results$genes = colnames(datExpr0)
}
save(wgcna_geo_results,file = 'rdata/wgcan_geo_result.RData')
{
MEnames <- colnames(wgcna_geo_results$full_MEs)
class <- wgcna_geo_results$full_module_labels
d2 <- rownames(class)
bindclass <- cbind(class,d2)
v <- as.data.frame(wgcna_geo_results$full_module_merged)
colnames(v) <- 'color'
v1 <- left_join(v, bindclass, by=c('color'='d2'))
wgcna_geo_results$cluster_n <- v1$cluster_number
wgcna_geo_results$data_n <- v1$data_number
e <- as.data.frame(cbind(wgcna_geo_results$genes,wgcna_geo_results$full_module_merged))
colnames(e) <- c('gene','color')
list_geo_prad <- list()
for (i in d2) {
  #name <- paste('tcga_me',i,sep="")
  a <- e %>%
    dplyr::filter(color == i)
  #list1 <- list(list1,paste('tcga_me',i,sep=""))
  list_geo_prad[[paste("ME",i,sep = "" )]] <- a
}
list_geo_prad[["MEgrey"]] <- NULL
#save(list_geo_prad,file = 'rdata/list_geo_prad.RData')
# load('data/list_geo_prad.RData')
class <- dplyr::filter(class, !cluster_number=='Unassigned')
names(list_geo_prad) <- class$data_number
list_geo_prad_num <- list_geo_prad
}
save(list_geo_prad_num,file = 'rdata/list_geo_prad_num.RData')
# Module membership ----
{nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
modNames = substring(names(wgcna_geo_results$full_MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, wgcna_geo_results$full_MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
namemms <- left_join(as.data.frame(modNames), bindclass, by=c('modNames'='d2'))
names(geneModuleMembership) = namemms$data_number
names(MMPvalue) = namemms$data_number
wgcna_geo_results$MM <- geneModuleMembership
wgcna_geo_results$MMpvalue <- MMPvalue
}
save(wgcna_geo_results,file = 'rdata/wgcna_geo_results.RData')
{
for(i in names(list_geo_prad_num)){
  
  MM <- tibble::rownames_to_column(wgcna_geo_results[["MM"]] ,  var = 'gene')
  mm <- left_join(list_geo_prad_num[[i]],MM,by=('gene'='gene'))
  mm <- mm[,c('gene','color',i)]
  colnames(mm) <- c('gene','color','MM')
  list_geo_prad_num[[i]] <- mm
}
}
save(list_geo_prad_num,file = 'rdata/list_geo_prad_num.RData')
# Intramodular connectivity ----
# connectivity_tcga <- intramodularConnectivity(adjacency, wgcna_tcga_results$data_n, scaleByMax = FALSE)
# 
# NC_WGCNA_geo <- fundamentalNetworkConcepts(adjacency, GS = NULL)
# save(NC_WGCNA_geo, file = 'rdata/NC_WGCNA_geo.RData')
load(file.path(prad_legacy_root, 'data/legacy_project/part3/m12/rdata/NC_WGCNA_geo.RData'))#signed
for (i in names(list_geo_prad_num)) {
  Connectivity <- NULL
  ScaledConnectivity <- NULL 
  ClusterCoef <- NULL
  MAR <- NULL
  GEO_NUM <- NULL
  for (j in list_geo_prad_num[[i]]$gene) {
    con <- NC_WGCNA_geo[["Connectivity"]][[j]]
    Connectivity <- c(Connectivity,con)
    sc <- NC_WGCNA_geo[["ScaledConnectivity"]][[j]]
    ScaledConnectivity <- c(ScaledConnectivity,sc)
    cl <- NC_WGCNA_geo[["ClusterCoef"]][[j]]
    ClusterCoef <- c(ClusterCoef,cl)
    ma <- NC_WGCNA_geo[["MAR"]][[j]]
    MAR <- c(MAR,ma)
    ge <- i
    GEO_NUM <- c(GEO_NUM,i)
  }
  bind <- cbind(list_geo_prad_num[[i]],Connectivity,ScaledConnectivity,ClusterCoef,MAR,GEO_NUM)
  list_geo_prad_num[[i]] <- bind
}
View(list_geo_prad_num)
save(list_geo_prad_num,file = 'rdata/list_geo_prad_num.RData')
 # rm(v,v1,moduleLabels,METree,MEs,mergedMEs,merge,MEList,
 # MEDiss,e,dissTOM,class,bindclass,a)
}
# Archived module overlap ----
if (FALSE) {
prad_legacy_root <- getwd()
{rm(list=ls())
wd <- file.path(prad_legacy_root, 'data/legacy_project/part3')
setwd(wd)
library(dplyr)
library(tibble)
library(stringr)
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/0.01overlap_18tcga_up_down.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/part2 network/bg.RData'))
}
two_gene_list_overlap<-function(gene_list_1,gene_list_2,bg_num){
  #note: gene_list_1 and gene_list_2 are come from bg_gene 
  over_gene<-intersect(gene_list_1,gene_list_2)
  p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
  result<-cbind(length(gene_list_1),length(gene_list_2),length(over_gene),p_hyper)
  return(result)
}
load('rdata/list_tcga_prad_num.RData')
load('rdata/list_geo_prad_num.RData')
names(list_tcga_prad_num)
names(list_geo_prad_num)
{overlap_tcga_up <- NULL
name_tcga <- NULL
for(i in names(list_tcga_prad_num)){
  
  gene_list_1 <- as.matrix(list_tcga_prad_num[[i]][1])
  gene_list_2 <- over_up
  o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
  
  if (o[,4]<0.05) {
    colname <- c(i,'deg_up',"overlaps","p_hyper")
    o <- rbind(colname,o)
    overlap_tcga_up <- rbind(overlap_tcga_up,o)
    name_tcga <- c(name_tcga,i)
  }
}
 View(overlap_tcga_up)
 length(name_tcga)}
{overlap_geo_up <- NULL
name_geo <- NULL
for(i in names(list_geo_prad_num)){
  
  gene_list_1 <- as.matrix(list_geo_prad_num[[i]][1])
  gene_list_2 <- over_up
  o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
  colname <- c(i,'deg_up',"overlaps","p_hyper")
  if (o[,4]<0.05) {
    o <- rbind(colname,o)
    overlap_geo_up <- rbind(overlap_geo_up,o)
    name_geo <- c(name_geo,i)
  }
}
View(overlap_geo_up)
length(name_geo)}
# Upregulated module overlap ----
{list_tcga <- as.numeric(str_replace(name_tcga,"TCGA-",""))
list_geo <- as.numeric(str_replace(name_geo,"GEO-",""))
overlap_up <- NULL
overtcga_up <- NULL
overgeo_up <- NULL
rownametcga <- NULL
dfupup <- NULL
for (i in list_tcga) {
  colnamegeo <- NULL
  ctcga <- NULL #p
  for (j in list_geo) {
    data1 <- as.matrix(list_tcga_prad_num[[paste('TCGA-',i,sep="")]][1])
    data2 <- as.matrix(list_geo_prad_num[[paste('GEO-',j,sep="")]][1])
    gene_list_1 <- data1[,1]
    gene_list_2 <- data2[,1]
    
    two_gene_list_overlap<-function(gene_list_1,gene_list_2,bg_num){
      
      over_gene<-intersect(gene_list_1,gene_list_2)
      p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
      result<-cbind(length(gene_list_1),length(gene_list_2),length(over_gene),p_hyper)
      save(over_gene,file = paste("result/up/",i,"_",j,".RData",sep = ""))
      return(result)
    }
    
    o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
    colname <- c(paste('TCGA-',i,sep = ""),paste('GEO-',j,sep = ""),"overlaps","p_hyper")
    
    colnamegeo <- c(colnamegeo,paste('GEO_',j,sep=""))
    ctcga <- cbind(ctcga,o[4])
    
    if (o[,4]<0.05) {
      o <- rbind(colname,o)
      overlap_up <- rbind(overlap_up,o)
      overtcga_up <- c(overtcga_up,i)
      overgeo_up <- c(overgeo_up,j) 
    }
  }
  dfupup <- rbind(dfupup,ctcga)
  rownametcga <- c(rownametcga,paste('TCGA_',i,sep=""))
  rownames(dfupup) <- rownametcga
  colnames(dfupup) <- colnamegeo  
}
View(overlap_up)
save(dfupup,file = 'rdata/heatmap_up_tcgageo.RData')
}
# Overlap with shared differentially expressed genes ----
{overlap_deg <- NULL
t_g <- NULL
for (i in paste(overtcga_up,'_',overgeo_up,sep = "")) {
  
  load(paste('result/up/',i,".RData",sep = ""))
  
  gene_list_1 <- over_gene
  gene_list_2 <- over_up
  
  two_gene_list_overlap<-function(gene_list_1,gene_list_2,bg_num){
    over_gene<-intersect(gene_list_1,gene_list_2)
    p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
    result<-cbind(length(gene_list_1),length(gene_list_2),length(over_gene),p_hyper)
    save(over_gene,file = paste("result/up/deg/",i,".RData",sep = ""))
    return(result)
  }
  
  o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
  colname <- c(paste('tcga_geo_',i,sep = ""),paste('deg_up',sep = ""),"overlaps","p_hyper")
  o <- rbind(colname,o)
  overlap_deg <- rbind(overlap_deg,o)
  i <- str_replace(i , '_', '')
  t_g <- as.numeric(c(t_g,i))
  
}
View(overlap_deg)
}
{
write.table(overlap_tcga_up,
            file="result/1 overlap_tcga_upn.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(overlap_geo_up,
            file="result/2 overlap_geo_upn.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(overlap_up,
            file="result/3 up_tcga_geo.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(overlap_deg,
            file="result/4 up_tcga_geo.txt",sep="\t",row.names=F,col.names=F,quote=F)
}
# Upregulated module Jaccard index ----
similarity_matrix_up = matrix(0, nrow = length(list_tcga), ncol = length(list_geo))
rownames(similarity_matrix_up) = paste('TCGA_',list_tcga,sep = "")
colnames(similarity_matrix_up) = paste('GEO_',list_geo,sep = "")
for (i in list_tcga) {
  for(j in list_geo){
    similarity_matrix_up[paste('TCGA_',i,sep = ""),paste('GEO_',j,sep = "")] = 
      length(intersect(list_tcga_prad_num[[paste('TCGA-',i,sep = "")]][['gene']],
                       list_geo_prad_num[[paste('GEO-',j,sep = "")]][['gene']]))/
      length(union(list_tcga_prad_num[[paste('TCGA-',i,sep = "")]][['gene']],
                   list_geo_prad_num[[paste('GEO-',j,sep = "")]][['gene']]))
   
  }
  
}
write.csv(similarity_matrix_up,file='result/11.1_jaccard_upup.csv')
save(similarity_matrix_up,file='rdata/11.1_jaccard_upup.RData')
{overlap_tcga_d <- NULL
  name_tcga <- NULL
  for(i in names(list_tcga_prad_num)){
    
    gene_list_1 <- as.matrix(list_tcga_prad_num[[i]][1])
    gene_list_2 <- over_down
    o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
    
    if (o[,4]<0.05) {
      colname <- c(i,'deg_down',"overlaps","p_hyper")
      o <- rbind(colname,o)
      overlap_tcga_d <- rbind(overlap_tcga_d,o)
      name_tcga <- c(name_tcga,i)
    }
  }
  View(overlap_tcga_d)
  length(name_tcga)}
{overlap_geo_d <- NULL
name_geo <- NULL
for(i in names(list_geo_prad_num)){
  
  gene_list_1 <- as.matrix(list_geo_prad_num[[i]][1])
  gene_list_2 <- over_down
  o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
  colname <- c(i,'deg_down',"overlaps","p_hyper")
  if (o[,4]<0.05) {
    o <- rbind(colname,o)
    overlap_geo_d <- rbind(overlap_geo_d,o)
    name_geo <- c(name_geo,i)
  }
}
View(overlap_geo_d)
length(name_geo)}
# Downregulated module overlap ----
{list_tcga <- as.numeric(str_replace(name_tcga,"TCGA-",""))
list_geo <- as.numeric(str_replace(name_geo,"GEO-",""))
overlap_down <- NULL
overtcga_down <- NULL
overgeo_down <- NULL
rownametcga <- NULL
dfdown <- NULL
for (i in list_tcga) {
  colnamegeo <- NULL
  ctcga <- NULL #p
  for (j in list_geo) {
    data1 <- as.matrix(list_tcga_prad_num[[paste('TCGA-',i,sep="")]][1])
    data2 <- as.matrix(list_geo_prad_num[[paste('GEO-',j,sep="")]][1])
    gene_list_1 <- data1[,1]
    gene_list_2 <- data2[,1]
    
    two_gene_list_overlap<-function(gene_list_1,gene_list_2,bg_num){
      
      over_gene<-intersect(gene_list_1,gene_list_2)
      p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
      result<-cbind(length(gene_list_1),length(gene_list_2),length(over_gene),p_hyper)
      save(over_gene,file = paste("result/down/",i,"_",j,".RData",sep = ""))
      return(result)
    }
    
    o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
    colname <- c(paste('TCGA-',i,sep = ""),paste('GEO-',j,sep = ""),"overlaps","p_hyper")
    
    colnamegeo <- c(colnamegeo,paste('GEO_',j,sep=""))
    ctcga <- cbind(ctcga,o[4])
    
    if (o[,4]<0.05) {
      o <- rbind(colname,o)
      overlap_down <- rbind(overlap_down,o)
      overtcga_down <- c(overtcga_down,i)
      overgeo_down <- c(overgeo_down,j) 
    }
  }
  dfdown <- rbind(dfdown,ctcga)
  rownametcga <- c(rownametcga,paste('TCGA_',i,sep=""))
  rownames(dfdown) <- rownametcga
  colnames(dfdown) <- colnamegeo  
}
View(overlap_down)
save(dfdown,file = 'rdata/heatmap_down_tcgageo.RData')
}
# Overlap with shared differentially expressed genes ----
{overlap_deg_d <- NULL
t_g <- NULL
for (i in paste(overtcga_down,'_',overgeo_down,sep = "")) {
  
  load(paste('result/down/',i,".RData",sep = ""))
  
  gene_list_1 <- over_gene
  gene_list_2 <- over_down
  
  two_gene_list_overlap<-function(gene_list_1,gene_list_2,bg_num){
    over_gene<-intersect(gene_list_1,gene_list_2)
    p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
    result<-cbind(length(gene_list_1),length(gene_list_2),length(over_gene),p_hyper)
    save(over_gene,file = paste("result/down/deg/",i,".RData",sep = ""))
    return(result)
  }
  
  o <- two_gene_list_overlap(gene_list_1,gene_list_2,bg_num)
  colname <- c(paste('tcga_geo_',i,sep = ""),paste('deg_down',sep = ""),"overlaps","p_hyper")
  o <- rbind(colname,o)
  overlap_deg_d <- rbind(overlap_deg_d,o)
  i <- str_replace(i , '_', '')
  t_g <- as.numeric(c(t_g,i))
  
}
View(overlap_deg_d)
}
{
  write.table(overlap_tcga_d,
              file="result/5 overlap_tcga_down.txt",sep="\t",row.names=F,col.names=F,quote=F)
  
  write.table(overlap_geo_d,
              file="result/6 overlap_geo_down.txt",sep="\t",row.names=F,col.names=F,quote=F)
  
  write.table(overlap_down,
              file="result/7 down_tcga_geo.txt",sep="\t",row.names=F,col.names=F,quote=F)
  
  write.table(overlap_deg_d,
              file="result/8 down_tcga_geo.txt",sep="\t",row.names=F,col.names=F,quote=F)
}
# Downregulated module Jaccard index ----
similarity_matrix_down = matrix(0, nrow = length(list_tcga), ncol = length(list_geo))
rownames(similarity_matrix_down) = paste('TCGA_',list_tcga,sep = "")
colnames(similarity_matrix_down) = paste('GEO_',list_geo,sep = "")
for (i in list_tcga) {
  for(j in list_geo){
    similarity_matrix_down[paste('TCGA_',i,sep = ""),paste('GEO_',j,sep = "")] = 
      length(intersect(list_tcga_prad_num[[paste('TCGA-',i,sep = "")]][['gene']],
                       list_geo_prad_num[[paste('GEO-',j,sep = "")]][['gene']]))/
      length(union(list_tcga_prad_num[[paste('TCGA-',i,sep = "")]][['gene']],
                   list_geo_prad_num[[paste('GEO-',j,sep = "")]][['gene']]))
    
  }
  
}
write.csv(similarity_matrix_down,file='result/11.2_jaccard_dd.csv')
save(similarity_matrix_down,file='rdata/11.2_jaccard_dd.RData')
}
# Archived module GO enrichment ----
if (FALSE) {
prad_legacy_root <- getwd()
{rm(list=ls())
wd <- file.path(prad_legacy_root, 'data/legacy_project/part3/')
setwd(wd)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library("BiocManager")
load('rdata/list_tcga_prad_num.RData')
load('rdata/list_geo_prad_num.RData')
load(file.path(prad_legacy_root, 'data/legacy_project/part2 network/bg.RData'))
}
## function
gofunction <- function(gene, bg){
  
  OrgDb <- 'org.Hs.eg.db'
  go1 <- enrichGO(gene=gene,
                  OrgDb='org.Hs.eg.db',
                  universe = bg,
                  keyType = "SYMBOL", 
                  ont="BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff=0.05,
                  pvalueCutoff=0.05)
  
  go_result <- go1@result
  go_result <- filter(go_result, pvalue < 0.05)
  return(go_result)
}
##tcga all go
go_tcga_result <- list()
for (i in names(list_tcga_prad_num)) {
  setwd(wd)
  gene <- as.vector(t(list_tcga_prad_num[[i]][1]))
  go_result <- gofunction(gene,bg)
  go_tcga_result[[i]] <- go_result
  write.table(go_result,
              file = paste('result/tcga_go/',i,'.txt',sep = ""),sep="\t",row.names=F,col.names=T,quote=F)
}
## geo all go
go_geo_result <- list()
for (i in names(list_geo_prad_num)) {
  setwd(wd)
  gene <- as.vector(t(list_geo_prad_num[[i]][1]))
  go_result <- gofunction(gene,bg)
  go_tcga_result[[i]] <- go_result
  write.table(go_result,
              file = paste('result/tcga_go/',i,'.txt',sep = ""),sep="\t",row.names=F,col.names=T,quote=F)
}
save(go_tcga_result,go_geo_result,file = 'rdata/4 go_result.RData')
}
# Archived upregulated module topology ----
if (FALSE) {
prad_legacy_root <- getwd()
{rm(list=ls())
rm(list=ls())
rm(list=ls())
setwd(file.path(prad_legacy_root, 'data/legacy_project/part3'))
{require(matrixStats)
require(data.table)
require(igraph)
require(Hmisc)
require(compiler)
require(reshape2)
require(gridExtra)
require(ggplot2)  
require(ggpubr)  
require(ggpmisc)  
require(ggrepel)  
  library(tibble)  
  }
load("rdata/list_geo_prad_num.RData")
load("rdata/list_tcga_prad_num.RData")
}
t_g <-t_g
bind_result <- list()
for (i in overtcga_up) {
  for (j in overgeo_up) {
    
    if (sum(duplicated(c(t_g,paste(i,j,sep = ""))))==1) {
      
        int_gene_tcga <- as.data.frame(list_tcga_prad_num[[paste('TCGA-',i,sep="")]])
        # list_name_tcga <- as.data.frame(c(list_name_tcga,paste('TCGA-',i,sep="")))
        int_gene_geo <- as.data.frame(list_geo_prad_num[[paste('GEO-',j,sep="")]])
        # list_name_geo <- as.data.frame(c(list_name_geo,paste('GEO-',j,sep="")))
        
        load(paste('result/up/deg/',paste(i,'_',j,sep = ""),".RData",sep = ""))
        overlap_gene <- over_gene
        
        dir.create(paste('result/',i,'_',j,'/',sep = ""))#new folder
        
      # Hub gene selection
      
        tcga <- arrange(int_gene_tcga, desc(Connectivity))[1:ceiling(dim(int_gene_tcga)[1]*0.15),]
        geo <- arrange(int_gene_geo, desc(Connectivity))[1:ceiling(dim(int_gene_geo)[1]*0.15),]
        bind <- left_join(tcga,geo, by=c('gene'='gene'))
        bind_deg <-left_join(as.data.frame(overlap_gene),bind,by=c('overlap_gene'='gene'))
        bind_deg <- na.omit(bind_deg)
        write.table(bind_deg,
                    file=paste('result/',i,'_',j,'/',i,j,'result_con.txt',sep = ""),sep="\t",row.names=F,col.names=T,quote=F)
        # View(bind_deg)
        bind_result[[paste(i,j,'0.15',sep = "_")]] <- bind_deg
        
   # >15%
      
        # # >0.8 
        tcga <- filter(int_gene_tcga, MM>0.8)
        geo <- filter(int_gene_geo, MM>0.8)
        bind <- left_join(tcga,geo, by=c('gene'='gene'))
        bind_deg2 <-left_join(as.data.frame(overlap_gene),bind,by=c('overlap_gene'='gene'))
        bind_deg2 <- na.omit(bind_deg2)
        # View(bind_deg2)
        write.table(bind_deg2,
                    file=paste('result/',i,'_',j,'/',i,j,'result_mm.txt',sep = ""),sep="\t",row.names=F,col.names=T,quote=F)
        bind_result[[paste(i,j,'mm',sep = "_")]] <- bind_deg2
       #0.8
      
      
    }
    
  }
  
}
View(bind_result)
# Module membership and connectivity data ----
load("rdata/list_geo_prad_num.RData")
load("rdata/list_tcga_prad_num.RData")
i <- 3
j <- 2
bind_result <- list()
    
      
      int_gene_tcga <- as.data.frame(list_tcga_prad_num[[paste('TCGA-',i,sep="")]])
      # list_name_tcga <- as.data.frame(c(list_name_tcga,paste('TCGA-',i,sep="")))
      int_gene_geo <- as.data.frame(list_geo_prad_num[[paste('GEO-',j,sep="")]])
      # list_name_geo <- as.data.frame(c(list_name_geo,paste('GEO-',j,sep="")))
      
      load(paste('result/up/deg/',paste(i,'_',j,sep = ""),".RData",sep = ""))
      overlap_gene <- over_gene
      
      dir.create(paste('result/',i,'_',j,'/',sep = ""))#new folder
      
      bind <- list()
      bind$bind_tcga_geo_mmall <- left_join(int_gene_tcga,int_gene_geo, by=c('gene'='gene'))
      bind$bind_tcga_geo_mmall <- na.omit(bind$bind_tcga_geo_mmall)
      bind$bind_tcga_geo_deg <-left_join(as.data.frame(overlap_gene),bind$bind_tcga_geo_mmall,by=c('overlap_gene'='gene'))
      bind$bind_plot_deg <- bind$bind_tcga_geo_deg[,c('overlap_gene','MM.x','MM.y')]  
      colnames(bind$bind_plot_deg) <- c('gene','TCGA MM','GEO MM')
      bind$bind_plot_deg2 <- bind$bind_tcga_geo_deg[,c('overlap_gene','Connectivity.x','Connectivity.y')]  
      colnames(bind$bind_plot_deg2) <- c('gene','TCGA Connectivity','GEO Connectivity')
bind$bind_plot_deg$car <- ''
lab <- c(27,33,26,25,38,29,34,35)
bind$bind_plot_deg$car[lab] <- bind$bind_plot_deg$gene[lab]
bind$bind_plot_deg2$car <- bind$bind_plot_deg$car
save(bind,file = 'rdata/dot_plot.RData')
t3_0 <- as.data.frame(list_tcga_prad_num[[paste('TCGA-','3',sep="")]])
g2 <- as.data.frame(list_geo_prad_num[[paste('GEO-','2',sep="")]])
load(paste('result/up/deg/',paste('3','_','2',sep = ""),".RData",sep = ""))
overlap_gene <- as.data.frame(over_gene)
bind_t3_0 <-arrange(left_join(overlap_gene , t3_0 , by=c('over_gene'='gene')),desc(Connectivity))
bind_g2 <-arrange(left_join(overlap_gene , g2 , by=c('over_gene'='gene')),desc(Connectivity))
View(bind_t3_0)
View(bind_g2)
genet3g2 <- c(bind_t3_0$over_gene[1:15],bind_g2$over_gene[1:15])
a <- as.data.frame(genet3g2[duplicated(genet3g2)])
top10_3 <- left_join(a , bind_t3_0 , by=c('genet3g2[duplicated(genet3g2)]'='over_gene'))
top10_2 <- left_join(a , bind_g2 , by=c('genet3g2[duplicated(genet3g2)]'='over_gene'))
write.csv(top10_3,file = 'result/13.1 gene sample info.txt')
write.csv(top10_2,file = 'result/13.2 gene sample info.txt')
write.table(top10_3,
            file="result/13.1 top8tcga-3 .txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(top10_2,
            file="result/13.2 top8geo-2.txt",sep="\t",row.names=F,col.names=F,quote=F)
}
# Archived downregulated module topology ----
if (FALSE) {
prad_legacy_root <- getwd()
{rm(list=ls())
  rm(list=ls())
  rm(list=ls())
  
  setwd(file.path(prad_legacy_root, 'data/legacy_project/part3/this！'))
  
  {require(matrixStats)
    require(data.table)
    require(igraph)
    require(Hmisc)
    require(compiler)
    require(reshape2)
    require(gridExtra)
    require(ggplot2)  
    require(ggpubr)  
    require(ggpmisc)  
    require(ggrepel)  
    
  }
  
  load("rdata/list_geo_prad_num.RData")
  load("rdata/list_tcga_prad_num.RData")
}
t_g <-t_g
bind_result <- list()
for (i in overtcga_down) {
  
  for (j in overgeo_down) {
    
    if (sum(duplicated(c(t_g,paste(i,j,sep = ""))))==1) {
      
      int_gene_tcga <- as.data.frame(list_tcga_prad_num[[paste('TCGA-',i,sep="")]])
      # list_name_tcga <- as.data.frame(c(list_name_tcga,paste('TCGA-',i,sep="")))
      int_gene_geo <- as.data.frame(list_geo_prad_num[[paste('GEO-',j,sep="")]])
      # list_name_geo <- as.data.frame(c(list_name_geo,paste('GEO-',j,sep="")))
      
      load(paste('result/down/deg/',paste(i,'_',j,sep = ""),".RData",sep = ""))
      overlap_gene <- over_gene
      
      dir.create(paste('result/down/',i,'_',j,'/',sep = ""))#new folder
      
      # Hub gene selection
      
      tcga <- arrange(int_gene_tcga, desc(Connectivity))[1:ceiling(dim(int_gene_tcga)[1]*0.15),]
      geo <- arrange(int_gene_geo, desc(Connectivity))[1:ceiling(dim(int_gene_geo)[1]*0.15),]
      bind <- left_join(tcga,geo, by=c('gene'='gene'))
      bind_deg <-left_join(as.data.frame(overlap_gene),bind,by=c('overlap_gene'='gene'))
      bind_deg <- na.omit(bind_deg)
      write.table(bind_deg,
                  file=paste('result/down/',i,'_',j,'/',i,j,'result_con.txt',sep = ""),sep="\t",row.names=F,col.names=T,quote=F)
      # View(bind_deg)
      bind_result[[paste(i,j,'0.15',sep = "_")]] <- bind_deg
      
      # >15%
      
      # # >0.8 
      tcga <- filter(int_gene_tcga, MM>0.8)
      geo <- filter(int_gene_geo, MM>0.8)
      bind <- left_join(tcga,geo, by=c('gene'='gene'))
      bind_deg2 <-left_join(as.data.frame(overlap_gene),bind,by=c('overlap_gene'='gene'))
      bind_deg2 <- na.omit(bind_deg2)
      # View(bind_deg2)
      write.table(bind_deg2,
                  file=paste('result/down/',i,'_',j,'/',i,j,'result_mm.txt',sep = ""),sep="\t",row.names=F,col.names=T,quote=F)
      bind_result[[paste(i,j,'mm',sep = "_")]] <- bind_deg2
      #0.8
      
      
    }
    
  }
  
}
View(bind_result)
# Module membership and connectivity data ----
load("rdata/list_geo_prad_num.RData")
load("rdata/list_tcga_prad_num.RData")
i <- 3
j <- 2
bind_result <- list()
int_gene_tcga <- as.data.frame(list_tcga_prad_num[[paste('TCGA-',i,sep="")]])
# list_name_tcga <- as.data.frame(c(list_name_tcga,paste('TCGA-',i,sep="")))
int_gene_geo <- as.data.frame(list_geo_prad_num[[paste('GEO-',j,sep="")]])
# list_name_geo <- as.data.frame(c(list_name_geo,paste('GEO-',j,sep="")))
load(paste('result/up/deg/',paste(i,'_',j,sep = ""),".RData",sep = ""))
overlap_gene <- over_gene
dir.create(paste('result/',i,'_',j,'/',sep = ""))#new folder
bind <- list()
bind$bind_tcga_geo_mmall <- left_join(int_gene_tcga,int_gene_geo, by=c('gene'='gene'))
bind$bind_tcga_geo_mmall <- na.omit(bind$bind_tcga_geo_mmall)
bind$bind_tcga_geo_deg <-left_join(as.data.frame(overlap_gene),bind$bind_tcga_geo_mmall,by=c('overlap_gene'='gene'))
bind$bind_plot_deg <- bind$bind_tcga_geo_deg[,c('overlap_gene','MM.x','MM.y')]  
colnames(bind$bind_plot_deg) <- c('gene','TCGA MM','GEO MM')
bind$bind_plot_deg2 <- bind$bind_tcga_geo_deg[,c('overlap_gene','Connectivity.x','Connectivity.y')]  
colnames(bind$bind_plot_deg2) <- c('gene','TCGA Connectivity','GEO Connectivity')
bind$bind_plot_deg$car <- ''
lab <- c(27,33,26,25,38,29,34,35)
bind$bind_plot_deg$car[lab] <- bind$bind_plot_deg$gene[lab]
bind$bind_plot_deg2$car <- bind$bind_plot_deg$car
save(bind,file = 'rdata/dot_plot.RData')
# Module topology plots ----
load('rdata/dot_plot.RData')
png("plot/dot_plot_mm.png",units="in", width=5, height=4.7,res=1000)
bind[["bind_plot_deg"]]%>%
  ggplot(aes(x=`TCGA MM`, y=`GEO MM`,label =car)) +
  geom_point(color='grey33')+
  geom_smooth(method = lm, level=0.99, color = 'black',size=0.5)+
  stat_cor(method = "spearman",label.x=-0.13, label.y=1,size=5)+
  geom_text_repel( ) +      #(box.padding  =  0.5 , max.overlaps  =  Inf)+
  geom_point(aes(0.6438509,0.96214012),color='red',size=2)+
  geom_point(aes(0.5731463,0.92767447),color='red',size=2)+
  geom_point(aes(0.8510341,0.87271877),color='red',size=2)+
  geom_point(aes(0.8864174,0.90407879),color='red',size=2)+
  geom_point(aes(0.6302074,0.90667241),color='red',size=2)+
  geom_point(aes(0.7527888,0.81644535),color='red',size=2)+  
  geom_point(aes(0.5534081,0.84470254),color='red',size=2)+
  geom_point(aes(0.8247785,0.81882123 ),color='red',size=2)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# ggtitle("TCGA-GEO-DEG MM")
dev.off() 
png("plot/dot_plot_connec.png",units="in", width=5, height=4.7 ,res=1000)
bind[["bind_plot_deg2"]]%>%
  ggplot(aes(x=`TCGA Connectivity`, y=`GEO Connectivity`,label =car)) +
  geom_point(color='grey33')+
  geom_smooth(method = lm, level=0.99, color = 'black',size=0.5)+
  stat_cor(method = "spearman",label.x=0, label.y=180,size=5)+
  geom_text_repel( ) +      #(box.padding  =  0.5 , max.overlaps  =  Inf)+
  geom_point(aes(80.13510,164.270907),color='red',size=2)+
  geom_point(aes(293.94036,142.104790),color='red',size=2)+
  geom_point(aes(64.05526,135.637182),color='red',size=2)+
  geom_point(aes(244.64729,108.904253),color='red',size=2)+
  geom_point(aes(82.00974,107.043533),color='red',size=2)+
  geom_point(aes(145.74099,101.118161),color='red',size=2)+  
  geom_point(aes(55.98441,84.840533),color='red',size=2)+
  geom_point(aes(202.70795,72.902622 ),color='red',size=2)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off() 
t3_0 <- as.data.frame(list_tcga_prad_num[[paste('TCGA-','3',sep="")]])
g2 <- as.data.frame(list_geo_prad_num[[paste('GEO-','2',sep="")]])
load(paste('result/up/deg/',paste('3','_','2',sep = ""),".RData",sep = ""))
overlap_gene <- as.data.frame(over_gene)
bind_t3_0 <-arrange(left_join(overlap_gene , t3_0 , by=c('over_gene'='gene')),desc(Connectivity))
bind_g2 <-arrange(left_join(overlap_gene , g2 , by=c('over_gene'='gene')),desc(Connectivity))
View(bind_t3_0)
View(bind_g2)
genet3g2 <- c(bind_t3_0$over_gene[1:15],bind_g2$over_gene[1:15])
a <- as.data.frame(genet3g2[duplicated(genet3g2)])
top10_3 <- left_join(a , bind_t3_0 , by=c('genet3g2[duplicated(genet3g2)]'='over_gene'))
top10_2 <- bind_g2
png("plot/gene line box plot.png",units="in", width=3.8, height=3,res=1000)
bind_t3_0[1:15,] %>%
  mutate(over_gene = fct_reorder(over_gene, Connectivity)) %>%
  ggplot( aes(x=over_gene, y=Connectivity)) +
  geom_point()+
  theme_bw()+
  geom_point(data=top10_3 ,aes(x=genet3g2[duplicated(genet3g2)],y=Connectivity), colour='red')+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  xlab("TCGA-3 genes") +
  ggtitle("  ")
dev.off() 
}
