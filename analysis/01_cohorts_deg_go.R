# TCGA and GEO differential expression ----
data_dir <- "data"
results <- file.path("result", "deg")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
library(DESeq2)
for (cohort in c("TCGA", "GEO")) {
  input <- new.env(parent = emptyenv())
  load(file.path(data_dir, paste0(cohort, "_cohort.Rdata")), input)
  tumor <- as.matrix(input[[paste0(cohort, "_prad_count")]])
  normal <- as.matrix(input[[paste0(cohort, "_normal_count")]])
  stopifnot(identical(rownames(tumor), rownames(normal)))
  colnames(tumor) <- paste0("p_", colnames(tumor))
  colnames(normal) <- paste0("n_", colnames(normal))
  counts <- cbind(tumor, normal)
  if (cohort == "GEO") counts <- round(counts)
  coldata <- data.frame(type = factor(rep(c("prad", "normal"), c(ncol(tumor), ncol(normal)))),
                        row.names = colnames(counts))
  model <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~type)
  fit <- DESeq(model)
  result <- as.data.frame(DESeq2::results(fit, contrast = c("type", "prad", "normal"),
    lfcThreshold = 0, altHypothesis = "greaterAbs", alpha = 0.05))
  result$gene <- rownames(result)
  result$direction <- ifelse(result$log2FoldChange < 0, "down", "up")
  write.csv(result, file.path(results, paste0(tolower(cohort), "_differential_expression.csv")), row.names = FALSE)
  significant <- result[which(result$padj < 0.01), ]
  write.csv(significant, file.path(results, paste0(tolower(cohort), "_degs.csv")), row.names = FALSE)
  message(cohort, ": ", nrow(counts), " genes, ", ncol(tumor), " tumor, ",
          ncol(normal), " normal, ", nrow(significant), " DEGs at FDR < 0.01.")
}
message("Results saved to ", results)
# GO enrichment and Venn diagram ----
data_dir <- "data"
results <- file.path("result", "go")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
plots <- "plot"
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
themelegendnone <- ggplot2::theme(legend.position = "none")
theme4_backline <- ggplot2::theme(panel.grid = ggplot2::element_blank())
library(patchwork)
library(grid)
load(file.path(data_dir, "GEO_cohort.Rdata"))
load(file.path(data_dir, "TCGA_cohort.Rdata"))
background_genes <- union(rownames(TCGA_prad_tpm),rownames(GEO_prad_tpm))
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
shared_up_genes <- intersect(rownames(GEOup_dif_exp_01),rownames(TCGAup_dif_exp_01))
shared_down_genes <- intersect(rownames(GEOdown_dif_exp_01),rownames(TCGAdown_dif_exp_01))
enrich_go_bp <- function( gene , background_genes){  
  library(dplyr) 
  library(org.Hs.eg.db) 
  library(clusterProfiler) 
  library(DOSE)
  library(ggplot2) 
  library(RColorBrewer)
  OrgDb <- 'org.Hs.eg.db'
  enrichment <- enrichGO(gene=gene,
                    OrgDb='org.Hs.eg.db',
                    universe = background_genes,
                    keyType = "SYMBOL", 
                    ont="BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff=0.05,
                    pvalueCutoff=0.05
  )
  enrichment_terms <- enrichment@result
  enrichment_terms <- filter(enrichment_terms, pvalue < 0.05)
  return(enrichment_terms)
}
deg_sets_by_direction <- list()
deg_sets_by_direction[['up_regulated']] <- shared_up_genes
deg_sets_by_direction[['down_regulated']] <-  shared_down_genes
library(dplyr)
go_all <- lapply(names(deg_sets_by_direction), function(i) {
  input_genes <- deg_sets_by_direction[[i]]
  enriched_terms <- enrich_go_bp(
    input_genes,
    background_genes = background_genes
  )
  enriched_terms %>%
    dplyr::filter(Count > 5) %>%
    dplyr::slice_min(order_by = pvalue, n = 20) %>%
    dplyr::mutate(
      cate = i,
      neglog10p = -log10(pvalue)
    )
}) %>%
  bind_rows()
go_all <- go_all %>% mutate(term_coverage = Count / as.numeric(sub("/.*", "", BgRatio)))
go_all$Description <- tools::toTitleCase(tolower(go_all$Description))
write.csv(go_all, file.path(results, "go_bp_enrichment_terms.csv"), row.names = FALSE)
go_up   <- go_all %>% dplyr::filter(cate == "up_regulated")
go_down <- go_all %>% dplyr::filter(cate == "down_regulated")
plot_go_dot <- function(df, title = NULL) {
  df <- df %>%
    arrange(pvalue) %>%
    mutate(Description = factor(Description,levels = rev(Description)  ))
  ggplot(df, aes(x = term_coverage,y = Description,size = Count,color = pvalue)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(name = "Gene count",range = c(1, 5)) +
    scale_color_gradient(low = "#475272",high = "#B5252E",name = "p-value",trans = "reverse") +
    labs(x = "GO term coverage",y = NULL,title = title) +
    theme_classic(base_size = 12) + theme4_backline
}
p_up <- plot_go_dot(go_up, title = "Upregulated genes")
p_up
p_down <- plot_go_dot(go_down, title = "Downregulated genes")
p_down
go_plot <- (p_up + themelegendnone) + p_down + plot_layout(heights = c(1,1)) 
ggsave(file.path(plots, "go_bp_enrichment.pdf"),plot = go_plot, width =16, height =20, units = "cm")
library("VennDiagram")
venn_gene_sets <- list()
venn_gene_sets[['TCGA up-regulated']] <- rownames(TCGAup_dif_exp_01)
venn_gene_sets[['GEO up-regulated']] <- rownames(GEOup_dif_exp_01)
venn_gene_sets[['TCGA down-regulated']] <- rownames(TCGAdown_dif_exp_01)
venn_gene_sets[['GEO down-regulated']] <- rownames(GEOdown_dif_exp_01)
venn_plot <- venn.diagram(x= venn_gene_sets, 
                  filename =NULL, 
                  disable.logging = TRUE,
                  height = 5000, width = 5000,resolution =3000, imagetype="png", 
                  col ="black", 
                  fill =c("white","white","white","white"),
                  alpha = 0.5, 
                  cex = 0.45,
                  cat.cex = 0.45,
                  lwd = 0.5,
                  margin         = 0.1,
                  cat.lwd = 0.5,
                  main.cex       = 0.5
)
pdf(file.path(plots, "deg_overlap_venn.pdf"), width = 1.7, height = 1.5)
grid.newpage()
grid.draw(venn_plot)
dev.off()
message("Results saved to ", results)
message("Plots saved to ", normalizePath(plots))
# Archived TCGA expression preprocessing ----
if (FALSE) {
prad_legacy_root <- getwd()
tpm_symbol_exp <- symbol_exp
save(tpm_symbol_exp,file=file.path(prad_legacy_root, 'data/legacy_project/tpm_normal.RData'))
tpm_prad_symbol_exp <- symbol_exp
save(tpm_prad_symbol_exp,file=file.path(prad_legacy_root, 'data/legacy_project/tpm_prad.RData'))
library(dplyr)
load("original.RData")
tpm_prad_symbol_exp <- as.data.frame(tpm_prad_symbol_exp)
df2 <- mutate(tpm_prad_symbol_exp, meanvalue = apply(tpm_prad_symbol_exp,1,mean)) # add new col 'meanvalue'
df12 <- filter(df2, meanvalue < 1)
v1 <- setdiff(rownames(tpm_symbol_exp),rownames(df12))
result_tpm_normal <- tpm_symbol_exp[v1,]
result_tpm_prad <- tpm_prad_symbol_exp[v1,]
result_count_normal <- no_symbol_exp[v1,]
result_count_prad <- prad_symbol_exp[v1,]
save(tpm_symbol_exp,tpm_prad_symbol_exp,no_symbol_exp,prad_symbol_exp, file = "original.RData")
save(result_tpm_normal,result_tpm_prad,result_count_normal,result_count_prad, file = "a_result.RData")
}
# Archived TCGA differential expression ----
if (FALSE) {
load("a_result.RData")
library('DESeq2')
library('stringr')
cn <- result_count_normal
cp <- result_count_prad
dimnames(cn) <-list(rownames(cn),t(as.list(str_c('n_', colnames(cn), sep=""))))# add colname'n_','p_'
cn <- as.matrix(cn)
dimnames(cp) <-list(rownames(cp),t(as.list(str_c('p_', colnames(cp), sep=""))))
cp <- as.matrix(cp)
data <- as.matrix(cbind(cp,cn))
group_type<-as.matrix(rbind(matrix("prad",dim(cp)[2],1),
                            matrix("normal",dim(cn)[2],1)))
colnames(group_type)<-"type"
count_data<-DESeqDataSetFromMatrix(countData=data,colData=group_type,design=~type)
deSeqAnalysis<-DESeq(count_data) 
res = results(deSeqAnalysis, 
              contrast=c("type","prad","normal"), 
              lfcThreshold = 0,
              altHypothesis="greaterAbs",
              alpha = 0.05) 
res1 <- as.data.frame(res)
res1 <- res1[which(res1 $ padj < 0.01),]
res1 <- mutate(res1, up_down = if_else(log2FoldChange < 0 , 'down', 'up')) #down4206 up4984 
table(res1$up_down)
dim(res1)
up_dif_exp <- res1[which(res1$up_down=='up'),]
down_dif_exp <- res1[which(res1$up_down=='down'),]
save(up_dif_exp, down_dif_exp, res1, file="dif_exp_gene_up_down.RData")
}
# Archived GEO expression preprocessing ----
if (FALSE) {
prad_legacy_root <- getwd()
library(dplyr)
library(survival)
library(tibble)
g <- file.path(prad_legacy_root, 'data/legacy_project/GEO/')
g5out <- file.path(prad_legacy_root, 'data/legacy_project/GEO/output/')
setwd(g)
# Gene identifier mapping ----
# cli<-as.data.frame(read.csv("GSE183019_processed_TPM.txt",header=F,sep="\t"))[,-2]
# cli <- t(cli)
cli[1,] <- str_replace(cli[1,], '!Sample_', '')
colnames(cli) <- cli[1,]
cli1<- as.data.frame(cli[-1,])
cli1$source_name_ch1 <- str_replace(cli1$source_name_ch1, 'prostate tissue normal sample from patient ', '')
cli2 <- cli1[,c(8,18,19,20)]
cli2$source_name_ch1 <- str_replace(cli2$source_name_ch1, 'buffy coat normal sample from patient ', '')
cli2$source_name_ch1 <- str_replace(cli2$source_name_ch1, 'prostate tissue cancer sample from patient ', '')
cli2$characteristics_ch1 <- str_replace(cli2$characteristics_ch1, '5-year disease status: ', '')
cli2$characteristics_ch1.1 <- str_replace(cli2$characteristics_ch1.1, 'sample type: ', '')
cli2$characteristics_ch1.2 <- str_replace(cli2$characteristics_ch1.2, 'tissue: ', '')
save(cli2,file = file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/clinicaldata.RData')) 
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/clinicaldata.RData'))
cli3 <- filter(cli2,characteristics_ch1.2=='prostate tissue')
cli3 <- cli3[,-c(2,4)]
colnames(cli3) <- c('id','state')
data <- as.data.frame(t(read.csv("GSE183019_processed_TPM.txt",header=F,sep="\t")[,-2]))
colnames(data) <- data[1,]
data2 <- (data[-1,])
bind <- left_join(cli3, data2, by=c('id'='GENE'))
bind2 <- as.data.frame(apply(bind,2,as.numeric))
rownames(bind2) <- bind$id
bind2$state <- bind$state
bind3 <- bind2[,-c(1,2)]
normal <- as.data.frame(t(filter(bind2,state=='normal')[,-c(1,2)]))
prad <- as.data.frame(t(filter(bind2,state=='cancer')[,-c(1,2)]))
prad1 <- mutate(prad, meanvalue = apply(prad,1,mean))
prad2 <- (filter(prad1,meanvalue > 1))[,-dim(prad1)[2]]
exp_gene <- rownames(prad2)
bind_final <- bind3[,exp_gene]
normal_exp_tpm <- normal[exp_gene,]
prad_exp_tpm <- prad2
save(bind_final,normal_exp_tpm,prad_exp_tpm,file=file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/>1express_normal_prad.RData'))
}
# Archived GEO differential expression ----
if (FALSE) {
prad_legacy_root <- getwd()
library(dplyr)
library(survival)
library(tibble)
g <- file.path(prad_legacy_root, 'data/legacy_project/GEO/')
g5out <- file.path(prad_legacy_root, 'data/legacy_project/GEO/output/')
setwd(g)
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/>1express_normal_prad.RData'))
count<-as.data.frame(read.csv("GSE183019_processed_counts.txt",header=F,sep="\t"))[,-2]
normal <- count
colnames(normal) <- normal[1,]
normal <- normal[,-1]
normal <- normal[-1,]
normal <- as.data.frame(apply(normal,2,as.numeric))
rownames(normal) <- count[c(2:dim(count)[1]),1]
normal0 <- normal[rownames(normal_exp_tpm),colnames(normal_exp_tpm)]
prad0 <- normal[rownames(prad_exp_tpm),colnames(prad_exp_tpm)]
save(prad0,file=file.path(prad_legacy_root, 'data/legacy_project/WGCNA/rdata/0 GEO original.RData'))
library('DESeq2')
library('stringr')
cn <- round(as.matrix(normal0))
cp <- round(as.matrix(prad0))
# which(is.na(cn),arr.ind = T)
dimnames(cn) <-list(rownames(cn),t(as.list(str_c('n_', colnames(cn), sep=""))))# add colname'n_','p_'
cn <- as.matrix(cn)
dimnames(cp) <-list(rownames(cp),t(as.list(str_c('p_', colnames(cp), sep=""))))
cp <- as.matrix(cp)
data <- as.matrix(cbind(cp,cn))
group_type<-as.matrix(rbind(matrix("prad",dim(cp)[2],1),
                            matrix("normal",dim(cn)[2],1)))
colnames(group_type)<-"type"
count_data<-DESeqDataSetFromMatrix(countData=data,colData=group_type,design=~type)
deSeqAnalysis<-DESeq(count_data) 
res = results(deSeqAnalysis, 
              contrast=c("type","prad","normal"), 
              lfcThreshold = 0,
              altHypothesis="greaterAbs",
              alpha = 0.05) 
res1 <- as.data.frame(res)
res0_05 <- res1[which(res1 $ padj < 0.05),]
res0_05 <- mutate(res0_05, up_down = if_else(log2FoldChange < 0 , 'down', 'up')) 
table(res0_05$up_down)
up_dif_exp_05 <- res0_05[which(res0_05$up_down=='up'),]
down_dif_exp_05 <- res0_05[which(res0_05$up_down=='down'),]
table(res0_05$up_down) #1846 2396
save(up_dif_exp_05,down_dif_exp_05,res1,res0_05, 
     file=file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/dif_exp_gene_up_down05.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/dif_exp_gene_up_down05.RData'))
}
# Archived shared signature genes ----
if (FALSE) {
prad_legacy_root <- getwd()
# Shared differentially expressed genes ----
setwd(file.path(prad_legacy_root, 'data/legacy_project/GEO'))
# background gene
load(file.path(prad_legacy_root, 'data/legacy_project/part2 network/bg.RData'))
# two overlap gene list
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/dif_exp_gene_up_down05.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/part1/output/2 DEG/dif_exp_gene_up_down001.RData'))
res01 <- filter(res0_05 , padj < 0.01)
res01up <- filter(res01 , up_down == 'up')
res01down <- filter(res01 , up_down == 'down')
gene_list_1 <- rownames(res01)
gene_list_2 <- rownames(res0_01)
list1_up <- rownames(res01up)
list1_down <- rownames(res01down)
list2_up <- rownames(up_dif_exp_01)
list2_down <- rownames(down_dif_exp_01)
#  function
two_gene_list_overlap<-function(gene_list_1,gene_list_2,list1_up,list1_down,list2_up,list2_down,bg_num){
  
  over_all<-as.matrix(intersect(gene_list_1,gene_list_2)) #jiao
  over_up<-as.matrix(intersect(list1_up,list2_up))
  over_down<-as.matrix(intersect(list1_down,list2_down))
  
  gene_list_1 <- as.matrix(gene_list_1)
  gene_list_2 <- as.matrix(gene_list_2)
  
  union_list<-as.matrix(unique(rbind(gene_list_1,gene_list_2)))
  union_list_num<-dim(union_list)[1]
  
  # p_binom=1-pbinom(dim(over_up)[1]+dim(over_down)[1]-1,dim(over_all)[1],0.5)
  
  p_hyper=1-phyper(dim(over_up)[1]+dim(over_down)[1]-1,length(gene_list_1),
                   bg_num-length(gene_list_1),length(gene_list_2))
  
  over_result <- t(as.matrix(c(length(gene_list_1),length(gene_list_2),dim(over_all)[1],
                               (dim(over_up)[1]+dim(over_down)[1])/dim(over_all)[1],
                               dim(over_up)[1],dim(over_down)[1],p_hyper)))
  return(over_result)
}
res_overlap_deg1_deg <- two_gene_list_overlap(gene_list_1,gene_list_2,list1_up,list1_down,list2_up,list2_down,bg_num)
colnames(res_overlap_deg1_deg)<-c("GEO18 deg","tcga deg","overlap","consistent","co_up","co_down","p_hyper")
#  save
save(res_overlap_deg1_deg,
     file = file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/0.01overlap_deg1_degtcga.RData'))
write.table(res_overlap_deg1_deg,
            file=file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/0.01overlap_deg1_degtcga.txt'),sep="\t",row.names=F,col.names=T,quote=F)
save(over_up,over_down,
     file = file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/0.01overlap_18tcga_up_down.RData'))
}
# Archived signature gene enrichment ----
if (FALSE) {
prad_legacy_root <- getwd()
setwd(file.path(prad_legacy_root, 'data/legacy_project/part3'))
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library("BiocManager")
load(file.path(prad_legacy_root, 'data/legacy_project/GEO/183output/0.01overlap_18tcga_up_down.RData'))
load(file.path(prad_legacy_root, 'data/legacy_project/part2 network/bg.RData'))
#GO
OrgDb <- 'org.Hs.eg.db'
go_up <- enrichGO(gene=over_up,
                  OrgDb='org.Hs.eg.db',
                  universe = bg,
                  keyType = "SYMBOL", 
                  ont="BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff=0.05,
                  pvalueCutoff=0.05,
)
r_go_up <- go_up@result
r_go_up <- filter(r_go_up, pvalue < 0.05)
go_down <- enrichGO(gene=over_down,
                  OrgDb='org.Hs.eg.db',
                  universe = bg,
                  keyType = "SYMBOL", 
                  ont="BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff=0.05,
                  pvalueCutoff=0.05,
)
r_go_down <- go_down@result
r_go_down <- filter(r_go_down, pvalue < 0.05)
write.table(r_go_up, file='result/8 DEG GO up.txt',sep="\t",row.names=F,col.names=T,quote=F)
write.table(r_go_down, file='result/8 DEG GO down.txt',sep="\t",row.names=F,col.names=T,quote=F)
# write.table(r_go_up,
#             file='result/8 DEG GO up.txt',sep="\t",row.names=F,col.names=T,quote=F)
# write.table(r_go_down,
#             file='result/8 DEG GO down.txt',sep="\t",row.names=F,col.names=T,quote=F)
r_go_up1 <- as.data.frame(r_go_up)
r_go_up1 <- arrange(r_go_up1, pvalue)
r_go_up1 <- r_go_up1[1:20,]
r_go_up1$pvaluelog <- abs(log10(as.numeric(r_go_up1$pvalue)))
#dotplot(go_down, showCategory = 20, title='down-regulated GO term')
r_go_down1 <- as.data.frame(r_go_down)
r_go_down1 <- arrange(r_go_down1, pvalue)
r_go_down1 <- r_go_down1[1:20,]
r_go_down1$pvaluelog <- abs(log10(as.numeric(r_go_down1$pvalue)))
up <- r_go_up1[,c(2,10)]
a <- c(rep('DEG up-regulated',20))
up <- cbind(up,a)
colnames(up) <- c('a','b','c')
  
down <- r_go_down1[,c(2,10)]
b <- c(rep('DEG down-regulated',20))
down <- cbind(down,b)
colnames(down) <- c('a','b','c')
  
bind <- rbind(down,up)
c <- c(1:40)
rownames(bind) <- c
GO_term_order <- factor(as.integer(rownames(bind)),labels=rev(bind$a))
png("result/sft/degGO2.png",units="in", width=8, height=10,res=900)
ggplot(data=bind, 
       aes(x=GO_term_order,y=rev(b), fill=c)) + 
  geom_bar(stat="identity", width=0.7) + 
  scale_fill_manual(values=c(  "#f68060","#89C6ED"))+
  
  coord_flip() +  
  guides(fill=FALSE)+
  ylab("-log10(pvalue)") + 
  xlab("") + 
  theme_bw()+
  theme(axis.text.y = element_text( size = 12),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off() 
bind <- arrange(bind, desc(c),desc(b) )
GO_term_order <- factor(as.integer(rownames(bind)),labels=rev(bind$a))
png("result/sft/degGO3.png",units="in", width=8, height=10,res=900)
ggplot(data=bind, 
       aes(x=GO_term_order,y=rev(b), fill=c)) + 
  geom_bar(stat="identity", width=0.7) + 
  scale_fill_manual(values=c(  "#f68060","#89C6ED"))+
  
  coord_flip() +  
  guides(fill=FALSE)+
  ylab("-log10(pvalue)") + 
  xlab("") + 
  theme_bw()+
  theme(axis.text.y = element_text( size = 12),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off() 
png("result/sft/deg_up.png",units="in", width=9, height=6,res=600)
r_go_up1 %>%
  mutate(Description = fct_reorder(Description, pvaluelog)) %>%
  ggplot(aes(x = Description, y = pvaluelog)) +
  geom_bar(stat="identity", fill="#f68060",width = 0.7) +
  coord_flip() +
  
  theme_classic()+
  theme(axis.text.y = element_text( size = 14))+
  labs(x = '',y = "-log10(FDR)") #,title = 'sig pathway of overlap gene (cox/km up-regulated term)'
dev.off() 
png("result/sft/deg_down.png",units="in", width=5, height=6,res=600)
r_go_down1 %>%
  mutate(Description = fct_reorder(Description, pvaluelog)) %>%
  ggplot(aes(x = Description, y = pvaluelog)) +
  geom_bar(stat="identity", fill="#89C6ED",width = 0.7) +
  coord_flip() +
  theme_classic()+
  theme(axis.text.y = element_text( size = 14))+
  labs(x = '',y = "-log10(FDR)") #,title = 'sig pathway of overlap gene (cox/km up-regulated term)'
dev.off() 
# 
# save(r_go_down,r_go_up,go_up,go_down,file = 'output/2.1 GO/0.01/GO_up_01.RData')
# 
# write.table(as.data.frame(r_go_up), file="output/2.1 GO/0.01/GO_up_01.txt",sep = "\t",row.names = F,col.names = T)
# write.table(as.data.frame(r_go_down), file="output/2.1 GO/0.01/GO_down_01.txt",sep = "\t",row.names = F,col.names = T)
}
