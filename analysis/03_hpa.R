# HPA expression plots ----
data_dir <- "data"
plots <- "plot"
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
themelegendnone <- ggplot2::theme(legend.position = "none")
themelegendtop <- ggplot2::theme(legend.position = "top")
themelegendright <- ggplot2::theme(legend.position = "right")
themetext90 <- ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
themexnone <- ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank())
theme2_noback <- ggplot2::theme(panel.grid = ggplot2::element_blank())
library(dplyr)
# Normal-tissue expression ----
rna_tissue_all <- read.delim(
  file.path(data_dir, "rna_tissue_consensus.tsv"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
rna_tissue <- rna_tissue_all %>% dplyr::filter(Gene.name %in%c('KLK2','KLK3','KLK4','HOXB13','MESP1','LRRC26','SPON2'))
library(ggplot2)
genes_to_plot <- c(
  "KLK2","KLK3","KLK4",
  "HOXB13","MESP1","LRRC26","SPON2"
)
rna_tissue$Tissue <- factor(rna_tissue$Tissue,levels = unique(rna_tissue$Tissue))
tissue_barplots <- list()
my_colors_51 <- c(
  "#1F3A5F", "#234A73", "#275B86", "#2B6C9A", "#2F7DAE",
  "#338EC1", "#3A9ED1", "#4AAAD7", "#5AB5DC", "#6AC0E1",
  "#7ACAE5", "#8AD3E9", "#9ADCEC", "#A9E3E6", "#B7E8DD",
  "#C4ECD2", "#D0F0C6", "#DBF3BA", "#E5F6AE", "#EFF8A2",
  "#F6F8A6", "#FAF4A8", "#FCEFA9", "#FDE9AA", "#FDE3AB",
  "#FDDCAB", "#FCD4A9", "#FBCDA7", "#FAC5A4", "#F8BDA1",
  "#F6B49D", "#F3AB99", "#B5252E", "#EC9991", "#E7908D",
  "#E28789", "#DC7E85", "#D67581", "#CF6C7D", "#C86479",
  "#C15C75", "#BA5471", "#B34C6D", "#AC4569", "#A43E65",
  "#9C3761", "#94305D", "#8C2A59", "#842455", "#7C1E51",
  "#74194D"
)
for (g in genes_to_plot) {
  df_g <- rna_tissue[rna_tissue$Gene.name == g, ]
  tissue_barplots[[g]] <- ggplot(df_g,aes(x = Tissue, y = nTPM,fill=Tissue)) +
    geom_bar(stat = "identity") +labs(title = g,x = NULL,y = "nTPM") +
    scale_fill_manual(values = my_colors_51, name='')+
    theme_bw() + theme2_noback+themetext90+themelegendnone
}
library(dplyr)
library(ggplot2)
library(dplyr)
library(stringr)
rna_tissue <- rna_tissue %>%
  mutate(Tissue = str_to_sentence(Tissue))
tissue_expression <- rna_tissue %>% mutate(log_nTPM = log10(nTPM + 1))
tissue_expression$Gene.name <- factor(tissue_expression$Gene.name,levels = genes_to_plot)
tissue_plot <- ggplot(tissue_expression, aes(x = Tissue, y = Gene.name)) +
    geom_tile(fill = "white",color = "grey90",linewidth = 0.3) +
  geom_point(aes(size = log_nTPM, fill = log_nTPM),shape  = 22,color  = "grey",
    stroke = 0.3) +
  scale_fill_gradientn(colors = c("#324856","#4D8AC6FF", "#F7F7F7", "#E67932FF", "#B2182B"),name = "log10(nTPM + 1)") +
  scale_size(range = c(0, 8),guide = "none") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(size=13,angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  )+themelegendtop+theme(legend.title = element_text(size = 13),
                         legend.text = element_text(size = 12))
ggsave(file.path(plots, "hpa_normal_tissue_expression.pdf"), plot = tissue_plot, width = 32, height = 11, units = "cm")
# Cancer expression ----
gene_id_map <- rna_tissue %>%
  distinct(Gene.name, Gene)
library(data.table)
gene_set <- unique(gene_id_map$Gene)
supplement_genes <- c("SPDEF", "DCXR")
supplement_ids <- rna_tissue_all$Gene[rna_tissue_all$Gene.name %in% supplement_genes]
all_gene_ids <- union(gene_set, supplement_ids)
read_cancer_expression <- function(path, gene_ids) {
  input <- file(path, open = "rt")
  on.exit(close(input))
  header <- readLines(input, n = 1L)
  gene_pattern <- paste0("^(?:", paste(gene_ids, collapse = "|"), ")\t")
  matching_lines <- list()
  repeat {
    lines <- readLines(input, n = 250000L)
    if (!length(lines)) break
    selected <- grep(gene_pattern, lines, value = TRUE, perl = TRUE)
    if (length(selected)) matching_lines[[length(matching_lines) + 1L]] <- selected
  }
  fread(text = paste(c(header, unlist(matching_lines)), collapse = "\n"),
        sep = "\t", select = c("Gene", "Cancer", "pTPM"))
}
rna_cancer_selected <- read_cancer_expression(
  file.path(data_dir, "rna_cancer_sample.tsv"), all_gene_ids
)
rna_cancer <- as.data.frame(rna_cancer_selected[Gene %in% gene_set])
rna_cancer_avg <- rna_cancer %>%group_by(Cancer, Gene) %>%
  summarise(mean_pTPM = mean(pTPM, na.rm = TRUE), .groups = "drop")
rna_cancer_avg <- left_join(rna_cancer_avg,gene_id_map,by='Gene')
mean_cancer_expression <- rna_cancer_avg %>%
  mutate(log_mean_pTPM = log10(mean_pTPM + 1))
tcga_mean_expression <- mean_cancer_expression %>%
  filter(str_detect(Cancer, "\\(TCGA\\)")) %>%
  mutate(Cancer = str_remove(Cancer, " \\(TCGA\\)"))
cancer_plot <- ggplot(tcga_mean_expression, aes(x = Cancer, y = Gene.name)) +
    geom_tile(fill = "white",color = "grey90",linewidth = 0.3) +
  geom_point(aes(size = log_mean_pTPM, fill = log_mean_pTPM),shape  = 22,color  = "grey",
             stroke = 0.3) +
  scale_fill_gradientn(colors = c("#F7F7F7", "#9ADCEC", "#3B8BC4", "#B2182B"),name = "log10(nTPM + 1)"  ) +
  scale_size(range = c(0.5, 12),guide = "none") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x  = element_text(size=13,angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(size=12),axis.title.x = element_blank(),axis.title.y = element_blank(),
    panel.grid   = element_blank())+
  theme(axis.text.x = element_text( size = 13, angle = 30,vjust = 1,hjust= 1 ))+theme(legend.title = element_text(size = 13),
                                                                                      legend.text = element_text(size = 12))
ggsave(file.path(plots, "hpa_mean_cancer_expression.pdf"), plot = cancer_plot, width = 32, height = 13, units = "cm")
genes_to_plot <- c(
  "KLK2","KLK3","KLK4",
  "HOXB13","MESP1","LRRC26","SPON2"
)
cancer_expression_with_symbols <- left_join(rna_cancer,gene_id_map,by='Gene')
sample_expression <- cancer_expression_with_symbols %>%
  mutate(log_pTPM = log10(pTPM + 1))
tcga_sample_expression <- sample_expression %>%
  filter(str_detect(Cancer, "\\(TCGA\\)")) %>%
  mutate(Cancer = str_remove(Cancer, " \\(TCGA\\)"))
my_colors_21 <- c(
  "#1F3A5F", "#275B86", "#2F7DAE", "#3A9ED1", "#5AB5DC",
  "#7ACAE5", "#9ADCEC", "#B7E8DD", "#D0F0C6", "#E5F6AE",
  "#F6F8A6", "#FDE3AB", "#FAC5A4", "#F3AB99", "#B5252E",
  "#E28789", "#D67581", "#C86479", "#BA5471", "#AC4569",
  "#9C3761"
)
my_colors_21_light <- c(
  "#60758C", "#6C83A0", "#7892B3", "#85A1C6", "#93B0D8",
  "#A1BFE2", "#B0CDE8", "#C0DAED", "#D0E7F1", "#E0F2F4",
  "#EEF5CC", "#F6EACC", "#F4D3C7", "#EFC1BB", "#B5252E",
  "#E7B8B3", "#D8A5A0", "#C78F8B", "#B47878", "#A26366",
  "#8F4F54"
)
cancer_boxplots <- list()
for(i in genes_to_plot){
  gene_expression <- tcga_sample_expression %>% dplyr::filter(Gene.name == i)
  cancer_boxplots[[i]] <- ggplot(gene_expression, aes(x = Cancer, y = log_pTPM)) +
      geom_point(aes(fill=Cancer,color=Cancer),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),alpha=1,
                 shape = 21, size = 0.1) +
    geom_boxplot(
      position = position_dodge(width = 0.8),width = 0.5,size = 0.4,colour = "black",fill = NA,
      outlier.shape = NA, alpha = 1)+
      scale_fill_manual(values =  my_colors_21_light) +
    scale_color_manual(values =  my_colors_21_light) +
      labs(title=i,x = NULL,y = "log10(pTPM + 1)",fill =NULL) +
      theme2_noback+themetext90+themelegendnone+theme(axis.title.x = element_text(size=15,face='bold',colour="black"),
                                                                         axis.text.x = element_text( size = 15 ,angle = 50,vjust = 1,hjust= 1 ),
                                                                         axis.title.y = element_text(size=15,face='bold',colour="black"),
                                                                         axis.text.y = element_text( size = 15,colour="black"),
                                                      plot.title = element_text( size = 13)
      )
}
library(patchwork)
hpa_panels <- wrap_plots(
  (cancer_boxplots[[1]]+themexnone),( cancer_boxplots[[2]]+themexnone), (cancer_boxplots[[3]]+themexnone), (cancer_boxplots[[4]]+themexnone),
 ( cancer_boxplots[[5]]+themexnone),( cancer_boxplots[[6]]+themexnone), (cancer_boxplots[[7]]+themexnone),
  ncol = 4
)
ggsave(file.path(plots, "hpa_cancer_expression_boxplots.pdf"), plot = hpa_panels, width  = 32,height = 13,units  = "cm")
legend_plot <- cancer_boxplots[[1]] + geom_point(aes(fill=Cancer,color=Cancer),position = position_jitterdodge(jitter.width = 2, dodge.width = 0.8),alpha=1,
                    size = 2)+themelegendright
ggsave(file.path(plots, "hpa_cancer_expression_legend.pdf"), plot = legend_plot, width  = 32,height = 13,units  = "cm")
# SPDEF and DCXR expression ----
rna_tissue <- rna_tissue_all %>% dplyr::filter(Gene.name %in% supplement_genes)
library(ggplot2)
genes_to_plot <- c('SPDEF','DCXR')
library(dplyr)
library(ggplot2)
library(dplyr)
library(stringr)
rna_tissue <- rna_tissue %>%
  mutate(Tissue = str_to_sentence(Tissue))
supplementary_tissue_expression <- rna_tissue %>% mutate(log_nTPM = log10(nTPM + 1))
supplementary_tissue_expression$Gene.name <- factor(supplementary_tissue_expression$Gene.name,levels = genes_to_plot)
supplement_tissue_plot <- ggplot(supplementary_tissue_expression, aes(x = Tissue, y = Gene.name)) +
  geom_tile(fill = "white",color = "grey90",linewidth = 0.3) +
  geom_point(aes(size = log_nTPM, fill = log_nTPM),shape  = 22,color  = "grey",
             stroke = 0.3) +
  scale_fill_gradientn(colors = c("#324856","#4D8AC6FF", "#F7F7F7", "#E67932FF", "#B2182B"),name = "log10(nTPM + 1)") +
  scale_size(range = c(0, 8),guide = "none") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(size=13,angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  )+themelegendtop+theme(legend.title = element_text(size = 13),
                         legend.text = element_text(size = 12))
ggsave(file.path(plots, "hpa_spdef_dcxr_normal_tissue_expression.pdf"), plot = supplement_tissue_plot, width = 32, height = 9, units = "cm")
gene_id_map <- rna_tissue %>%
  distinct(Gene.name, Gene)
rna_cancer <- as.data.frame(rna_cancer_selected[Gene %in% supplement_ids])
genes_to_plot <- c('SPDEF','DCXR')
cancer_expression_with_symbols <- left_join(rna_cancer,gene_id_map,by='Gene')
supplementary_sample_expression <- cancer_expression_with_symbols %>%
  mutate(log_pTPM = log10(pTPM + 1))
tcga_supplementary_expression <- supplementary_sample_expression %>%
  filter(str_detect(Cancer, "\\(TCGA\\)")) %>%
  mutate(Cancer = str_remove(Cancer, " \\(TCGA\\)"))
my_colors_21 <- c(
  "#1F3A5F", "#275B86", "#2F7DAE", "#3A9ED1", "#5AB5DC",  "#7ACAE5", "#9ADCEC", "#B7E8DD", "#D0F0C6", "#E5F6AE",  "#F6F8A6", "#FDE3AB", "#FAC5A4", "#F3AB99", "#B5252E",
  "#E28789", "#D67581", "#C86479", "#BA5471", "#AC4569",  "#9C3761")
my_colors_21_light <- c(
  "#60758C", "#6C83A0", "#7892B3", "#85A1C6", "#93B0D8",  "#A1BFE2", "#B0CDE8", "#C0DAED", "#D0E7F1", "#E0F2F4",
  "#EEF5CC", "#F6EACC", "#F4D3C7", "#EFC1BB", "#B5252E",  "#E7B8B3", "#D8A5A0", "#C78F8B", "#B47878", "#A26366",  "#8F4F54")
supplementary_boxplots <- list()
for(i in genes_to_plot){
  gene_expression <- tcga_supplementary_expression %>% dplyr::filter(Gene.name == i)
  supplementary_boxplots[[i]] <- ggplot(gene_expression, aes(x = Cancer, y = log_pTPM)) +
    geom_point(aes(fill=Cancer,color=Cancer),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),alpha=1,shape = 21, size = 0.1) +
    geom_boxplot(position = position_dodge(width = 0.8),width = 0.5,size = 0.4,colour = "black",fill = NA,outlier.shape = NA, alpha = 1)+
    scale_fill_manual(values =  my_colors_21_light) +  scale_color_manual(values =  my_colors_21_light) +
    labs(title=i,x = NULL,y = "log10(pTPM + 1)",fill =NULL) +
    theme2_noback+themetext90+themelegendnone+theme(axis.title.x = element_text(size=13,face='bold',colour="black"),
                                                    axis.text.x = element_text( size = 13 ,angle = 90,vjust = 0.5,hjust= 1 ),
                                                    axis.title.y = element_text(size=13,face='bold',colour="black"),
                                                    axis.text.y = element_text( size = 13,colour="black"),
                                                    plot.title = element_text( size = 13)
    )
}
library(patchwork)
hpa_panels <- wrap_plots(
  (supplementary_boxplots[[1]]),( supplementary_boxplots[[2]]),
  ncol = 2
)
ggsave(file.path(plots, "hpa_spdef_dcxr_cancer_expression_boxplots.pdf"), plot = hpa_panels, width  = 28,height = 23,units  = "cm")
legend_plot <- supplementary_boxplots[[1]] + geom_point(aes(fill=Cancer,color=Cancer),position = position_jitterdodge(jitter.width = 2, dodge.width = 0.8),alpha=1,
                   size = 2)+themelegendright
ggsave(file.path(plots, "hpa_spdef_dcxr_cancer_expression_legend.pdf"), plot = legend_plot, width  = 32,height = 13,units  = "cm")
message("Plots saved to ", normalizePath(plots))
