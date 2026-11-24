# CMap compound ranking ----
data_dir <- "data"
results <- file.path("result", "cmap")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
read_gctx <- function(file, field) rhdf5::h5read(file, field)
read_column_metadata <- function(file, fields) {
  setNames(as.data.frame(lapply(fields, function(field)
    read_gctx(file, paste0("/0/META/COL/", field)))), fields)
}
drug_file <- file.path(data_dir, "PC3_cp.gctx")
knockout_file <- file.path(data_dir, "PC3_ko.gctx")
genetic_files <- c(file.path(data_dir, "PC3_sh.gctx"), knockout_file[file.exists(knockout_file)])
targets <- c("KLK3","MESP1","SPDEF","HOXB13","DCXR")
drug_metadata <- read_column_metadata(drug_file, c("id","cell_iname","pert_type","pert_id","cmap_name","pert_idose","pert_itime"))
drug_gene_ids <- read_gctx(drug_file,"/0/META/ROW/id")
matrix_path <- "/0/DATA/0/matrix"
genetic_inputs <- lapply(genetic_files, function(file) {
  meta <- read_column_metadata(file, c("id","cell_iname","pert_type","cmap_name"))
  keep <- which(meta$cell_iname == "PC3" & meta$pert_type %in% c("trt_sh","trt_xpr") &
                  meta$cmap_name %in% targets)
  if (!length(keep)) return(NULL)
  genes <- read_gctx(file,"/0/META/ROW/id")
  stopifnot(!anyDuplicated(genes), all(drug_gene_ids %in% genes))
  values <- rhdf5::h5read(file,matrix_path,index=list(NULL,keep))
  list(meta=meta[keep,], matrix=values[match(drug_gene_ids, genes),,drop=FALSE])
})
meta <- do.call(rbind,lapply(genetic_inputs,`[[`,"meta"))
genetic_signature_matrix <- do.call(cbind,lapply(genetic_inputs,`[[`,"matrix"))
stopifnot(nrow(meta)>0, !anyDuplicated(meta$id))
missing <- setdiff(targets, meta$cmap_name)
write.csv(meta,file.path(results,"cmap_genetic_signatures.csv"),row.names=FALSE)
writeLines(missing,file.path(results,"cmap_missing_targets.txt"))
batch_size <- 512L
drug_signature_rows <- which(drug_metadata$cell_iname=="PC3" & drug_metadata$pert_type=="trt_cp" & drug_metadata$pert_id!="" & drug_metadata$pert_id!="-666")
scores <- matrix(NA_real_,nrow(drug_metadata),nrow(meta))
for (start in seq(1,length(drug_signature_rows),by=batch_size)) {
  columns <- drug_signature_rows[start:min(start+batch_size-1L,length(drug_signature_rows))]
  values <- rhdf5::h5read(drug_file,matrix_path,index=list(NULL,columns))
  scores[columns,] <- cor(values,genetic_signature_matrix,method="spearman")
}
for (target in intersect(targets,meta$cmap_name)) {
  signatures <- which(meta$cmap_name==target)
  best_per_drug <- lapply(signatures,function(i) {
    ordered <- drug_signature_rows[order(-scores[drug_signature_rows,i],drug_metadata$id[drug_signature_rows],na.last=NA)]
    rows <- ordered[!duplicated(drug_metadata$pert_id[ordered])]
    data.frame(drug_metadata[rows,],rho=scores[rows,i])
  })
  median_score <- vapply(best_per_drug,function(d) median(d$rho),numeric(1))
  winner <- which.max(median_score)
  output <- best_per_drug[[winner]]
  output$rank <- seq_len(nrow(output))
  output$target <- target
  output$genetic_signature <- meta$id[signatures[winner]]
  write.csv(output,file.path(results,paste0("cmap_", tolower(target), "_compound_ranking.csv")),row.names=FALSE)
  write.csv(data.frame(meta[signatures,], median_score),file.path(results,paste0("cmap_", tolower(target), "_signature_scores.csv")),row.names=FALSE)
}
message("Results saved to ", results)
