# PRAD drug repositioning analysis

This repository accompanies the manuscript **Discovery of Therapeutic Targets and Drugs for Prostate Adenocarcinoma by Drug Repositioning Approach**. The study uses TCGA and GEO transcriptomic data to identify candidate therapeutic targets in prostate adenocarcinoma and compares PC3 perturbation signatures from the Connectivity Map (CMap) to prioritize drugs for repositioning.

## Code files

- [01_cohorts_deg_go.R](analysis/01_cohorts_deg_go.R) — Differential expression and Gene Ontology enrichment analyses for the TCGA and GEO cohorts.
- [02_wgcna.R](analysis/02_wgcna.R) — Overlap analysis of saved WGCNA modules.
- [03_hpa.R](analysis/03_hpa.R) — Candidate gene expression plots using Human Protein Atlas data.
- [04_cmap.R](analysis/04_cmap.R) — PC3 perturbation signature comparison and compound ranking.
- [05_supplementary_tables.R](analysis/05_supplementary_tables.R) — Exports supplementary tables.
