# ==================================================================================
# RNA-seq DGE Analysis: E6.5 vs E8.5 (FPKM-normalized data)
# ==================================================================================

# Load packages
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(limma)
library(tidyverse)
library(pheatmap)
library(enrichR)
library(gridExtra)

# ==================================================================================
# SECTION 1: Load Gene Annotation
# ==================================================================================
txdb <- makeTxDbFromGFF("data/Mus_musculus.GRCm39.115.gtf")
genes <- genes(txdb)
gtf <- import("data/Mus_musculus.GRCm39.115.gtf")
gene_info <- gtf[gtf$type == "gene"]
gene_map <- data.frame(
  ensembl_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  stringsAsFactors = FALSE
)
gene_map <- gene_map[!duplicated(gene_map$ensembl_id), ]
rownames(gene_map) <- gene_map$ensembl_id

# ==================================================================================
# SECTION 2: Import BigWig Files and Extract FPKM
# ==================================================================================
convert_to_symbols <- function(ensembl_ids, gene_map) {
  ensembl_ids_clean <- sub("\\.\\d+$", "", ensembl_ids)
  symbols <- gene_map$gene_name[match(ensembl_ids_clean, gene_map$ensembl_id)]
  symbols[is.na(symbols) | symbols == ""] <- ensembl_ids[is.na(symbols) | symbols == ""]
  if (any(duplicated(symbols))) {
    dup_symbols <- unique(symbols[duplicated(symbols)])
    for (sym in dup_symbols) {
      idx <- which(symbols == sym)
      symbols[idx] <- paste0(sym, "_", ensembl_ids[idx])
    }
  }
  numeric_symbols <- grepl("^[0-9]+$", symbols)
  symbols[numeric_symbols] <- paste0("GENE_", symbols[numeric_symbols])
  return(symbols)
}

bw_files <- list.files(path = "data/rna_bw", pattern = "\\.bw$", full.names = TRUE)
sample_names <- gsub("\\.bw$", "", basename(bw_files))
bw_test <- import(bw_files[1], format = "BigWig")
bw_chroms <- as.character(unique(seqnames(bw_test)))
gtf_chroms <- as.character(unique(seqnames(genes)))

if (!any(bw_chroms %in% gtf_chroms)) {
  if (all(grepl("^chr", bw_chroms)) && !any(grepl("^chr", gtf_chroms))) {
    seqlevels(genes) <- paste0("chr", seqlevels(genes))
    seqlevels(genes)[seqlevels(genes) == "chrMT"] <- "chrM"
  } else if (!any(grepl("^chr", bw_chroms)) && all(grepl("^chr", gtf_chroms))) {
    seqlevels(genes) <- gsub("^chr", "", seqlevels(genes))
    seqlevels(genes)[seqlevels(genes) == "M"] <- "MT"
  }
}

fpkm_matrix <- matrix(0, nrow = length(genes), ncol = length(bw_files))
rownames(fpkm_matrix) <- names(genes)
colnames(fpkm_matrix) <- sample_names

for (i in seq_along(bw_files)) {
  message("Processing sample ", i, " of ", length(bw_files), ": ", sample_names[i])
  bw_cov <- import(bw_files[i], format = "BigWig", as = "RleList")
  common_chroms <- intersect(seqlevels(genes), names(bw_cov))
  genes_subset <- keepSeqlevels(genes, common_chroms, pruning.mode = "coarse")
  seqlevels(genes_subset) <- common_chroms
  bw_cov <- bw_cov[common_chroms]
  genes_by_chr <- split(genes_subset, seqnames(genes_subset))
  gene_views <- Views(bw_cov, ranges(genes_by_chr))
  gene_means_list <- lapply(gene_views, viewMeans)
  gene_means <- unlist(gene_means_list)
  all_gene_names <- unlist(lapply(genes_by_chr, names))
  names(gene_means) <- all_gene_names
  fpkm_matrix[all_gene_names, i] <- gene_means
}

rownames(fpkm_matrix) <- sub("\\.\\d+$", "", rownames(fpkm_matrix))
rownames(fpkm_matrix) <- convert_to_symbols(rownames(fpkm_matrix), gene_map)
write.table(fpkm_matrix, "fpkm_matrix_all_samples.tsv", sep = "\t", quote = FALSE, col.names = NA)

# ==================================================================================
# SECTION 3: Load Metadata
# ==================================================================================
coldata <- read.csv("metadata.csv", stringsAsFactors = FALSE)
coldata$stage <- factor(coldata$stage, levels = c("E6.5", "E8.5"))
coldata$condition <- factor(coldata$condition, levels = c("WT", "KO"))
coldata$lineage <- factor(coldata$lineage)

# ==================================================================================
# SECTION 4: Filter Low Expression Genes
# ==================================================================================
keep <- rowMeans(fpkm_matrix) >= 1
fpkm_filtered <- fpkm_matrix[keep, ]
message(paste("Retained", nrow(fpkm_filtered), "genes after filtering"))

# Log2-transform FPKM with small offset
log2_fpkm <- log2(fpkm_filtered + 0.1)

# ==================================================================================
# SECTION 5: Differential Expression Analysis with limma
# ==================================================================================
design <- model.matrix(~ stage, data = coldata)
fit <- lmFit(log2_fpkm, design)

# eBayes with trend=TRUE for mean-variance modeling of log2(FPKM)
fit <- eBayes(fit, trend = TRUE)

res_E85_vs_E65 <- topTable(fit, coef = "stageE8.5", number = Inf, adjust.method = "BH", sort.by = "none")
colnames(res_E85_vs_E65)[colnames(res_E85_vs_E65) == "logFC"] <- "log2FoldChange"
colnames(res_E85_vs_E65)[colnames(res_E85_vs_E65) == "adj.P.Val"] <- "padj"
colnames(res_E85_vs_E65)[colnames(res_E85_vs_E65) == "P.Value"] <- "pvalue"

# ==================================================================================
# SECTION 6: Extract Significant Genes
# ==================================================================================
get_sig_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  sig <- res[!is.na(res$padj) & res$padj < padj_cutoff & abs(res$log2FoldChange) > lfc_cutoff, ]
  return(sig)
}

sig_E85_vs_E65 <- get_sig_genes(res_E85_vs_E65)
sig_E85_vs_E65_UP <- sig_E85_vs_E65[sig_E85_vs_E65$log2FoldChange > 0, ]
sig_E85_vs_E65_DOWN <- sig_E85_vs_E65[sig_E85_vs_E65$log2FoldChange < 0, ]

message(paste("Found", nrow(sig_E85_vs_E65_UP), "upregulated genes"))
message(paste("Found", nrow(sig_E85_vs_E65_DOWN), "downregulated genes"))

# Save gene lists
up_genes <- sig_E85_vs_E65_UP$ID
down_genes <- sig_E85_vs_E65_DOWN$ID
writeLines(if(length(up_genes) > 0) up_genes else character(0), "E85_vs_E65_UP.txt")
writeLines(if(length(down_genes) > 0) down_genes else character(0), "E85_vs_E65_DOWN.txt")
write.csv(res_E85_vs_E65, "E85_vs_E65_all_results.csv")

# ==================================================================================
# SECTION 7: EnrichR Pathway Enrichment
# ==================================================================================
dbs <- listEnrichrDbs()
dbs_selection <- unique(dbs$libraryName)[grep("Mouse", unique(dbs$libraryName))]

plot_top_enrichments <- function(enriched_list, top_n = 20, filename_plot = NULL, filename_csv = NULL, plot_title = NULL) {
  
  df <- bind_rows(lapply(names(enriched_list), function(db) {
    res <- enriched_list[[db]]
    if (!is.null(res) && nrow(res) > 0) { 
      res$Database <- db
      return(res)
    }
  }))
  
  if (nrow(df) == 0) { 
    message("No enriched terms found.")
    return(NULL)
  }
  
  df_top <- df %>% arrange(Adjusted.P.value) %>% slice_head(n = top_n)
  
  if (!is.null(filename_csv)) write.csv(df_top, filename_csv, row.names = FALSE)
  
  df_top$Term_short <- sapply(df_top$Term, function(x) {
    if (nchar(x) > 60) paste0(substr(x, 1, 57), "...") else x
  })
  
  title_text <- ifelse(is.null(plot_title), paste("Top", top_n, "Enriched Mouse Terms"), plot_title)
  
  p <- ggplot(df_top, aes(x = reorder(Term_short, -log10(Adjusted.P.value)), 
                          y = -log10(Adjusted.P.value), fill = Database)) +
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_bw(base_size = 12) +
    labs(title = title_text, x = "", y = "-log10(Adjusted P-value)") + 
    theme(
      legend.position = "left",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 10)
    )
  
  # Save plot with larger size and higher resolution
  if (!is.null(filename_plot)) {
    ggsave(
      filename = filename_plot,
      plot = p,
      width = 16,
      height = 12,
      dpi = 300
    )
  }
  
  return(p)
}

mouse_genes_UP <- if(file.exists("E85_vs_E65_UP.txt")) readLines("E85_vs_E65_UP.txt") else character(0)
mouse_genes_UP <- mouse_genes_UP[nchar(mouse_genes_UP) > 0]
mouse_genes_DOWN <- if(file.exists("E85_vs_E65_DOWN.txt")) readLines("E85_vs_E65_DOWN.txt") else character(0)
mouse_genes_DOWN <- mouse_genes_DOWN[nchar(mouse_genes_DOWN) > 0]

if(length(mouse_genes_UP) > 0) {
  message("Running enrichment analysis for upregulated genes...")
  enriched_UP <- enrichr(mouse_genes_UP, dbs_selection)
  plot_top_enrichments(
    enriched_UP, top_n = 20,
    filename_plot = "enrichr_up_top20_enrichment_large.png",
    filename_csv = "enrichr_up_top20_enrichment.csv",
    plot_title = "Top 20 Enriched Pathways (UP-regulated Genes, E8.5 vs E6.5)"
  )
} else { message("No upregulated genes found - skipping enrichment analysis") }

if(length(mouse_genes_DOWN) > 0) {
  message("Running enrichment analysis for downregulated genes...")
  enriched_DOWN <- enrichr(mouse_genes_DOWN, dbs_selection)
  plot_top_enrichments(
    enriched_DOWN, top_n = 20,
    filename_plot = "enrichr_down_top20_enrichment_large.png",
    filename_csv = "enrichr_down_top20_enrichment.csv",
    plot_title = "Top 20 Enriched Pathways (DOWN-regulated Genes, E8.5 vs E6.5)"
  )
} else { message("No downregulated genes found - skipping enrichment analysis") }
