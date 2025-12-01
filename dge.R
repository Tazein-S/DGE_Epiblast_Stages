# ==================================================================================
# RNA-seq DGE Analysis: E6.5 vs E8.5 
# ==================================================================================
# Load up all the packages we need
library(rtracklayer)        # for reading GTF and bigWig files
library(GenomicFeatures)    # for making TxDb objects
library(GenomicRanges)      # for working with genomic coordinates
library(DESeq2)             # for DGE
library(tidyverse)          # data wrangling and plotting
library(pheatmap)           # for heatmaps
library(enrichR)            # pathway enrichment
library(RDAVIDWebService)   # DAVID annotation tool
library(gridExtra)          # for arranging plots

# ==================================================================================
# SECTION 1: Load Gene Annotation
# ==================================================================================
# First we need to load the GTF file and get all the gene info

# Make a transcript database from the Ensembl GTF
txdb <- makeTxDbFromGFF("data/Mus_musculus.GRCm39.115.gtf")

# Grab the gene coordinates
genes <- genes(txdb)

# Load the full GTF again to get gene names and IDs
gtf <- import("data/Mus_musculus.GRCm39.115.gtf")
gene_info <- gtf[gtf$type == "gene"]

# Make a lookup table: Ensembl ID -> Gene Symbol
gene_map <- data.frame(
  ensembl_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  stringsAsFactors = FALSE
)

# Get rid of duplicates and set up row names
gene_map <- gene_map[!duplicated(gene_map$ensembl_id), ]
rownames(gene_map) <- gene_map$ensembl_id

# ==================================================================================
# SECTION 2: Import BigWig Files and Count Genes
# ==================================================================================
# Now let's read in the bigWig files and figure out how much each gene is expressed

# Find all the bigWig files
bw_files <- list.files(
  path = "data/rna_bw",
  pattern = "\\.bw$",
  full.names = TRUE
)

# Get sample names from the file names
sample_names <- gsub("\\.bw$", "", basename(bw_files))

# --- Fix chromosome naming if needed ---
# Sometimes bigWig files say "chr1" and GTFs say "1", or vice versa
# Let's check and fix that
bw_test <- import(bw_files[1], format = "BigWig")
bw_chroms <- as.character(unique(seqnames(bw_test)))
gtf_chroms <- as.character(unique(seqnames(genes)))

if (!any(bw_chroms %in% gtf_chroms)) {
  # If bigWig has "chr" but GTF doesn't, add "chr" to GTF
  if (all(grepl("^chr", bw_chroms)) && !any(grepl("^chr", gtf_chroms))) {
    seqlevels(genes) <- paste0("chr", seqlevels(genes))
    seqlevels(genes)[seqlevels(genes) == "chrMT"] <- "chrM"
  }
  # If GTF has "chr" but bigWig doesn't, remove "chr" from GTF
  else if (!any(grepl("^chr", bw_chroms)) && all(grepl("^chr", gtf_chroms))) {
    seqlevels(genes) <- gsub("^chr", "", seqlevels(genes))
    seqlevels(genes)[seqlevels(genes) == "M"] <- "MT"
  }
}

# Set up an empty count matrix (rows = genes, columns = samples)
gene_counts <- matrix(0, nrow = length(genes), ncol = length(bw_files))
rownames(gene_counts) <- names(genes)
colnames(gene_counts) <- sample_names

# --- Loop through each bigWig file and count up gene expression ---
# This basically sums up the coverage over each gene
for (i in seq_along(bw_files)) {
  # Load the bigWig data
  bw_gr <- import(bw_files[i], format = "BigWig")
  
  # Find where genes and bigWig regions overlap
  overlaps <- findOverlaps(genes, bw_gr)
  
  if (length(overlaps) > 0) {
    # Get the genes and bigWig chunks that overlap
    gene_hits <- genes[queryHits(overlaps)]
    bw_hits <- bw_gr[subjectHits(overlaps)]
    
    # How much do they overlap?
    overlap_widths <- width(pintersect(gene_hits, bw_hits))
    
    # Weight the coverage by how much they overlap
    weighted_scores <- bw_hits$score * overlap_widths
    
    # Add up all the weighted scores for each gene
    gene_sums <- tapply(weighted_scores, queryHits(overlaps), sum)
    
    # Store the counts
    gene_counts[as.integer(names(gene_sums)), i] <- gene_sums
  }
}

# ==================================================================================
# SECTION 3: Load Metadata and Filter to Just E6.5 and E8.5
# ==================================================================================
# We only want to compare the last two developmental stages

# Load the sample info
coldata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Keep only E6.5 and E8.5 samples
coldata <- coldata[coldata$stage %in% c("E6.5", "E8.5"), ]

# Make sure the factors are in the right order
coldata$stage <- factor(coldata$stage, levels = c("E6.5", "E8.5"))
coldata$condition <- factor(coldata$condition, levels = c("WT", "KO"))
coldata$lineage <- factor(coldata$lineage)

# Match up the count matrix with the filtered metadata
gene_counts <- gene_counts[, coldata$sampleName]

# ================================================================
# SAVE COUNTS MATRIX
# ================================================================
write.table(
  gene_counts,
  file = "counts_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ==================================================================================
# HELPER FUNCTION: Convert Ensembl IDs to Gene Symbols
# ==================================================================================
# Because "ENSMUSG00000..." is way less readable than actual gene names

convert_to_symbols <- function(ensembl_ids, gene_map) {
  # Look up the gene symbols
  symbols <- gene_map$gene_name[match(ensembl_ids, gene_map$ensembl_id)]
  
  # If there's no symbol, just keep the Ensembl ID
  symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
  
  # If we have duplicate gene names, add numbers to make them unique
  if (any(duplicated(symbols))) {
    dup_symbols <- unique(symbols[duplicated(symbols)])
    for (sym in dup_symbols) {
      idx <- which(symbols == sym)
      symbols[idx] <- paste0(sym, "_", seq_along(idx))
    }
  }
  
  return(symbols)
}

# ==================================================================================
# SECTION 4: DESeq2 Time!
# ==================================================================================
# This is where the actual differential expression analysis happens

# Create the DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(gene_counts),
  colData = coldata,
  design = ~ stage
)

# Filter out genes with super low counts (keeps things cleaner)
# Only keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run the whole DESeq2 pipeline (normalization, dispersion, testing)
dds <- DESeq(dds)

# Transform the data for visualization (variance stabilizing transformation)
vsd <- vst(dds, blind = FALSE)

# Get the results comparing E8.5 vs E6.5
res_E85_vs_E65 <- results(
  dds,
  contrast = c("stage", "E8.5", "E6.5"),
  alpha = 0.05
)

# ==================================================================================
# SECTION 5: Extract Significantly Differentially Expressed Genes
# ==================================================================================
# Filter genes based on adjusted p-value and log2 fold change thresholds

# Function to extract significant genes
get_sig_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  # Filter by adjusted p-value and absolute log2 fold change
  sig <- res[
    !is.na(res$padj) &
      res$padj < padj_cutoff &
      abs(res$log2FoldChange) > lfc_cutoff,
  ]
  
  # Convert to data frame
  sig_df <- as.data.frame(sig)
  
  return(sig_df)
}

# Get all significant genes (still with Ensembl IDs)
sig_E85_vs_E65 <- get_sig_genes(res_E85_vs_E65)

# Convert Ensembl IDs to gene symbols for the significant genes
sig_gene_symbols <- convert_to_symbols(rownames(sig_E85_vs_E65), gene_map)
rownames(sig_E85_vs_E65) <- sig_gene_symbols

# Split into upregulated and downregulated genes
sig_E85_vs_E65_UP <- sig_E85_vs_E65[sig_E85_vs_E65$log2FoldChange > 0, ]
sig_E85_vs_E65_DOWN <- sig_E85_vs_E65[sig_E85_vs_E65$log2FoldChange < 0, ]

# ==================================================================================
# SECTION 6: Save Gene Lists for Downstream Analysis
# ==================================================================================
# Export lists of up- and down-regulated genes to text files

# Extract gene names
up_genes <- rownames(sig_E85_vs_E65_UP)
down_genes <- rownames(sig_E85_vs_E65_DOWN)

# Save to files
if (length(up_genes) > 0) {
  writeLines(up_genes, "E85_vs_E65_UP.txt")
}

if (length(down_genes) > 0) {
  writeLines(down_genes, "E85_vs_E65_DOWN.txt")
}

# ==================================================================================
# SECTION 7: EnrichR Pathway Enrichment Analysis
# ==================================================================================
# Perform functional enrichment analysis using mouse-specific databases

# Get list of available EnrichR databases
dbs <- listEnrichrDbs()

# Select mouse-specific databases
dbs_selection <- unique(dbs$libraryName)[grep("Mouse", unique(dbs$libraryName))]

# --- Function to plot top enrichment results ---
plot_top_enrichments <- function(enriched_list,
                                 top_n = 20,
                                 filename_plot = NULL,
                                 filename_csv = NULL,
                                 plot_title = NULL) {
  # Combine results from all databases
  df <- bind_rows(lapply(names(enriched_list), function(db) {
    res <- enriched_list[[db]]
    if (!is.null(res) && nrow(res) > 0) {
      res$Database <- db
      return(res)
    }
  }))
  
  # Check if any enriched terms found
  if (nrow(df) == 0) {
    message("No enriched terms found.")
    return(NULL)
  }
  
  # Select top N most significant terms
  df_top <- df %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = top_n)
  
  # Save results to CSV if requested
  if (!is.null(filename_csv)) {
    write.csv(df_top, filename_csv, row.names = FALSE)
  }
  
  # Shorten long term names for better visualization
  df_top$Term_short <- sapply(df_top$Term, function(x) {
    if (nchar(x) > 60) {
      paste0(substr(x, 1, 57), "...")
    } else {
      x
    }
  })
  
  # Create plot title
  title_text <- ifelse(
    is.null(plot_title),
    paste("Top", top_n, "Enriched Mouse Terms"),
    plot_title
  )
  
  # Generate bar plot
  p <- ggplot(
    df_top,
    aes(
      x = reorder(Term_short, -log10(Adjusted.P.value)),
      y = -log10(Adjusted.P.value),
      fill = Database
    )
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      title = title_text,
      x = "",
      y = "-log10(Adjusted P-value)"
    ) +
    theme(legend.position = "bottom")
  
  # Save plot if requested
  if (!is.null(filename_plot)) {
    ggsave(filename_plot, p, width = 10, height = 6, dpi = 300)
  }
  
  return(p)
}

# --- Run enrichment analysis ---
# Read gene lists
mouse_genes_UP <- readLines("E85_vs_E65_UP.txt")
mouse_genes_DOWN <- readLines("E85_vs_E65_DOWN.txt")

# Perform enrichment analysis
enriched_UP <- enrichr(mouse_genes_UP, dbs_selection)
enriched_DOWN <- enrichr(mouse_genes_DOWN, dbs_selection)

# --- Plot upregulated gene enrichments ---
plot_top_enrichments(
  enriched_UP,
  top_n = 20,
  filename_plot = "enrichr_up_top20_enrichment.png",
  filename_csv = "enrichr_up_top20_enrichment.csv",
  plot_title = "Top 20 Enriched Pathways (UP-regulated Genes, E8.5 vs E6.5)"
)

# --- Plot downregulated gene enrichments ---
plot_top_enrichments(
  enriched_DOWN,
  top_n = 20,
  filename_plot = "enrichr_down_top20_enrichment.png",
  filename_csv = "enrichr_down_top20_enrichment.csv",
  plot_title = "Top 20 Enriched Pathways (DOWN-regulated Genes, E8.5 vs E6.5)"
)

# ==================================================================================
# SECTION 8: PCA and Volcano Plot Visualization
# ==================================================================================
# Generate plots to visualize sample relationships and differential expression

# --- Function to create PCA plot ---
plot_pca <- function(vsd, coldata_subset, title, filename) {
  # Perform PCA and extract data
  pca_data <- plotPCA(
    vsd,
    intgroup = c("stage", "lineage"),
    returnData = TRUE
  )
  
  # Get variance explained by each PC
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  # Create PCA plot
  p <- ggplot(
    pca_data,
    aes(x = PC1, y = PC2, color = stage, shape = lineage)
  ) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw(base_size = 14) +
    labs(title = title) +
    scale_color_manual(values = c("E6.5" = "#00BA38", "E8.5" = "#619CFF"))
  
  # Save plot
  ggsave(filename, p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

# Generate PCA plot
plot_pca(vsd, coldata, "PCA: E6.5 vs E8.5", "PCA_E65_vs_E85.png")

# --- Function to create volcano plot ---
plot_volcano <- function(res, title, filename) {
  # Convert results to data frame and add gene symbols
  res_df <- as.data.frame(res)
  res_df$gene <- convert_to_symbols(rownames(res_df), gene_map)
  
  # Classify genes by regulation status
  res_df$regulation <- "Not significant"
  res_df$regulation[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange > 1
  ] <- "Upregulated"
  res_df$regulation[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange < -1
  ] <- "Downregulated"
  
  # Create volcano plot
  p <- ggplot(
    res_df,
    aes(x = log2FoldChange, y = -log10(padj), color = regulation)
  ) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Upregulated" = "#D32F2F",
      "Downregulated" = "#1976D2",
      "Not significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw(base_size = 14) +
    labs(title = title)
  
  # Save plot
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# Generate volcano plot
plot_volcano(res_E85_vs_E65, "E8.5 vs E6.5", "volcano_E85_vs_E65.png")

# ==================================================================================
# SECTION 9: DAVID Functional Annotation
# ==================================================================================
# Load and process DAVID enrichment results (from DAVID website analysis)

# Load DAVID results for up- and down-regulated genes
DAVID_up <- read.csv("DAVID_up.csv", stringsAsFactors = FALSE)
DAVID_down <- read.csv("DAVID_down.csv", stringsAsFactors = FALSE)

# Filter for significant terms (FDR < 0.05)
sig_DAVID_up <- DAVID_up[DAVID_up$FDR < 0.05, ]
sig_DAVID_down <- DAVID_down[DAVID_down$FDR < 0.05, ]

# Add regulation direction labels
sig_DAVID_up$DGE <- "Upregulated"
sig_DAVID_down$DGE <- "Downregulated"

# Combine up and down results
sig_DAVID_combined <- suppressWarnings(rbind(sig_DAVID_up, sig_DAVID_down))

# ==================================================================================