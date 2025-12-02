# ==================================================================================
# RNA-seq DGE Analysis: PCA + Labeled Volcano Plot + EnrichR
# ==================================================================================

# Load packages
library(rtracklayer)      #importing GTF/GFF files and BigWig files
library(GenomicFeatures)  #creating TxDb objects from gene annotations
library(GenomicRanges)    #handling genomic ranges and intervals
library(limma)            #differential expression analysis
library(tidyverse)        #data manipulation and visualization
library(enrichR)          #pathway enrichment analysis

# ==================================================================================
# STEP 1: Load GTF, Convert Gene ID -> Gene Symbols
# ==================================================================================

#create a transcript database object from GTF annotation file
txdb <- makeTxDbFromGFF("data/Mus_musculus.GRCm39.115.gtf")

#extract genomic coordinates from the database
genes <- genes(txdb)

#import the full GTF file with all metadata
gtf <- import("data/Mus_musculus.GRCm39.115.gtf")

#filter the gtf to keep only gene-level entries 
gene_info <- gtf[gtf$type == "gene"]

#create a map that connects the Ensembl IB -> gene symbols 
gene_map <- data.frame(
  ensembl_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  stringsAsFactors = FALSE
)

#remove the duplicate Ensembl IDs
gene_map <- gene_map[!duplicated(gene_map$ensembl_id), ]

#set the Ensembl IDs as rownames for easy lookup
rownames(gene_map) <- gene_map$ensembl_id

#convert Ensembl IDs to gene symbols
convert_to_symbols <- function(ensembl_ids, gene_map) {
  #remove version numbers from Ensembl IDs (e.g., ENSMUSG00000000001.5 -> ENSMUSG00000000001)
  ensembl_ids_clean <- sub("\\.\\d+$", "", ensembl_ids)
  
  #match Ensembl IDs to gene symbols using the gene_map
  symbols <- gene_map$gene_name[match(ensembl_ids_clean, gene_map$ensembl_id)]
  
  #for IDs without symbols, keep the original Ensembl ID
  symbols[is.na(symbols) | symbols == ""] <- ensembl_ids[is.na(symbols) | symbols == ""]
  
  #handle duplicate gene symbols by appending the Ensembl ID (so there are two distinct symbols)
  if(any(duplicated(symbols))) {
    dup_symbols <- unique(symbols[duplicated(symbols)])
    for(sym in dup_symbols){
      idx <- which(symbols == sym)
      symbols[idx] <- paste0(sym, "_", ensembl_ids[idx])
    }
  }
  
  #prepend "GENE_" to symbols that are purely numeric (to avoid issues later on)
  numeric_symbols <- grepl("^[0-9]+$", symbols)
  symbols[numeric_symbols] <- paste0("GENE_", symbols[numeric_symbols])
  
  return(symbols)
}

# ==================================================================================
# STEP 2: Import BigWig Files and Extract FPKM
# ==================================================================================

#find all the bigwig files in the rna_bw directory
bw_files <- list.files(path="data/rna_bw", pattern="\\.bw$", full.names=TRUE)

#extract the sample names from the files (for the FPKM matrix)
sample_names <- gsub("\\.bw$", "", basename(bw_files))

#create an empty matrix to store FPKM values -> rows are genes, columns are samples
fpkm_matrix <- matrix(0, nrow=length(genes), ncol=length(bw_files))
rownames(fpkm_matrix) <- names(genes)
colnames(fpkm_matrix) <- sample_names

#filling in the empty FPKM matrix
#loop through every bigwig
for(i in seq_along(bw_files)){
  #message("Processing sample ", i, " of ", length(bw_files), ": ", sample_names[i])
  
  #import the bigwig coverage data as an RleList (run-length encoded list)
  bw_cov <- import(bw_files[i], format="BigWig", as="RleList")
  
  #find the chromosomes in the gene annotation and the bigwig
  common_chroms <- intersect(seqlevels(genes), names(bw_cov))
  
  #subset the genes to only include the common chromosomes
  genes_subset <- keepSeqlevels(genes, common_chroms, pruning.mode="coarse")
  seqlevels(genes_subset) <- common_chroms
  
  #now, subset the bigwig coverage to onlyinclude common chromosomes
  bw_cov <- bw_cov[common_chroms]
  
  #split the genes by chromosome (for efficient processing)
  genes_by_chr <- split(genes_subset, seqnames(genes_subset))
  
  #create a Views object to map the coverage to gene regions
  gene_views <- Views(bw_cov, ranges(genes_by_chr))
  
  #calculate mean coverage (FPKM) for each gene
  gene_means_list <- lapply(gene_views, viewMeans)
  gene_means <- unlist(gene_means_list)
  
  #extract gene names from the split gene list
  all_gene_names <- unlist(lapply(genes_by_chr, names))
  names(gene_means) <- all_gene_names
  
  #store the calculated FPKM values in the matrix
  fpkm_matrix[all_gene_names, i] <- gene_means
}

#convert rownames from Ensembl IDs to gene symbols
rownames(fpkm_matrix) <- convert_to_symbols(sub("\\.\\d+$", "", rownames(fpkm_matrix)), gene_map)

# ==================================================================================
# STEP 3: Load Metadata & Adjust to Fit Limma
# ==================================================================================

#read in the metadata
coldata <- read.csv("metadata.csv", stringsAsFactors=FALSE)

#set the sample name as rownames for matching to FPKM
rownames(coldata) <- coldata$sample

#convert stage, condition, and lineage to factors 
coldata$stage <- factor(coldata$stage, levels=c("E6.5","E8.5"))
coldata$condition <- factor(coldata$condition)
coldata$lineage <- factor(coldata$lineage)

#reorder FPKM matrix columns to match the order in metadata
fpkm_matrix <- fpkm_matrix[, rownames(coldata)]

# ==================================================================================
# STEP 4: Filter the Genes and then Log2 Transform
# ==================================================================================

#filter out lowly expressed genes (keep genes with mean FPKM >= 1)
keep <- rowMeans(fpkm_matrix) >= 1
fpkm_filtered <- fpkm_matrix[keep, ]

#Log2 transform FPKM values (adding 0.1 to avoid log(0))
log2_fpkm <- log2(fpkm_filtered + 0.1)

#remove genes with NA or infinite values after transformation
log2_fpkm <- log2_fpkm[!apply(log2_fpkm, 1, function(x) any(is.na(x)|is.infinite(x))), ]

# ==================================================================================
# STEP 5: PCA 
# ==================================================================================

#perform the pca 
pca_res <- prcomp(t(log2_fpkm), scale.=FALSE)

#calculate variance based on pca
pca_var <- round(100*(pca_res$sdev^2 / sum(pca_res$sdev^2)),1)

#create a df for pc1 and pc2 for each sample
pca_data <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

#merge the pca coord with the metadata we loaded in step 3
pca_data <- merge(pca_data, coldata, by.x="Sample", by.y="row.names", all.x=TRUE)

#now, create a pca plot colored by stage and shaped by condition (should all be the same condition)
pca_plot <- ggplot(pca_data, aes(x=PC1, y=PC2, color=stage, shape=condition)) +
  geom_point(size=4, alpha=0.8) +
  labs(title="PCA: PC1 vs PC2",
       x=paste0("PC1 (",pca_var[1],"%)"),
       y=paste0("PC2 (",pca_var[2],"%)")) +
  theme_bw(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"))

#save the pca plot as a png file 
ggsave("PCA_PC1_PC2.png", pca_plot, width=10, height=8, dpi=300)

# ==================================================================================
# STEP 6: Diff Expression via limma
# ==================================================================================

#create a design matrix with stage as the variable of interest
# E6.5 is the reference level (intercept), E8.5 coefficient represents the difference
design <- model.matrix(~stage, data=coldata)

#fit linear model for each gene
fit <- lmFit(log2_fpkm, design)

#apply empirical Bayes smoothing to standard errors (required for FPKM)
#trend=TRUE assumes a mean-variance trend, robust=TRUE downweights outlier genes
fit <- eBayes(fit, trend=TRUE, robust=TRUE)

#extract results for the E8.5 vs E6.5 comparison
res <- topTable(fit, coef="stageE8.5", number=Inf, adjust.method="BH", sort.by="none")

#now, fix up the format of the dge
#add gene names as a column
res$gene <- rownames(res)

#rename columns to standard nomenclature
colnames(res)[colnames(res)=="logFC"] <- "log2FoldChange"
colnames(res)[colnames(res)=="adj.P.Val"] <- "padj"
colnames(res)[colnames(res)=="P.Value"] <- "pvalue"

#identify significantly upregulated genes (padj < 0.05, log2FC > 1) (value from the paper)
sig_UP <- res[res$padj<0.05 & res$log2FoldChange>1, "gene"]

#identify significantly downregulated genes (padj < 0.05, log2FC < -1)
sig_DOWN <- res[res$padj<0.05 & res$log2FoldChange< -1, "gene"]

#save gene lists to text files for EnrichR analysis
writeLines(sig_UP, "E85_vs_E65_UP.txt")
writeLines(sig_DOWN, "E85_vs_E65_DOWN.txt")

# ==================================================================================
# STEP 7: Labeled Volcano Plot from DGE
# ==================================================================================

#create a copy of results for volcano plot
volcano_data <- res

#classify the genes into UP, DOWN, or not significant (NS)
volcano_data$significance <- "NS"
volcano_data$significance[volcano_data$padj<0.05 & volcano_data$log2FoldChange>1] <- "UP"
volcano_data$significance[volcano_data$padj<0.05 & volcano_data$log2FoldChange< -1] <- "DOWN"

#set up the factor levels for proper ordering in legend
volcano_data$significance <- factor(volcano_data$significance, levels=c("DOWN","NS","UP"))

#calculate -log10(adjusted p-value) for y-axis
volcano_data$neg_log10_padj <- -log10(volcano_data$padj)

#handle infinite values (very small p-values) by capping at max + 10
max_log_p <- max(volcano_data$neg_log10_padj[!is.infinite(volcano_data$neg_log10_padj)], na.rm=TRUE)
volcano_data$neg_log10_padj[is.infinite(volcano_data$neg_log10_padj)] <- max_log_p + 10

#select top 10 upregulated genes by absolute log2FC for labeling
top_up <- volcano_data %>% filter(significance=="UP") %>% arrange(desc(abs(log2FoldChange))) %>% slice_head(n=10)

#select top 10 downregulated genes by absolute log2FC for labeling
top_down <- volcano_data %>% filter(significance=="DOWN") %>% arrange(desc(abs(log2FoldChange))) %>% slice_head(n=10)

#combine top genes for labeling on plot
top_genes <- rbind(top_up, top_down)

#create volcano plot
volcano_plot <- ggplot(volcano_data, aes(x=log2FoldChange, y=neg_log10_padj, color=significance)) +
  geom_point(alpha=0.6, size=1.5) +
  # Color scheme: blue for downregulated, grey for NS, red for upregulated
  scale_color_manual(values=c("DOWN"="#2166ac","NS"="grey60","UP"="#b2182b")) +
  # Add vertical lines at log2FC thresholds (±1)
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey30", alpha=0.5) +
  # Add horizontal line at significance threshold (padj = 0.05)
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey30", alpha=0.5) +
  # Add gene labels for top differentially expressed genes
  geom_text(data=top_genes, aes(label=gene), size=3, hjust=-0.1, vjust=0.5, color="black", show.legend=FALSE) +
  labs(title="Volcano Plot: E8.5 vs E6.5 (Top Genes Labeled)",
       x="log2 Fold Change",
       y="-log10(adjusted p-value)",
       color="Regulation") +
  theme_bw(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"))

#save volcano plot as PNG file
ggsave("volcano_plot_E85_vs_E65_labeled.png", volcano_plot, width=12, height=10, dpi=300)

# ==================================================================================
# STEP 8: EnrichR Pathway Enrichment
# ==================================================================================

#get list of all available EnrichR databases
dbs <- listEnrichrDbs()

#filter to include only mouse-specific databases
dbs_selection <- unique(dbs$libraryName)[grep("Mouse", unique(dbs$libraryName))]

#function to plot top enriched pathways from EnrichR results
plot_top_enrichments <- function(enriched_list, top_n = 20, filename_plot = NULL, filename_csv = NULL, plot_title = NULL) {
  
  #combine results from all databases into a single data frame
  df <- bind_rows(lapply(names(enriched_list), function(db) {
    res <- enriched_list[[db]]
    # Only include databases with results
    if (!is.null(res) && nrow(res) > 0) { 
      res$Database <- db
      return(res)
    }
  }))
  
  #check if any enriched terms were found
  if (is.null(df) || nrow(df) == 0) { 
    message("No enriched terms found.")
    return(NULL)
  }
  
  #select top N most significant terms
  df_top <- df %>% arrange(Adjusted.P.value) %>% slice_head(n = top_n)
  
  #save results to CSV if filename provided
  if (!is.null(filename_csv)) write.csv(df_top, filename_csv, row.names = FALSE)
  
  #shorten long pathway names to 60 characters for better display
  df_top$Term_short <- sapply(df_top$Term, function(x) {
    if (nchar(x) > 60) paste0(substr(x, 1, 57), "...") else x
  })
  
  #set plot title (use default if not provided)
  title_text <- ifelse(is.null(plot_title), paste("Top", top_n, "Enriched Mouse Terms"), plot_title)
  
  #create horizontal bar plot of enriched pathways
  p <- ggplot(df_top, aes(x = reorder(Term_short, -log10(Adjusted.P.value)), 
                          y = -log10(Adjusted.P.value), fill = Database)) +
    geom_bar(stat = "identity") + 
    coord_flip() +  # Horizontal bars for easier reading
    theme_bw(base_size = 12) +
    labs(title = title_text, x = "", y = "-log10(Adjusted P-value)") + 
    theme(
      legend.position = "left",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 10)
    )
  
  #save plot if filename provided
  if (!is.null(filename_plot)) {
    ggsave(filename = filename_plot, plot = p, width = 16, height = 12, dpi = 300)
  }
  
  return(p)
}

#store upregulated gene list
mouse_genes_UP <- sig_UP

#store downregulated gene list
mouse_genes_DOWN <- sig_DOWN

#run enrichment analysis for upregulated genes
if(length(mouse_genes_UP) > 0) {
  message("Running enrichment analysis for upregulated genes...")
  # Submit genes to EnrichR using all mouse databases
  enriched_UP <- enrichr(mouse_genes_UP, dbs_selection)
  # Plot and save top 20 enriched pathways
  plot_top_enrichments(enriched_UP, top_n=20,
                       filename_plot="enrichr_up_top20_enrichment_large.png",
                       filename_csv="enrichr_up_top20_enrichment.csv",
                       plot_title="Top 20 Enriched Pathways (UP-regulated Genes, E8.5 vs E6.5)")
} else { 
  message("No upregulated genes found - skipping enrichment analysis") 
}

#run enrichment analysis for downregulated genes
if(length(mouse_genes_DOWN) > 0) {
  message("Running enrichment analysis for downregulated genes...")
  # Submit genes to EnrichR using all mouse databases
  enriched_DOWN <- enrichr(mouse_genes_DOWN, dbs_selection)
  # Plot and save top 20 enriched pathways
  plot_top_enrichments(enriched_DOWN, top_n=20,
                       filename_plot="enrichr_down_top20_enrichment_large.png",
                       filename_csv="enrichr_down_top20_enrichment.csv",
                       plot_title="Top 20 Enriched Pathways (DOWN-regulated Genes, E8.5 vs E6.5)")
} else { 
  message("No downregulated genes found - skipping enrichment analysis") 
}