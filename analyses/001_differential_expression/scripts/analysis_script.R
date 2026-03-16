# Load required packages
library(DESeq2)
library(GEOquery)
library(pheatmap)
library(EnhancedVolcano)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Read counts file into R
counts <- read.delim("data/GSE81089_readcounts_featurecounts.tsv",row.names = 1)

# Extract samples names from counts
sample_titles <- colnames(counts)

# From the experiment design, I know that samples ending in T are tumour, and
# ending in N are normal...

# Get last character of each column name
last_char <- substring(colnames(counts), nchar(colnames(counts)))

# Assign condition based on last character
condition <- ifelse(last_char == "T", "NSCLC", "Normal")

# Convert to factor
condition <- factor(condition)

# Build colData
colData <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# Check for NAs
sum(is.na(counts))  # how many NA values

# Replace all NA values with 0
counts[is.na(counts)] <- 0

# Verify
sum(is.na(counts))  # should now be 0

# -----------------------------------------------------------------------

# Create DEseq2 dataset for analysis
dds <- DESeqDataSetFromMatrix(
  counts,
  colData,
  design = ~ condition
)

# Remove genes with very low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
results <- results(dds)
head(results)

# Use adjusted p-value < 0.05
sig_genes <- results[which(results$padj < 0.05), ]
nrow(sig_genes)

# ---------------------------------------------------------

# Convert ENSG to SYMBOL
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(results),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Add a symbol column to your results table
results$symbol <- gene_symbols
results$symbol <- ifelse(is.na(results$symbol), rownames(results), results$symbol)

# ---------------------------------------------------------------

# Basic volcano plot
EnhancedVolcano(
  results,
  lab = results$symbol,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,           # significance threshold
  FCcutoff = 1,             # log2 fold change threshold
  pointSize = 2.0,
  labSize = 3.0
)
