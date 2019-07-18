#!/usr/bin/env Rscript

# Changing the default R behaviour that all strings are treated as factors when reading the table
options(stringsAsFactors = FALSE)

# Reading all command line arguments
args <- commandArgs(trailingOnly=TRUE)

# This script is used to plot a heatmap with reads associated to a particular transcript

# Loading heatmap3 library used to plot everything
library(heatmap3)

# Reading input table with all the normalized (TMM normalization) read counts for reads associated with a particular transcript
input_df_path <- args[1]
out_dendro <- args[2]
out_no_dendro <- args[3]

# Test path
#input_df_path <- "//tgac-user-data/tgac_hpc_work/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_novel_exons/transcript_positions/all_reads/all_reads_samples_normalized.txt"

# Reading table with read counts
input_df <- read.table(input_df_path, header = T, sep = "\t", comment.char = "")

# Removing splice pattern column from data frame (first column)
input_df <- input_df[,-1]

# Log10 transformation of the normalized read counts
input_df[,2:ncol(input_df)] = log10(input_df[,2:ncol(input_df)]+0.000001)

# Using transcript_id column as row names and removing this column from data frame
rownames(input_df) <- input_df$Transcript_ID
input_df <- input_df[,-1]

# Transposing the data frame (results in generation of matrix)
input_df <- t(input_df)

# Generating the plot with dendrogram
pdf(
  out_dendro, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
  )
heatmap3(
  t(input_df), # Transposed-transposed (i.e. native) matrix of normalized read counts
  margins = c(8.5, 0), # Margins sizes (in inches)
  showColDendro=T, # Printing out the column dendrogram
  showRowDendro=T, # Printing out the row dendrogram
  scale='none', # Do not scale anything
  distfun = dist, # Use 'dist' as a function to calculate the distance among rows (instead of 'function(x) as.dist(1 - cor(t(x), use = "pa"))')
  cexRow = 0.2 # Decrease size of row labels
  )
dev.off()

# Generating the plot without dendrogram (but with reordering rows and columns as with dendrogram)
pdf(
  out_no_dendro, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
heatmap3(
  t(input_df), # Transposed-transposed (i.e. native) matrix of normalized read counts
  margins = c(8.5, 0), # Margins sizes (in inches)
  showColDendro=F, # Do not show column dendrogram
  showRowDendro=F, # Do not show row dendrogram
  scale='none', # Do not scale anything
  distfun = dist, # Use 'dist' as a function to calculate the distance among rows (instead of 'function(x) as.dist(1 - cor(t(x), use = "pa"))')
  cexRow = 0.2 # Decrease size of row labels
  )
dev.off()

# Removing everything
rm(list = ls(all.names = T))
