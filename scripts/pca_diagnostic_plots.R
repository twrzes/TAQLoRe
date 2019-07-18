#!/usr/bin/env Rscript

# Changing the default R behaviour that all strings are treated as factors when reading the table
options(stringsAsFactors = FALSE)

# Reading all command line arguments
args <- commandArgs(trailingOnly=TRUE)

# This script generates a number of diagnostic plots:
# out_sample_distance: plot with sample distance on log10-transformed normalized read values
# out_pc1_pc2: PCA plot with PC1 and PC2
# out_pc1_pc3: PCA plot with PC1 and PC3
# out_pc2_pc3: PCA plot with PC2 and PC3
# out_scree_plot: Plot with percentages of explained variance by each PC
# out_contrib_pc1: Plot with percentages of transcript contributions to PC1
# out_contrib_pc2: Plot with percentages of transcript contributions to PC2
# out_contrib_pc3: Plot with percentages of transcript contributions to PC3
# out_corr_pc1: XY plot with relation (correlation) of sample PC1 loadings and number of reads
# out_corr_pc2: XY plot with relation (correlation) of sample PC2 loadings and number of reads
# out_corr_pc3: XY plot with relation (correlation) of sample PC3 loadings and number of reads

# Input files
input_df_path <- args[1]
splicing_patterns_dir <- args[2]
distinct_string <- args[3]
barcode_to_sample_mappings_path <- args[4]
# Output files
out_sample_distance <- args[5]
out_pc1_pc2 <- args[6]
out_pc1_pc3 <- args[7]
out_pc2_pc3 <- args[8]
out_scree_plot <- args[9]
out_contrib_pc1 <- args[10]
out_contrib_pc2 <- args[11]
out_contrib_pc3 <- args[12]
out_corr_pc1 <- args[13]
out_corr_pc2 <- args[14]
out_corr_pc3 <- args[15]


# Importing barcode-to-sample-name mappings
sample_names_barcodes <- read.table(barcode_to_sample_mappings_path, header=F, sep="\t", comment.char = "")
# Generating a vector with paths to splicing patterns
file_list <- list.files(splicing_patterns_dir, distinct_string, full.names = T)
# Generate a named vector for length of lines
read_counts <- c()
# Iterate over sample names
for (i in sample_names_barcodes[,3]) {
  # Path to file to count lines
  file_to_count_lines <- c()
  # Obtain file name
  for (j in file_list) {
    # If sequencing run and barcode are in file name, append the path to the variable
    if (grepl(sample_names_barcodes[sample_names_barcodes[,3]==i,][,1], j, fixed = T) & grepl(sample_names_barcodes[sample_names_barcodes[,3]==i,][,2], j, fixed = T)) {
      file_to_count_lines <- c(file_to_count_lines, j)
    }
  }
  rm(j)
  # Read line count
  line_count <- length(readLines(file_to_count_lines))-1
  # Sub-vector of names
  sub_vector_to_append <- c(line_count)
  names(sub_vector_to_append) <- i
  # Append read counts
  read_counts <- c(read_counts, sub_vector_to_append)
  # Delete unnecessary variables
  rm(file_to_count_lines, line_count, sub_vector_to_append)
}
rm(i, file_list, sample_names_barcodes)

# Reading table with read counts
input_df <- read.table(input_df_path, header = T, sep = "\t", comment.char = "")

# Removing splice pattern column from data frame (first column)
input_df <- input_df[,-1]

# Log10 transformation of the normalized read counts
input_df[,2:ncol(input_df)] <- log10(input_df[,2:ncol(input_df)]+0.000001)

# Using transcript_id column as row names and removing this column from data frame
rownames(input_df) <- input_df$Transcript_ID
input_df <- input_df[,-1]

# Transposing the data frame (results in generation of matrix)
input_df <- t(input_df)

# Calculating distance matrix for samples
distance_matrix <- as.matrix(dist(input_df))

library(heatmap3)

# Plotting the heatmap of sample distances
pdf(
  out_sample_distance, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
heatmap3(
  distance_matrix, # Transposed-transposed (i.e. native) matrix of normalized read counts
  margins = c(10, 0.5), # Margins sizes (in inches)
  showColDendro=T, # Printing out the column dendrogram
  showRowDendro=T, # Printing out the row dendrogram
  scale='none', # Do not scale anything
  distfun = dist,  # Use 'dist' as a function to calculate the distance among rows (instead of 'function(x) as.dist(1 - cor(t(x), use = "pa"))')
  col = colorRampPalette(c("white",
                     "navy"))(1024)
)
dev.off()
rm(distance_matrix)
library(ggplot2)
library(ggfortify)

# Creating annotations for PCA plots
annotated_df <- as.data.frame(input_df)

# String with all information
string_to_split <- rownames(annotated_df)

# Sequencing run
sequencing_run <- c()
# Sequencing run is the first element in the string (underscore-delimited)
for (i in string_to_split) sequencing_run <- c(sequencing_run, strsplit(i, '_')[[1]][1])

# Sample name
sample_name <- c()
# Sample name is the second element in the string (underscore-delimited)
for (i in string_to_split) sample_name <- c(sample_name, strsplit(i, '_')[[1]][2])

# Brain region
brain_region <- c()
# Brain region is the third element in the string (underscore-delimited)
for (i in string_to_split) brain_region <- c(brain_region, strsplit(i, '_')[[1]][3])
rm(i, string_to_split)

# Adding the month of the sequencing run
annotated_df['sequencing_run'] = as.factor(sequencing_run)
# Adding sample names
annotated_df['sample_name'] = as.factor(sample_name)
# Adding sequenced regions of a brain
annotated_df['brain_region'] = as.factor(brain_region)

# Creating named vector for alpha values in PCA plots
alpha_values_for_sequencing_run <- c()

# Starting alpha_value
alpha_value_decrement <- 1/length(unique(sequencing_run))

#################### TEST!!!!!!!!!!!!!!!!
for (i in 1:length(unique(sequencing_run))) {
  alpha_value_sample <- 1-((i-1)*alpha_value_decrement)
  names(alpha_value_sample) <- unique(sequencing_run)[i]
  alpha_values_for_sequencing_run <- c(alpha_values_for_sequencing_run, alpha_value_sample)
  rm(alpha_value_sample)
}
rm(i, sequencing_run, sample_name, brain_region, alpha_value_decrement)




# Doing PCA
pca_df <- as.data.frame(input_df)
p <- prcomp(pca_df, scale=T)


# For loop for all PC variants (PC1 vs. PC2, PC1 vs. PC3, PC2 vs. PC3)

for (i in c("PC1_vs_PC2", "PC1_vs_PC3", "PC2_vs_PC3")) {
  if (i=="PC1_vs_PC2") {
    plot.data <- ggplot2::fortify(p, data = annotated_df)
    variance_explained_one <- p$sdev[1]^2 / sum(p$sdev^2)
    variance_explained_two <- p$sdev[2]^2 / sum(p$sdev^2)
    label_x <- paste0("PC1", " (", round(variance_explained_one * 100, 2), "%)")
    label_y <- paste0("PC2", " (", round(variance_explained_two * 100, 2), "%)")
    lam <- p$sdev[c(1, 2)]
    lam <- lam * sqrt(nrow(plot.data))
    plot.data[, c("PC1", "PC2")] <- t(t(plot.data[, c("PC1", "PC2")]) / lam)
    pdf(
      out_pc1_pc2, # Name of the file
      width=11.7, # Width in inches
      height=8.3 # Height in inches
    )
    print(ggplot(plot.data,aes(x=PC1, y=PC2, color=brain_region, shape=sample_name)) + theme_minimal() + geom_point(size=3, aes(fill=brain_region, alpha=sequencing_run)) + geom_point(size=3) + scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25)) + scale_alpha_manual(values=alpha_values_for_sequencing_run) + xlab(label_x) + ylab(label_y) + theme(legend.text=element_text(size=8), legend.title=element_text(size=8), axis.title = element_text(size=8), axis.text = element_text(size=8)))
    dev.off()
    rm(plot.data, variance_explained_one, variance_explained_two, label_x, label_y, lam)
  } else if (i=="PC1_vs_PC3") {
    plot.data <- ggplot2::fortify(p, data = annotated_df)
    variance_explained_one <- p$sdev[1]^2 / sum(p$sdev^2)
    variance_explained_two <- p$sdev[3]^2 / sum(p$sdev^2)
    label_x <- paste0("PC1", " (", round(variance_explained_one * 100, 2), "%)")
    label_y <- paste0("PC3", " (", round(variance_explained_two * 100, 2), "%)")
    lam <- p$sdev[c(1, 3)]
    lam <- lam * sqrt(nrow(plot.data))
    plot.data[, c("PC1", "PC3")] <- t(t(plot.data[, c("PC1", "PC3")]) / lam)
    pdf(
      out_pc1_pc3, # Name of the file
      width=11.7, # Width in inches
      height=8.3 # Height in inches
    )
    print(ggplot(plot.data,aes(x=PC1, y=PC3, color=brain_region, shape=sample_name)) + theme_minimal() + geom_point(size=3, aes(fill=brain_region, alpha=sequencing_run)) + geom_point(size=3) + scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25)) + scale_alpha_manual(values=alpha_values_for_sequencing_run) + xlab(label_x) + ylab(label_y) + theme(legend.text=element_text(size=8), legend.title=element_text(size=8), axis.title = element_text(size=8), axis.text = element_text(size=8)))
    dev.off()
    rm(plot.data, variance_explained_one, variance_explained_two, label_x, label_y, lam)
  } else if (i=="PC2_vs_PC3") {
    plot.data <- ggplot2::fortify(p, data = annotated_df)
    variance_explained_one <- p$sdev[2]^2 / sum(p$sdev^2)
    variance_explained_two <- p$sdev[3]^2 / sum(p$sdev^2)
    label_x <- paste0("PC2", " (", round(variance_explained_one * 100, 2), "%)")
    label_y <- paste0("PC3", " (", round(variance_explained_two * 100, 2), "%)")
    lam <- p$sdev[c(2, 3)]
    lam <- lam * sqrt(nrow(plot.data))
    plot.data[, c("PC2", "PC3")] <- t(t(plot.data[, c("PC2", "PC3")]) / lam)
    pdf(
      out_pc2_pc3, # Name of the file
      width=11.7, # Width in inches
      height=8.3 # Height in inches
    )
    print(ggplot(plot.data,aes(x=PC2, y=PC3, color=brain_region, shape=sample_name)) + theme_minimal() + geom_point(size=3, aes(fill=brain_region, alpha=sequencing_run)) + geom_point(size=3) + scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25)) + scale_alpha_manual(values=alpha_values_for_sequencing_run) + xlab(label_x) + ylab(label_y) + theme(legend.text=element_text(size=8), legend.title=element_text(size=8), axis.title = element_text(size=8), axis.text = element_text(size=8)))
    dev.off()
    rm(plot.data, variance_explained_one, variance_explained_two, label_x, label_y, lam)
  } else print("WRONG!!!!")
}

rm(annotated_df, alpha_values_for_sequencing_run)
# Calculating PCA using different function
library("FactoMineR")
p2 <- PCA(pca_df, graph=FALSE)
rm(pca_df)
# Visualizing PCA object
library("factoextra")

# Plotting percentages of explained variance by each PC
pdf(
  out_scree_plot, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
fviz_screeplot(
  p2, # Name of PCA object
  ncp=10 # Name of principal components to plot
)
dev.off()

# Plotting percentages of transcript contributions to PC1
pdf(
  out_contrib_pc1, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
fviz_contrib(
  p2, # Name of the object with PCA results
  choice = "var", # Plotting transcript contriubutions (i.e. not sample contributions)
  axes = 1, # Name of PC to plot
  top=20 # Number of contributions to plot
)
dev.off()

# Plotting percentages of transcript contributions to PC2
pdf(
  out_contrib_pc2, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
fviz_contrib(
  p2, # Name of the object with PCA results
  choice = "var", # Plotting transcript contriubutions (i.e. not sample contributions)
  axes = 2, # Name of PC to plot
  top=20 # Number of contributions to plot
)
dev.off()

# Plotting percentages of transcript contributions to PC3
pdf(
  out_contrib_pc3, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
fviz_contrib(
  p2, # Name of the object with PCA results
  choice = "var", # Plotting transcript contriubutions (i.e. not sample contributions)
  axes = 3, # Name of PC to plot
  top=20 # Number of contributions to plot
)
dev.off()

rm(p2)

# Extract the sample loadings (i.e. PC values for all sample)
sample_loadings <- as.data.frame(p$x)
rm(p)
# Add the number of reads to a data frame
sample_loadings['num_reads'] <- unname(read_counts)
rm(read_counts)

# Plot correlation between PC1 and number of reads
to_plot <- sample_loadings[,c(1,ncol(sample_loadings))]
pdf(
  out_corr_pc1, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
ggplot(
  to_plot, # Name of the data object
  aes(
    x=PC1, # Name of the variable used for x axis
    y=num_reads # Name of the variable used for y axis
    )
) +
geom_point() + # Add points to the plot
geom_smooth(method=lm) + # Add correlation curve to the plot (with standard error confidence intervals)
annotate( # Add correlation coefficient to the plot
  "text",
  colour = "red", # Colour of the text
  x = 5, # Position of the text in x axis
  y = 2400, # Position of the text in y axis
  label = paste("correlation =", cor(to_plot$PC1, to_plot$num_reads)) # What to print (i.e. correlation={correlation  coefficient})
)
dev.off()

# Plot correlation between PC2 and number of reads
to_plot <- sample_loadings[,c(2,ncol(sample_loadings))]
pdf(
  out_corr_pc2, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
ggplot(
  to_plot, # Name of the data object
  aes(
    x=PC2, # Name of the variable used for x axis
    y=num_reads # Name of the variable used for y axis
  )
) +
  geom_point() + # Add points to the plot
  geom_smooth(method=lm) + # Add correlation curve to the plot (with standard error confidence intervals)
  annotate( # Add correlation coefficient to the plot
    "text",
    colour = "red", # Colour of the text
    x = 0, # Position of the text in x axis
    y = 2400, # Position of the text in y axis
    label = paste("correlation =", cor(to_plot$PC2, to_plot$num_reads)) # What to print (i.e. correlation={correlation  coefficient})
  )
dev.off()

# Plot correlation between PC2 and number of reads
to_plot <- sample_loadings[,c(3,ncol(sample_loadings))]
pdf(
  out_corr_pc3, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
ggplot(
  to_plot, # Name of the data object
  aes(
    x=PC3, # Name of the variable used for x axis
    y=num_reads # Name of the variable used for y axis
  )
) +
  geom_point() + # Add points to the plot
  geom_smooth(method=lm) + # Add correlation curve to the plot (with standard error confidence intervals)
  annotate( # Add correlation coefficient to the plot
    "text",
    colour = "red", # Colour of the text
    x = 1, # Position of the text in x axis
    y = 2400, # Position of the text in y axis
    label = paste("correlation =", cor(to_plot$PC3, to_plot$num_reads)) # What to print (i.e. correlation={correlation  coefficient})
  )
dev.off()

rm(list = ls(all.names = T))
