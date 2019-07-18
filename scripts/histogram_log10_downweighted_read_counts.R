#!/usr/bin/env Rscript

# Changing the default R behaviour that all strings are treated as factors when reading the table
options(stringsAsFactors = FALSE)

# Reading all command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Arguments
# Input file with table with downweighted counts
input_file_path <- args[1]
# Sample name string
sample_name <- args[2]
# Path to the PDF file with output plot
output_plot_path <- args[3]

# Loading table with down-weighted counts
input_table <- read.table(file = input_file_path, header = T, sep = "\t", comment.char = "")
# Leaving only transcript splice pattern and specific sample name column
input_table <- input_table[,c("X.Transcript_splicing_pattern", sample_name)]
# Log10 downweighted counts
input_table[,sample_name] <- log10(input_table[,sample_name]+0.0001)

# Loading required library
library(ggplot2)

# Writing the plot to PDF file
pdf(
  output_plot_path, # Name of the file
  width=11.7, # Width in inches
  height=8.3 # Height in inches
)
print(ggplot( # Using 'print' function to print to the ggplot2
  input_table, # Data table name
  aes_string( # Using aes_string - I can use variables with strings
    x = paste0("reorder(", "X.Transcript_splicing_pattern", ", -", sample_name, ")"), # Using reorder function to reorder x values by descending y values
    y = sample_name)) # Name of y axis variable
  + geom_bar(stat = "identity") # Generating barplot
  + theme_classic() # Using classic theme
  + theme( # Removing x axis labels
    axis.title.x=element_blank(), # Remove x axis title
    axis.text.x=element_blank(), # Remove x axis text
    axis.ticks.x=element_blank())) # Remove x axits
dev.off() # Closing file handle

## END OF THE SCRIPT
