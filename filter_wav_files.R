# Clear workspace
rm(list = ls())

# Load necessary libraries
library(tidyverse)

# Define paths
source_folder <- "/Users/dan/Documents/birdnet_large_files/unfiltered_bacs_files"
destination_folder <- "/Users/dan/Documents/birdnet_large_files/filtered_bacs_files"
csv_file <- "data/val_bacs_for_filtering.csv"

# Read the CSV file
file_list <- read.csv(csv_file, stringsAsFactors = FALSE)

# Trim spaces from the filenames in the CSV
file_list$full_id <- str_trim(file_list$full_id)

# List all .wav files in the source folder
wav_files <- list.files(source_folder, pattern = "\\.wav$", full.names = TRUE)

# Filter .wav files based on the CSV list
filtered_files <- wav_files[basename(wav_files) %in% file_list$full_id]

# Create the destination folder if it doesn't exist
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

# Move the filtered .wav files to the new folder
file.copy(filtered_files, destination_folder, overwrite = TRUE)

# Optionally, you can remove the files from the source folder
# file.remove(filtered_files)
