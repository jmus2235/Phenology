# R script to concatenate GEE HLS EVI CSV files from 2016-2024
# Author: Generated for GEE data processing
# Date: 2025-06-24

library(dplyr)
library(readr)

wd = "C:/Users/jmusinsky/Documents/R_Scripts/HLS/data"
setwd(wd)

# Function to concatenate CSV files for a given pattern
concatenate_gee_files <- function(file_pattern = NULL, input_dir = ".", output_dir = ".") {
  
  # If no pattern provided, try to detect from files in directory
  if (is.null(file_pattern)) {
    # List all CSV files matching the general pattern
    all_files <- list.files(input_dir, pattern = "Mean_HLS_EVI_\\d{4}_8day_composite_.*\\.csv$", full.names = TRUE)
    
    if (length(all_files) == 0) {
      stop("No files found matching the expected pattern. Please check your directory or specify file_pattern.")
    }
    
    # Extract the suffix pattern from the first file
    first_file <- basename(all_files[1])
    # Extract everything after "composite_"
    suffix_match <- regmatches(first_file, regexpr("composite_.*\\.csv$", first_file))
    file_pattern <- gsub("composite_", "", gsub("\\.csv$", "", suffix_match))
    
    cat("Detected file pattern suffix:", file_pattern, "\n")
  }
  
  # Generate expected filenames for years 2016-2024
  years <- 2016:2025
  expected_files <- paste0("Mean_HLS_EVI_", years, "_8day_composite_", file_pattern, ".csv")
  expected_paths <- file.path(input_dir, expected_files)
  
  # Check which files exist
  existing_files <- expected_paths[file.exists(expected_paths)]
  missing_files <- expected_files[!file.exists(expected_paths)]
  
  if (length(missing_files) > 0) {
    cat("Warning: The following expected files are missing:\n")
    cat(paste(missing_files, collapse = "\n"), "\n\n")
  }
  
  if (length(existing_files) == 0) {
    stop("No files found for the specified pattern: ", file_pattern)
  }
  
  cat("Found", length(existing_files), "files to concatenate:\n")
  cat(paste(basename(existing_files), collapse = "\n"), "\n\n")
  
  # Read and combine all files
  cat("Reading and concatenating files...\n")
  combined_data <- data.frame()
  
  for (file_path in existing_files) {
    cat("Processing:", basename(file_path), "\n")
    
    # Read the CSV file
    current_data <- read_csv(file_path, show_col_types = FALSE)
    
    # Add year column for reference (extract from filename)
    year <- as.numeric(regmatches(basename(file_path), regexpr("\\d{4}", basename(file_path))))
    current_data$year <- year
    
    # Combine with existing data
    if (nrow(combined_data) == 0) {
      combined_data <- current_data
    } else {
      combined_data <- bind_rows(combined_data, current_data)
    }
  }
  
  # Sort by year and then by composite_date if available
  if ("composite_date" %in% names(combined_data)) {
    combined_data <- combined_data %>%
      arrange(year, composite_date)
  } else {
    combined_data <- combined_data %>%
      arrange(year)
  }
  
  # Generate output filename
  output_filename <- paste0("Mean_HLS_EVI_2016_2025_8day_composite_", file_pattern, ".csv")
  output_path <- file.path(output_dir, output_filename)
  
  # Write the combined data
  cat("Writing combined data to:", output_filename, "\n")
  write_csv(combined_data, output_path)
  
  # Print summary
  cat("\nSummary:\n")
  cat("- Total rows in combined file:", nrow(combined_data), "\n")
  cat("- Total columns:", ncol(combined_data), "\n")
  cat("- Years included:", paste(sort(unique(combined_data$year)), collapse = ", "), "\n")
  cat("- Output file:", output_path, "\n")
  
  return(combined_data)
}

# Example usage:
# Method 1: Auto-detect pattern from files in current directory
# combined_data <- concatenate_gee_files()

# Method 2: Specify the pattern explicitly (everything after "composite_")
combined_data <- concatenate_gee_files(file_pattern = "v2C_TOOL_FB_SH")

# Method 3: Specify custom input and output directories
# combined_data <- concatenate_gee_files(file_pattern = "BART_FB", 
#                                       input_dir = "path/to/input", 
#                                       output_dir = "path/to/output")

# Run the function (uncomment the line below and modify as needed)
# combined_data <- concatenate_gee_files(file_pattern = "BART_FB")