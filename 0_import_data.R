# Load necessary packages
packages <- c('tidyverse',
              'rstudioapi')
lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})

library(tidyverse)
library(rstudioapi)

# Get the directory of the current script
script_dir <- normalizePath(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import data from CSV files (ensure these are in same directory as current script)
inpain.df <- read_csv(file.path(script_dir, 'inpain.csv'))
painfree.df <- read_csv(file.path(script_dir, 'painfree.csv'))

# Define a function to handle missing data
handle_missing_data <- function(df) {
  df <- df %>% mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  return(df)
}

# Handle missing data
inpain.df <- handle_missing_data(inpain.df)
painfree.df <- handle_missing_data(painfree.df)

# Save each dataframe to an RDS file
saveRDS(inpain.df, "inpain_df.rds")
saveRDS(painfree.df, "painfree_df.rds")

# Check if any items common to both datasets have already been identified
if (file.exists("num_common_items.rds")) {
  prev_num_common_items <- readRDS("num_common_items.rds")
  inpain.df <- determine_factors$inpain_df
} else {
  prev_num_common_items <- 0
}

# Function to check whether the stopping condition is met
check_condition <- function(prev_num_common_items, num_common_items) {
  if(prev_num_common_items == num_common_items) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Ensure other R scripts are in same directory as current script
script_1 <- file.path(script_dir, "1_determine_number_factors.R")
script_2 <- file.path(script_dir, "2_run_bootstrap_efa.R")
script_3 <- file.path(script_dir, "3_analyze_bootstrap_efa.R")
script_4 <- file.path(script_dir, "4_multi_group_cfa.R")

# Loop over scripts 1 and 2 until condition is met
while(TRUE){
  source(script_1)
  source(script_2)
  source(script_3)
  
  num_common_items <- readRDS("num_common_items.rds")
  
  # Check if the condition is met
  if (check_condition(prev_num_common_items, num_common_items)) {
    break
  }
  prev_num_common_items <- num_common_items
}

# Once condition is met, run script 3
source(script_4)
