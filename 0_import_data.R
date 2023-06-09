# Load necessary packages
packages <- c('tidyverse')
lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})

# Define a function to handle missing data
handle_missing_data <- function(df) {
  df <- df %>% mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  return(df)
}

# Import datasets
inpain.df <- read_csv('/Users/davidevans/Library/CloudStorage/OneDrive-Personal/My Projects/UBham/Pain-agnositc Qs in painful subjects/EFA/GPT/inpain.csv')
painfree.df <- read_csv('/Users/davidevans/Library/CloudStorage/OneDrive-Personal/My Projects/UBham/Pain-agnositc Qs in painful subjects/EFA/GPT/painfree.csv')

# Handle missing data
inpain.df <- handle_missing_data(inpain.df)
painfree.df <- handle_missing_data(painfree.df)

# Save important variables into RDS file
saveRDS(list(inpain_df = inpain.df, 
             painfree_df = painfree.df),
        file = "import_data.rds")