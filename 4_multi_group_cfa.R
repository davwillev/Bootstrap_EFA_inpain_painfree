# Load necessary packages
packages <- c('tidyverse', 
              'lavaan')
lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})

# Load data from RDS files
inpain.df <- readRDS("inpain_df.rds")
painfree.df <- readRDS("painfree_df.rds")

inpain_summary <- readRDS("inpain_summary.rds")
painfree_summary <- readRDS("painfree_summary.rds")
combined_summary <- readRDS("combined_summary.rds")
factors <- readRDS("factors.rds")
inpain_efa <- readRDS("inpain_efa.rds")
painfree_efa <- readRDS("painfree_efa.rds")
combined_summary <- readRDS("combined_summary.rds")
num_common_items <- readRDS("num_common_items.rds")
common_items_list <- readRDS("common_items_list.rds")


# Add 'group' variable to each dataframe
inpain.df$group <- 'inpain'
painfree.df$group <- 'painfree'

# Combine the two dataframes
combined.df <- rbind(inpain.df, painfree.df)

# Function to create CFA model
create_cfa_model <- function(items, factor_name) {
  model <- paste0(factor_name, ' =~ ', paste(items, collapse = ' + '))
  return(model)
}

# Create CFA model
model <- create_cfa_model(common_items_list, 'F1')
print(model)

# Perform multi-group CFA
fit <- cfa(model, data = combined.df, group = "group")

# Check for measurement invariance across the two groups
fit.configural <- cfa(model, data = combined.df, group = "group")
fit.metric <- update(fit.configural, group.equal = "loadings")
fit.scalar <- update(fit.metric, group.equal = "intercepts")

# CFA evaluation
cfa_summary <- summary(fit, fit.measures = TRUE)
print(cfa_summary)

# Save CFA summary to text file
capture.output(summary(fit, fit.measures = TRUE), file = paste0(script_dir, "/cfa_summary.txt"))
