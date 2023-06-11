# Load necessary packages
packages <- c('tidyverse', 
              'psych', 
              'GPArotation', 
              'EFA.dimensions')
lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})

# Define the upper and lower thresholds to filter items based on their factor loadings
upper_threshold = 0.3
lower_threshold = -0.3

# Load data from RDS
inpain_efa <- readRDS("inpain_efa.rds")
painfree_efa <- readRDS("painfree_efa.rds")
determine_factors <- readRDS("determine_factors.rds")
n_factors_inpain <- determine_factors$n_factors_inpain
n_factors_painfree <- determine_factors$n_factors_painfree

# Analyze loadings and create heatmaps
analyze_loadings <- function(loadings_list, group_name) {
  loadings_df <- do.call(rbind, lapply(loadings_list, function(x) {
    if(!is.null(x)){
      loadings = x$loadings_df
      if(all(is.na(loadings))){
        print(paste("All values in loadings are NA for group", group_name))
        return(NULL)
      }
      df = data.frame(Item = rownames(loadings), loadings)
      df <- df %>% pivot_longer(cols = -Item, names_to = "Factor", values_to = "Loading")
      return(df)
    }
  }))
  loadings_df$Group <- group_name
  print("Printing loadings_df after adding Group:")
  print(head(loadings_df))
  
  loadings_summary <- loadings_df %>%
    group_by(Item, Factor, Group) %>%
    summarise(mean = mean(Loading, na.rm = TRUE), 
              sd = sd(Loading, na.rm = TRUE), 
              lower = mean - qt(0.975, length(Loading)-1)*sd/sqrt(length(Loading)),  # Lower bound of 95% CI
              upper = mean + qt(0.975, length(Loading)-1)*sd/sqrt(length(Loading)),  # Upper bound of 95% CI
              .groups = 'drop')
  
  print("Printing loadings_summary after summarising:")
  print(head(loadings_summary))
  
  return(loadings_summary)
}

# Analyze loadings for both datasets
inpain_summary <- analyze_loadings(inpain_efa, "In Pain")
saveRDS(inpain_summary, "inpain_summary.rds")

painfree_summary <- analyze_loadings(painfree_efa, "Pain Free")
saveRDS(painfree_summary, "painfree_summary.rds")

# Function to generate heatmap based on loadings_summary for a specific factor
generate_heatmap <- function(loadings_summary, factor_name) {
  # Filter the summary for the specific factor
  loadings_summary <- loadings_summary %>% filter(Factor == factor_name)
  
  print("Printing loadings_summary after filtering for factor:")
  print(head(loadings_summary))
  
  # Generate the heatmap
  heatmap_plot <- ggplot(loadings_summary, aes(x = Group, y = Item, fill = mean)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste0(round(mean, 2), " (", round(lower, 2), "-", round(upper, 2), ")")), color = "black") +  # Modify this line
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Loading") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = paste("Average Factor Loadings for", factor_name), x = "Group", y = "Item")
  
  return(heatmap_plot)
}

# Combine the summaries
combined_summary <- rbind(inpain_summary, painfree_summary)
saveRDS(combined_summary, "combined_summary.rds")

# Get the list of factors
factors <- unique(combined_summary$Factor)
saveRDS(factors, "factors.rds")

# Generate a heatmap for each factor
for (factor in factors) {
  heatmap <- generate_heatmap(combined_summary, factor)
  print(heatmap)
}

generate_combined_scree_plot <- function(loadings_summary) {
  # Aggregate mean loadings for each factor and group
  scree_data <- loadings_summary %>%
    group_by(Factor, Group) %>%
    summarise(mean_loading = mean(mean, na.rm = TRUE),
              lower = mean(lower, na.rm = TRUE),  # Lower bound of 95% CI
              upper = mean(upper, na.rm = TRUE),  # Upper bound of 95% CI
              .groups = 'drop') %>%
    arrange(Group, desc(mean_loading)) %>%
    mutate(Factor = factor(Factor, levels = Factor))  # To keep the factors in descending order in the plot
  
  # Generate the scree plot
  scree_plot <- ggplot(scree_data, aes(x = Factor, y = mean_loading, color = Group)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    labs(title = "Scree Plot for Combined Data", x = "Factor", y = "Mean Loading") +
    theme_minimal()
  
  return(scree_plot)
}

# Generate scree plot for the combined summary
combined_scree_plot <- generate_combined_scree_plot(combined_summary)
print(combined_scree_plot)

average_residuals <- function(residmat_list) {
  # Get the number of items from the dimensions of the first residual matrix
  n_items <- dim(residmat_list[[1]])[1]
  
  # Initialize a matrix to hold the sums of residuals
  sum_residmat <- matrix(0, nrow = n_items, ncol = n_items)
  
  # Loop over the list of residual matrices
  for(i in seq_along(residmat_list)) {
    residmat <- residmat_list[[i]]
    
    # Check if residmat is a matrix and has the same dimensions as sum_residmat
    if(!is.matrix(residmat) || any(dim(residmat) != dim(sum_residmat))) {
      warning(paste0("Element ", i, " is not a matrix or doesn't have the same dimensions as the sum_residmat. Skipping this matrix."))
      next
    }
    
    # Add the absolute values of residuals to the sum matrix
    sum_residmat <- sum_residmat + abs(residmat)
  }
  
  # Divide by the number of iterations to get the average
  avg_residmat <- sum_residmat / length(residmat_list)
  
  # Return the average residual matrix
  return(avg_residmat)
}

# Define a function to convert a matrix to a data frame
matrix_to_df <- function(mat) {
  df <- as.data.frame(mat)
  df <- df %>% 
    rownames_to_column(var = "Item1") %>% 
    pivot_longer(cols = -Item1, names_to = "Item2", values_to = "AverageResidual")
  return(df)
}

# Define a function to average residuals across iterations
analyze_procrustes_results <- function(procrustes_list, group_name) {
  
  # Extract congruence values from each element of the procrustes_list and store in a vector
  congruence_vec <- sapply(procrustes_list, function(x) x$congruence)
  
  # Extract RMSR (Root Mean Square Residual) values from each element of the procrustes_list and store in a vector
  rmsr_vec <- sapply(procrustes_list, function(x) x$rmsr)
  
  # Create a summary data frame with the mean and standard deviation of congruence and RMSR values for the group
  summary_df <- data.frame(
    Group = group_name,
    Congruence_Mean = mean(congruence_vec, na.rm = TRUE),
    Congruence_SD = sd(congruence_vec, na.rm = TRUE),
    RMSR_Mean = mean(rmsr_vec, na.rm = TRUE),
    RMSR_SD = sd(rmsr_vec, na.rm = TRUE)
  )
  
  # Extract residual matrix from each element of the procrustes_list and store in a list
  residmat_list <- lapply(procrustes_list, function(x) x$residmat)
  
  # Call the function to calculate average residuals for the group
  avg_residmat <- average_residuals(residmat_list)
  
  # Convert the average residual matrix to a data frame using the matrix_to_df function
  # The resulting data frame has one row for each item pair, with columns for the two items and the average residual
  avg_residmat_df <- matrix_to_df(avg_residmat)
  
  # Add a column to the data frame to indicate the group
  avg_residmat_df$Group <- group_name
  
  # Append the data frame of average residuals to the summary data frame
  # This results in a single data frame that contains both the summary statistics and the average residuals for the group
  summary_df <- bind_rows(summary_df, avg_residmat_df)
  
  # Return a list that contains both data frames
  results <- list(
    summary = summary_df,
    avg_residmat = avg_residmat_df
  )
  return(results)
}

procrustes_results_inpain <- analyze_procrustes_results(inpain_efa, "inpain")
procrustes_results_painfree <- analyze_procrustes_results(painfree_efa, "painfree")

# Store the summary and average residual data frames in separate variables
inpain_summary <- procrustes_results_inpain$summary
inpain_avg_residmat <- procrustes_results_inpain$avg_residmat

painfree_summary <- procrustes_results_painfree$summary
painfree_avg_residmat <- procrustes_results_painfree$avg_residmat

# Filter items in a data frame based on factor loadings
filter_factor_items <- function(df1, df2, upper_threshold, lower_threshold) {
  
  # Define a helper function that will filter and rename columns in a data frame.
  filter_and_rename <- function(df, threshold, mean_name, ci_name) {
    # Filter the data frame based on the mean column values and threshold.
    # We also create a new CI column here that combines lower and upper values.
    df <- df %>% 
      filter(mean > upper_threshold | mean < lower_threshold) %>%
      mutate(CI = paste0("[", round(lower, 3), ", ", round(upper, 3), "]")) %>%
      select(Item, Factor, mean, CI)
    
    # Rename the mean and CI columns for clarity.
    names(df)[3:4] <- c(mean_name, ci_name)
    return(df)
  }
  
  # Define required columns.
  required_columns <- c("Item", "Factor", "mean", "lower", "upper")
  # Check if both df1 and df2 contain required columns, else throw an error.
  if (!all(required_columns %in% names(df1)) | !all(required_columns %in% names(df2))) {
    stop("Both df1 and df2 should contain columns: Item, Factor, mean, lower, upper")
  }
  
  # Use the helper function to filter and rename columns in df1 and df2.
  items1 <- filter_and_rename(df1, upper_threshold, "mean.df1", "CI.df1")
  items2 <- filter_and_rename(df2, lower_threshold, "mean.df2", "CI.df2")
  
  # Identify items that are common in both summaries.
  common_items <- items1 %>% inner_join(items2, by = c("Item", "Factor"))
  
  # If there are no common items, create an empty data frame with the same column names.
  if (nrow(common_items) == 0) {
    common_items <- data.frame(Item=character(), Factor=character(), mean.df1=numeric(), CI.df1=character(), mean.df2=numeric(), CI.df2=character())
  }
  
  # Return the common_items data frame.
  return(common_items)
}

# Filter items and save
common_items <- filter_factor_items(inpain_summary, painfree_summary, upper_threshold, lower_threshold)
print(common_items)
saveRDS(common_items, "common_items.rds")

# Get the list of common items
common_items_list <- unique(common_items$Item)
num_common_items <- (length(common_items_list))
saveRDS(num_common_items, "num_common_items.rds")

# Filter the original data frames
inpain.df <- inpain.df[, common_items_list]
painfree.df <- painfree.df[, common_items_list]

# Save each dataframe to an RDS file
saveRDS(inpain.df, "inpain_df.rds")
saveRDS(painfree.df, "painfree_df.rds")
