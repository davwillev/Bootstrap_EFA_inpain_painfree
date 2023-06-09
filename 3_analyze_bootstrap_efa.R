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

# Define the absolute threshold to filter items based on their factor loadings
absolute_threshold = 0.3

# Load data from RDS
iteration_number <- readRDS("iteration_number.rds")
inpain_efa <- readRDS("inpain_efa.rds")
painfree_efa <- readRDS("painfree_efa.rds")
determine_factors <- readRDS("determine_factors.rds")
inpain.df <- determine_factors$inpain_df
painfree.df <- determine_factors$painfree_df
n_factors_inpain <- determine_factors$n_factors_inpain
n_factors_painfree <- determine_factors$n_factors_painfree

analyze_loadings <- function(loadings_list, group_name) {
  # Loop over each element of the loadings_list and perform some transformations
  loadings_df <- do.call(rbind, lapply(loadings_list, function(x) {
    # If the element is not NULL, proceed with the transformations
    if(!is.null(x)){
      # Extract the loadings dataframe from the list element
      loadings = x$loadings_df
      # If all the values in loadings are NA, print a warning and return NULL
      if(all(is.na(loadings))){
        print(paste("All values in loadings are NA for group", group_name))
        return(NULL)
      }
      # Create a dataframe from the loadings dataframe with item names as a separate column
      df = data.frame(Item = rownames(loadings), loadings)
      # Convert the dataframe from wide format to long format with a separate row for each item and factor
      df <- df %>% pivot_longer(cols = -Item, names_to = "Factor", values_to = "Loading")
      # Return the transformed dataframe
      return(df)
    }
  }))
  
  # Add a new column to the dataframe to indicate the group
  loadings_df$Group <- group_name
  print("Printing loadings_df after adding Group:")
  print(head(loadings_df))
  
  # Summarise the loadings dataframe by calculating mean, sd, and 95% confidence interval for each item and factor
  loadings_summary <- loadings_df %>%
    group_by(Item, Factor, Group) %>%
    summarise(mean = mean(Loading, na.rm = TRUE), 
              sd = sd(Loading, na.rm = TRUE), 
              lower = mean - qt(0.975, length(Loading)-1)*sd/sqrt(length(Loading)),  # Lower bound of 95% CI
              upper = mean + qt(0.975, length(Loading)-1)*sd/sqrt(length(Loading)),  # Upper bound of 95% CI
              .groups = 'drop')
  
  print("Printing loadings_summary after summarising:")
  print(head(loadings_summary))
  
  # Return the summary dataframe
  return(loadings_summary)
}

# Analyze loadings for both datasets
inpain_summary <- analyze_loadings(inpain_efa, "In Pain")
saveRDS(inpain_summary, "inpain_summary.rds")

painfree_summary <- analyze_loadings(painfree_efa, "Pain Free")
saveRDS(painfree_summary, "painfree_summary.rds")

# Function to generate heatmap based on loadings_summary for a specific factor
generate_heatmap <- function(loadings_summary, factor_name, iteration_number) {
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
    labs(title = paste("Average Factor Loadings for", factor_name, "- Iteration", iteration_number), x = "Group", y = "Item")
  
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
  heatmap <- generate_heatmap(combined_summary, factor, iteration_number)
  
  # Create a filename using the iteration number and factor
  filename <- paste0("heatmap_", factor, "_", iteration_number, ".pdf")
  
  # Print and save the heatmap to a PDF
  print(heatmap)
  ggsave(filename, plot = heatmap, device = "pdf", width = 8.3, height = 11.7)
}

analyze_eigenvalues <- function(efa_list, group_name) {
  # Extract the eigenvalues from each EFA result and convert them to a data frame
  eigenvalues_list <- lapply(efa_list, function(x) {
    if(!is.null(x)) {
      eigenvalues <- x$eigenvalues
      if(all(is.na(eigenvalues))){
        print(paste("All values in eigenvalues are NA for group", group_name))
        return(NULL)
      }
      # Define factor names as MR1, MR2 etc. based on the length of eigenvalues
      factor_names <- factor(paste0("MR", seq_along(eigenvalues)), 
                             levels = paste0("MR", seq_len(length(eigenvalues))))
      
      df = data.frame(Factor = factor_names, Eigenvalue = eigenvalues)
      return(df)
    }
  })
  
  # Combine all eigenvalues dataframes into one
  eigenvalues_df <- do.call(rbind, eigenvalues_list)
  
  # Add a new column to the dataframe to indicate the group
  eigenvalues_df$Group <- group_name
  
  # Summarise the eigenvalues dataframe by calculating mean, sd, and 95% confidence interval for each factor
  eigenvalues_summary <- eigenvalues_df %>%
    group_by(Factor, Group) %>%
    summarise(mean = mean(Eigenvalue, na.rm = TRUE), 
              sd = sd(Eigenvalue, na.rm = TRUE), 
              lower = mean - qt(0.975, length(Eigenvalue)-1)*sd/sqrt(length(Eigenvalue)),  # Lower bound of 95% CI
              upper = mean + qt(0.975, length(Eigenvalue)-1)*sd/sqrt(length(Eigenvalue)),  # Upper bound of 95% CI
              .groups = 'drop')
  
  return(eigenvalues_summary)
}

# Define a function that creates a scree plot from mean eigenvalues
generate_scree_plot <- function(eigenvalues_summary, iteration_number) {
  # Generate the scree plot
  scree_plot <- ggplot(eigenvalues_summary, aes(x = Factor, y = mean, color = Group)) +
    geom_line(aes(group = Group)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6, color = "black", size = 0.3) +
    labs(title = paste("Scree Plot - Iteration", iteration_number), x = "Factor", y = "Eigenvalue") +
    theme_minimal()
  
  # Create a filename using the iteration number
  filename <- paste0("scree_plot_", iteration_number, ".pdf")
  
  # Print and save the plot to a PDF
  print(scree_plot)
  ggsave(filename, plot = scree_plot, device = "pdf")
  
  return(scree_plot)
}

# Analyze eigenvalues for both datasets
inpain_eigenvalues_summary <- analyze_eigenvalues(inpain_efa, "In Pain")
painfree_eigenvalues_summary <- analyze_eigenvalues(painfree_efa, "Pain Free")

# Combine both summaries into one dataframe
combined_eigenvalues_summary <- rbind(inpain_eigenvalues_summary, painfree_eigenvalues_summary)

# Generate scree plot
scree_plot <- generate_scree_plot(combined_eigenvalues_summary, iteration_number)

# Define a function to create a plot of mean factor loadings
generate_mean_loadings_plot <- function(loadings_summary, iteration_number) {
  # Aggregate mean loadings for each factor and group
  loadings_data <- loadings_summary %>%
    group_by(Factor, Group) %>%
    summarise(mean_loading = mean(mean, na.rm = TRUE),
              lower = mean(lower, na.rm = TRUE),  # Lower bound of 95% CI
              upper = mean(upper, na.rm = TRUE),  # Upper bound of 95% CI
              .groups = 'drop') %>%
    arrange(Group, desc(mean_loading))
  
  # Check for and handle duplicated factor levels
  loadings_data <- loadings_data %>%
    mutate(Factor = factor(Factor, levels = unique(Factor))) # The levels are in the order they first appear in scree_data
  
  # Generate the plot of mean factor loadings
  loadings_plot <- ggplot(loadings_data, aes(x = Factor, y = mean_loading, color = Group)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    labs(title = paste("Mean Factor Loadings for Combined Data - Iteration", iteration_number), x = "Factor", y = "Mean Loading") +
    theme_minimal()
  
  # Create a filename using the iteration number
  filename <- paste0("loadings_plot_", iteration_number, ".pdf")
  
  # Print and save the plot to a PDF
  print(loadings_plot)
  ggsave(filename, plot = loadings_plot, device = "pdf")
  
  return(loadings_plot)
}

# Generate the plot for the combined summary
combined_loadings_plot <- generate_mean_loadings_plot(combined_summary, iteration_number)

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
inpain_procrustes_summary <- procrustes_results_inpain$summary
inpain_avg_residmat <- procrustes_results_inpain$avg_residmat

painfree_procrustes_summary <- procrustes_results_painfree$summary
painfree_avg_residmat <- procrustes_results_painfree$avg_residmat

# Print the procrustes summaries
print(inpain_procrustes_summary)
print(painfree_procrustes_summary)
print(inpain_avg_residmat)
print(painfree_avg_residmat)

# Define a function to filter and rename columns in a data frame
filter_and_rename <- function(df, threshold) {
  df <- df %>%
    filter(abs(mean) >= threshold) %>%
    mutate(
      CI = paste0("[", round(lower, 3), ", ", round(upper, 3), "]"),
      direction = ifelse(mean < 0, "negative", "positive")
    )
  return(df)
}

# Define a function to check cross loadings
check_cross_loadings <- function(items, threshold) {
  items <- items %>%
    mutate(over_threshold = if_else(abs(mean) > threshold, 1, 0)) %>%
    group_by(Item) %>%
    summarise(total_cross_loadings = sum(over_threshold)) %>%
    filter(total_cross_loadings > 1)
  
  return(items)
}

# Define a function to get common factor items
filter_factor_items <- function(df1, df2, threshold) {
  # Filter and rename columns for each data frame
  items1 <- filter_and_rename(df1, threshold)
  items2 <- filter_and_rename(df2, threshold)
  
  # Check for cross loadings
  cross_loaded_items1 <- check_cross_loadings(items1, threshold)
  cross_loaded_items2 <- check_cross_loadings(items2, threshold)
  
  # Filter out cross loaded items
  items1 <- filter(items1, !Item %in% cross_loaded_items1$Item)
  items2 <- filter(items2, !Item %in% cross_loaded_items2$Item)
  
  # Get common items
  common_items <- inner_join(items1, items2, by = c("Item", "Factor"))
  
  # Renaming columns
  common_items <- common_items %>%
    rename(mean.df1 = mean.x,
           CI.df1 = CI.x,
           direction.df1 = direction.x,
           mean.df2 = mean.y,
           CI.df2 = CI.y,
           direction.df2 = direction.y)
  
  return(common_items)
}

# Get common items
common_items <- filter_factor_items(inpain_summary, painfree_summary, absolute_threshold)
print(common_items)
saveRDS(common_items, "common_items.rds")

# Function to save the number of common items and the common items list
save_common_items_info <- function(common_items) {
  common_items_list <- unique(common_items$Item)
  num_common_items <- length(common_items_list)
  
  # Save these to RDS
  saveRDS(num_common_items, "num_common_items.rds")
  saveRDS(common_items_list, "common_items_list.rds")
}

# Save common items info
save_common_items_info(common_items)

# Function to filter the original data frames
filter_original_dfs <- function(df1, df2, common_items) {
  common_items_list <- unique(common_items$Item)
  df1 <- df1[, common_items_list, drop = FALSE]
  df2 <- df2[, common_items_list, drop = FALSE]
  
  return(list(df1, df2))
}

# Filter the original data frames and save them
filtered_dfs <- filter_original_dfs(inpain.df, painfree.df, common_items)
inpain.df <- filtered_dfs[[1]]
painfree.df <- filtered_dfs[[2]]

# Save each data frame to an RDS file
saveRDS(inpain.df, "inpain_df.rds")
saveRDS(painfree.df, "painfree_df.rds")
