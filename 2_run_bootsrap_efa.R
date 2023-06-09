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

# Set seed for reproducibility
set.seed(123)

# Define the number of bootstrap iterations
n_iterations <- 1000

# Define the upper and lower thresholds to filter items based on their factor loadings
upper_threshold = 0.6
lower_threshold = -0.6

# Load variables from last script
determine_factors <- readRDS("determine_factors.rds")

n_factors_inpain <- determine_factors$n_factors_inpain
n_factors_painfree <- determine_factors$n_factors_painfree
inpain_suggested_factors <- determine_factors$inpain_suggested_factors
painfree_suggested_factors <- determine_factors$painfree_suggested_factors
inpain.df <- determine_factors$inpain_df
painfree.df <- determine_factors$painfree_df

# Run EFA with the fixed number of factors to establish the factor structure
efa_full_inpain <- fa(inpain.df, nfactors = n_factors_inpain, rotate = "oblimin", fm = "minres")
efa_full_painfree <- fa(painfree.df, nfactors = n_factors_painfree, rotate = "oblimin", fm = "minres")

# Function to perform EFA
perform_efa <- function(data, target_loadings, iter) {
  tryCatch({
    fa(r = data, nfactors = ncol(target_loadings), rotate = "oblimin", fm = "minres")
  }, error = function(e) {
    print(paste("EFA failed at iteration", iter, ":", e$message))
    return(NULL)
  })
}

# Function to perform Procrustes rotation
perform_procrustes <- function(f, target_loadings, iter) {
  # Perform Procrustes rotation on the factor loadings
  procrustes_result <- tryCatch({
    EFA.dimensions::PROCRUSTES(loadings = f, target = target_loadings, type = 'oblique', verbose = TRUE)
  }, error = function(e) {
    print(paste("ERROR in Procrustes rotation at iteration ", iter, ": ", e$message))
    return(NULL)
  }, warning = function(w) {
    print(paste("WARNING in Procrustes rotation at iteration ", iter, ": ", w$message))
    return(NULL)
  })
  
  # If Procrustes rotation was successful, proceed with the analysis
  if(!is.null(procrustes_result)){
    # Print the factor loadings after Procrustes rotation
    print(paste("Factor loadings after Procrustes rotation at iteration ", iter, ":"))
    print(procrustes_result$loadingsPROC)
    
    # Return the rotated loadings
    return(procrustes_result$loadingsPROC)
  } else {
    return(NULL)
  }
}

# Define a function to perform EFA and extract the suggested number of factors
efa_func <- function(data, target_loadings, n_iterations) {
  iter = 1
  repeat {
    # Print the iteration number
    print(paste("Starting iteration ", iter))
    
    # Resample the data with replacement
    resample <- data[sample(1:nrow(data), replace = TRUE), ]
    
    # Run EFA on the resampled data and handle any errors that may occur
    efa <- perform_efa(resample, target_loadings, iter)
    
    # If EFA was successful, proceed with the analysis
    if(!is.null(efa)){
      # Perform Procrustes rotation on the factor loadings
      efa$loadings <- perform_procrustes(efa$loadings, target_loadings, iter)
      
      # If Procrustes rotation was successful, proceed with the analysis
      if(!is.null(efa$loadings)){
        # Convert the factor loadings to a data frame and return the result
        loadings_df <- as.data.frame(efa$loadings)
        rownames(loadings_df) <- rownames(efa$loadings)
        colnames(loadings_df) <- colnames(efa$loadings)
        return(loadings_df)
      }
    }
    # Increment the iteration counter
    iter <- iter + 1
    
    # If the maximum number of iterations has been reached, print a message and break the loop
    if(iter > n_iterations){
      print(paste("EFA finished at iteration", iter - 1))
      break
    }
  }
  # If the function reaches this point without returning a result, return NULL
  return(NULL)
}

# Define a function to bootstrap EFA
bootstrap_efa <- function(data, n_iterations, target_loadings) {
  results_list <- vector("list", n_iterations)
  for (i in 1:n_iterations) {
    print(paste("Bootstrap EFA iteration ", i))
    efa_results <- efa_func(data, target_loadings, n_iterations)
    if(!is.null(efa_results)){
      results_list[[i]] <- efa_results
    }
  }
  # Filter out NULL values from the list
  results_list <- Filter(Negate(is.null), results_list)
  return(results_list)
}

# Bootstrap EFA with Procrustes rotation to the established structure
inpain_efa <- bootstrap_efa(inpain.df, n_iterations, efa_full_inpain$loadings)
painfree_efa <- bootstrap_efa(painfree.df, n_iterations, efa_full_painfree$loadings)

saveRDS(inpain_efa, "inpain_efa.rds")
saveRDS(painfree_efa, "painfree_efa.rds")

# Analyze loadings and create heatmaps
analyze_loadings <- function(loadings_list, group_name) {
  loadings_df <- do.call(rbind, lapply(loadings_list, function(x) {
    if(!is.null(x)){
      loadings = x
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

saveRDS(inpain.df, "inpain_df.rds")
saveRDS(painfree.df, "painfree_df.rds")
