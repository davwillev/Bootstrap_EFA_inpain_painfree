# Load necessary packages
packages <- c('tidyverse', 
              'psych')
lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})

# Set seed for reproducibility
set.seed(123)

# Define the number of bootstrap iterations
n_iterations <- 10 # Very slow but change to 1000 for actual analysis

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

# Determine the number of factors to retain using MAP criterion
bootstrap_map <- function(data, n_iterations) {
  # Determine the maximum number of factors
  n_factors <- min(nrow(data), ncol(data))
  
  # Initialize an array to hold the number of factors suggested in each bootstrap sample
  suggested_factors <- vector("numeric", n_iterations)
  
  # Iterate over bootstrap samples
  for (i in 1:n_iterations) {  
    # Resample data with replacement
    indices <- sample(1:nrow(data), replace = TRUE)
    resample <- data[indices, ]
    
    # Initialize an array to hold MAP values for each number of factors
    map_values <- rep(Inf, n_factors)
    
    # Compute MAP for each number of factors
    for (j in 1:n_factors) {
      # Compute the Velicer's MAP criterion for 'j' factors
      vss <- try(VSS(resample, n = j, fm = "minres", plot = FALSE), silent = TRUE)
      
      # Check if VSS computation was successful
      if(!inherits(vss, "try-error")) {
        map_values[j] <- vss$map
      }
    }
    
    # The optimal number of factors is the one that minimizes the MAP value
    suggested_factors[i] <- which.min(map_values)
    
    # Save intermediate output
    saveRDS(list(map_values = map_values, suggested_factors = suggested_factors),
            file = paste0("bootstrap_", i, "_map.rds"))
    
    # Print progress information
    print(paste("Iteration", i, ": Retaining", suggested_factors[i], "factors"))
  }
  
  # The function returns the number of factors suggested in each bootstrap sample
  return(suggested_factors)
}

# Run bootstrap MAP analysis
inpain_suggested_factors <- bootstrap_map(inpain.df, n_iterations)
painfree_suggested_factors <- bootstrap_map(painfree.df, n_iterations)

# Save final output
saveRDS(list(inpain_suggested_factors = inpain_suggested_factors,
             painfree_suggested_factors = painfree_suggested_factors),
        file = "final_bootstrap_map.rds")

# Create a combined dataframe to compare datasets
combined_df <- rbind(
  data.frame(Group = "Inpain", SuggestedFactors = inpain_suggested_factors),
  data.frame(Group = "Painfree", SuggestedFactors = painfree_suggested_factors)
)

# Plot histograms of number of factors side by side
ggplot(combined_df, aes(x = SuggestedFactors)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(breaks = seq(floor(min(combined_df$SuggestedFactors)), ceiling(max(combined_df$SuggestedFactors)), by = 1)) +
  facet_grid(. ~ Group) +
  labs(title = "Suggested Number of Factors", x = "Number of Factors", y = "Count")

# Function to decide on the number of factors to retain
decide_on_factors <- function(suggested_factors) {
  factor_frequencies <- table(suggested_factors)
  sorted_frequencies <- sort(factor_frequencies, decreasing = TRUE)
  if (length(sorted_frequencies) > 1 && sorted_frequencies[1] >= 1.2 * sorted_frequencies[2]) {
    return(as.integer(names(sorted_frequencies)[1]))
  } else {
    return(as.integer(round(median(suggested_factors))))
  }
}

# Decide on the number of factors to retain
n_factors_inpain <- decide_on_factors(inpain_suggested_factors)
n_factors_painfree <- decide_on_factors(painfree_suggested_factors)

# Save important variables into RDS file
saveRDS(list(n_factors_inpain = n_factors_inpain,
             n_factors_painfree = n_factors_painfree,
             inpain_suggested_factors = inpain_suggested_factors,
             painfree_suggested_factors = painfree_suggested_factors,
             inpain_df = inpain.df, 
             painfree_df = painfree.df),
        file = "determine_factors.rds")
