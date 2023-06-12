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
boot_iterations <- 1000

# Load data from RDS
iteration_number <- readRDS("iteration_number.rds")
determine_factors <- readRDS("determine_factors.rds")
inpain.df <- determine_factors$inpain_df
painfree.df <- determine_factors$painfree_df
n_factors_inpain <- determine_factors$n_factors_inpain
n_factors_painfree <- determine_factors$n_factors_painfree
inpain_suggested_factors <- determine_factors$inpain_suggested_factors
painfree_suggested_factors <- determine_factors$painfree_suggested_factors

# Run EFA with a fixed number of factors to establish the factor structure
efa_full_inpain <- fa(inpain.df, nfactors = n_factors_inpain, rotate = "oblimin", fm = "minres")
efa_full_painfree <- fa(painfree.df, nfactors = n_factors_painfree, rotate = "oblimin", fm = "minres")

# Print initial EFA
print(efa_full_inpain)
print(efa_full_painfree)

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
perform_procrustes <- function(f, target_loadings, iter, rotation_type = 'oblique', verbose = TRUE) {
  # Perform Procrustes rotation on the factor loadings
  procrustes_result <- tryCatch({
    EFA.dimensions::PROCRUSTES(loadings = f, target = target_loadings, type = rotation_type, verbose = verbose)
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
    if(verbose) {
      print(paste("Factor loadings after Procrustes rotation at iteration ", iter, ":"))
      print(procrustes_result$loadingsPROC)
    }
    # Return the full Procrustes result instead of just the rotated loadings
    return(procrustes_result)
  } else {
    return(NULL)
  }
}

# Define a function to perform EFA and extract the suggested number of factors
efa_func <- function(data, target_loadings, boot_iterations, rotation_type = 'oblique', verbose = TRUE) {
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
      procrustes_result <- perform_procrustes(efa$loadings, target_loadings, iter, rotation_type, verbose)
      
      # If Procrustes rotation was successful, proceed with the analysis
      if(!is.null(procrustes_result)){
        # Convert the factor loadings to a data frame and return the result
        loadings_df <- as.data.frame(procrustes_result$loadingsPROC)
        rownames(loadings_df) <- rownames(procrustes_result$loadingsPROC)
        colnames(loadings_df) <- colnames(procrustes_result$loadingsPROC)
        
        # Include additional results from Procrustes rotation in the output
        result_list <- list(
          'loadings_df' = loadings_df,
          'eigenvalues' = efa$values,
          'congruence' = procrustes_result$congruence,
          'rmsr' = procrustes_result$rmsr,
          'residmat' = procrustes_result$residmat
        )
        return(result_list)
      }
    }
    # Increment the iteration counter
    iter <- iter + 1
    
    # If the maximum number of iterations has been reached, print a message and break the loop
    if(iter > boot_iterations){
      print(paste("EFA finished at iteration", iter - 1))
      break
    }
  }
  # If the function reaches this point without returning a result, return NULL
  return(NULL)
}

# Define a function to bootstrap EFA
bootstrap_efa <- function(data, boot_iterations, target_loadings) {
  results_list <- vector("list", boot_iterations)
  for (i in 1:boot_iterations) {
    print(paste("Bootstrap EFA iteration ", i))
    efa_results <- efa_func(data, target_loadings, boot_iterations)
    if(!is.null(efa_results)){
      results_list[[i]] <- list(
        'loadings_df' = efa_results$loadings_df,
        'eigenvalues' = efa_results$eigenvalues,
        'congruence' = efa_results$congruence,
        'rmsr' = efa_results$rmsr,
        'residmat' = efa_results$residmat
      )
    }
  }
  # Filter out NULL values from the list
  results_list <- Filter(Negate(is.null), results_list)
  return(results_list)
}

# Bootstrap EFA with Procrustes rotation to the established structure
inpain_efa <- bootstrap_efa(inpain.df, boot_iterations, efa_full_inpain$loadings)
painfree_efa <- bootstrap_efa(painfree.df, boot_iterations, efa_full_painfree$loadings)

# Save each dataframe to an RDS file
saveRDS(inpain_efa, "inpain_efa.rds")
saveRDS(painfree_efa, "painfree_efa.rds")
