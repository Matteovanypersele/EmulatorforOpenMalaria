library(hetGP)
library(roxygen2)

#' Train a HetGP model for multiple target variables
#'
#' This function trains a separate HetGP model for each specified target variable (Y_col) using the provided training data.
#'
#' @param X_cols A character vector of predictor variable names (features).
#' @param Y_cols A character vector of target variable names (responses).
#' @param train_data A data frame containing the training data, with columns matching X_cols and Y_cols.
#'
#' @return A named list of HetGP models, with names corresponding to the Y_cols.
#' @examples
#' \dontrun{
#' hetgp_model <- train_hetGP(X_cols, Y_cols, train_data)
#' }
 #fonction to train hetgp
 train_hetGP <- function(X_cols, Y_cols, train_data) {
   models <- lapply(Y_cols, function(Y_col) {
     start_time <- Sys.time()
     hetgp_model <- mleHetGP(X = as.matrix(train_data[, X_cols]), 
                             Z = train_data[[Y_col]], 
                             covtype = "Matern5_2")
     end_time <- Sys.time()
     training_time <- difftime(end_time, start_time, units = "secs")
     cat("Training time for HetGP model for", Y_col, ":", training_time, "seconds\n")
     return(hetgp_model)
   })
   
   names(models) <- Y_cols
   return(models)
 }

#' Train HetGP models across multiple folds of training data
#'
#' This function trains HetGP models for each fold in a cross-validation setup.
#'
#' @param train_data_folds A list of data frames, each containing the training data for one fold.
#' @param X_cols A character vector of predictor variable names (features).
#' @param Y_cols A character vector of target variable names (responses).
#' @param num_folds An integer specifying the number of folds. Default is 10.
#'
#' @return A list of trained HetGP models for each fold.
#' @examples
#' \dontrun{
#' all_models <- train_models(train_data_folds, X_cols, Y_cols, num_folds = 10)
#' }
train_models <- function(train_data_folds, X_cols, Y_cols, num_folds = 10) {
  all_models <- list()
  
  for (i in 1:num_folds) {
    cat("Training fold:", i, "\n")
    
    # Train the HetGP models
    hetgp_models <- train_hetGP(X_cols, Y_cols, train_data_folds[[paste0("train_data_fold_", i)]])
    
    # Save the models for this fold
    all_models[[paste0("fold_", i)]] <- hetgp_models
  }
  
  return(all_models)
}
