library(dplyr)
library(ggplot2)
library(hetGP)
library(roxygen2)


#' Predict using HetGP models
#'
#' This function generates predictions using trained HetGP models for specified target variables.
#'
#' @param hetgp_models A named list of trained HetGP models, with names corresponding to the target variables (Y_cols).
#' @param X_cols A character vector of predictor variable names (features).
#' @param Y_cols A character vector of target variable names (responses).
#' @param test_data A data frame containing the test data, with columns matching X_cols and Y_cols.
#'
#' @return A named list of data frames, each containing the actual values, predicted values, and prediction intervals for a target variable.
#' @examples
#' \dontrun{
#'predictionsnewdata <- predict_hetGP(hetgp_models, X_cols, Y_cols, test_data)
#' }
predict_hetGP <- function(hetgp_models, X_cols, Y_cols, test_data) {
  if (!all(c(X_cols, Y_cols) %in% colnames(test_data))) {
    stop("Some of the specified columns do not exist in the data.")
  }
  
  predictions_list <- lapply(Y_cols, function(Y_col) {
    hetgp_model <- hetgp_models[[Y_col]]
    hetgp_predictions <- predict(hetgp_model, x = as.matrix(test_data[, X_cols]))
    predicted_hetgp_values <- hetgp_predictions$mean
    lower_hetgp_bounds <- predicted_hetgp_values - 1.96 * sqrt(hetgp_predictions$sd2)
    upper_hetgp_bounds <- predicted_hetgp_values + 1.96 * sqrt(hetgp_predictions$sd2)
    
    results <- data.frame(
      Actual = test_data[[Y_col]],
      Predicted_HetGP = predicted_hetgp_values,
      Lower_HetGP = lower_hetgp_bounds,
      Upper_HetGP = upper_hetgp_bounds,
      SD2 = hetgp_predictions$sd2
    )
    
    return(results)
  })
  names(predictions_list) <- Y_cols
  return(predictions_list)
}

#' Evaluate HetGP model predictions
#'
#' This function evaluates the predictions made by HetGP models, calculating performance metrics such as RMSE, MAE, MAPE, and R-squared.
#'
#' @param predictions_list A named list of data frames, each containing actual values, predicted values, and prediction intervals for a target variable.
#'
#' @return A named list of evaluation metrics for each target variable, including RMSE, MAE, MAPE, and R-squared.
#' @examples
#' \dontrun{
#'evaluationnewdata <- evaluate_hetGP(predictionsnewdata)
#' }
evaluate_hetGP <- function(predictions_list) {
  evaluation_list <- lapply(predictions_list, function(results) {
    rmse_hetgp <- sqrt(mean((results$Actual - results$Predicted_HetGP)^2))
    mae_hetgp <- mean(abs(results$Actual - results$Predicted_HetGP))
    
    # Exclude rows where Actual is zero for MAPE calculation
    mape_data <- subset(results, Actual != 0)
    mape_hetgp <- mean(abs((mape_data$Actual - mape_data$Predicted_HetGP) / mape_data$Actual)) * 100
    
    results$MAPE <- abs((results$Actual - results$Predicted_HetGP) / results$Actual) * 100
    results$RMSE <- abs(results$Actual - results$Predicted_HetGP)
    ss_tot <- sum((results$Actual - mean(results$Actual))^2)
    ss_res <- sum((results$Actual - results$Predicted_HetGP)^2)
    rsquared_hetgp <- 1 - (ss_res / ss_tot)
    
    return(list(results = results, rmse = rmse_hetgp, mae = mae_hetgp, mape = mape_hetgp, rsquared = rsquared_hetgp))
  })
  names(evaluation_list) <- names(predictions_list)
  return(evaluation_list)
}

#' Make predictions across multiple folds using HetGP models
#'
#' This function generates predictions for multiple cross-validation folds using HetGP models.
#'
#' @param all_models A list of trained HetGP models for each fold.
#' @param test_data_folds A list of data frames containing the test data for each fold.
#' @param X_cols A character vector of predictor variable names (features).
#' @param Y_cols A character vector of target variable names (responses).
#' @param num_folds An integer specifying the number of folds. Default is 10.
#'
#' @return A list of prediction results for each fold.
#' @examples
#' \dontrun{
#' all_predictions <- predict_models(all_models, test_data_folds, X_cols, Y_cols, num_folds = 10)
#' }
predict_models <- function(all_models, test_data_folds, X_cols, Y_cols, num_folds = 10) {
  all_predictions <- list()
  
  for (i in 1:num_folds) {
    cat("Predicting for fold:", i, "\n")
    
    # Make predictions with the HetGP models
    predictions <- predict_hetGP(all_models[[paste0("fold_", i)]], X_cols, Y_cols, test_data_folds[[paste0("test_data_fold_", i)]])
    
    # Save the predictions for this fold
    all_predictions[[paste0("fold_", i)]] <- predictions
  }
  
  return(all_predictions)
}

#' Evaluate predictions across multiple folds
#'
#' This function evaluates predictions made across multiple cross-validation folds using HetGP models.
#'
#' @param all_predictions A list of prediction results for each fold.
#' @param num_folds An integer specifying the number of folds. Default is 10.
#'
#' @return A list of evaluation metrics for each fold.
#' @examples
#' \dontrun{
#' all_evaluations <- evaluate_predictions(all_predictions, num_folds = 10)
#' }
evaluate_predictions <- function(all_predictions, num_folds = 10) {
  all_evaluations <- list()
  
  for (i in 1:num_folds) {
    cat("Evaluating fold:", i, "\n")
    
    # Evaluate the predictions
    evaluations <- evaluate_hetGP(all_predictions[[paste0("fold_", i)]])
    
    # Save the evaluation results for this fold
    all_evaluations[[paste0("fold_", i)]] <- evaluations
  }
  
  return(all_evaluations)
}




# Conversion des résultats d'évaluation en data frame pour ggplot
evaluation_df <- do.call(rbind, lapply(names(evaluationnewdata), function(name) {
  eval_result <- evaluationnewdata[[name]]
  data.frame(
    Y_col = name,
    R_squared = eval_result$rsquared,
    RMSE = eval_result$rmse,
    MAE = eval_result$mae,
    MAPE = eval_result$mape
  )
}))
write.csv(evaluation_df, "evaluationnewdata.csv", row.names = FALSE)


# Graphique pour R_squared
ggplot(evaluation_df, aes(x = Y_col, y = R_squared)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = expression(R^2~"Values for Each Response Variable"), x = "Y Columns", y = expression(R^2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Graphique pour MAPE
ggplot(evaluation_df, aes(x = Y_col, y = MAPE)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = expression(MAPE~"Values for Each Response Variable"), x = "Y Columns", y = "MAPE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Graphique pour RMSE
ggplot(evaluation_df, aes(x = Y_col, y = RMSE)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = expression(RMSE~"Values for EIR8ageALL"), x = "Y Columns", y = "RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


results <- evaluationnewdata$prevalenceRate_2024$results
table(results$Predicted_HetGP)
p <- ggplot(results, aes(x = Actual, y = Predicted_HetGP)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c() +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(title = paste("Actual vs Predicted with Heteroscedastic Gaussian Process Model for prevalencerate2024 newdata"),
       x = paste("Actual prevalencerate2024"),
       y = paste("Predicted prevalencerate2024")) +
  theme_minimal()
print(p)
# Extract the results for the specific Y_col, e.g., nTreatments1_2027
p <- ggplot(results, aes(x = Actual, y = MAPE)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c() +
  labs(title = "MAPE vs Actual prevalencerate 2027 for EIR8Age05",
       x = "Actual",
       y = "MAPE (%)") +
  theme_minimal()
print(p)

# Create the plot
p <- ggplot(results, aes(x = Actual, y = RMSE)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c() +
  labs(title = "AbsoluteError vs Actual prevalencerate 2027 for EIR8Age05",
       x = "Actual",
       y = "AE") +
  theme_minimal()

# Display the plot
print(p)

