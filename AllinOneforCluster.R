# Load necessary libraries
library(hetGP)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(sensitivity)
library(purrr)

# Accept the dataset index as an argument
args <- commandArgs(trailingOnly = TRUE)
dataset_index <- as.integer(args[1])

# List of datasets
datasets <- c(
  "dataset1_EIR8_age_05.rds",
  "dataset2_EIR8_age_all.rds",
  "dataset3_EIR22_age_05.rds",
  "dataset4_EIR22_age_all.rds",
  "dataset5_EIR88_age_05.rds",
  "dataset6_EIR88_age_all.rds",
  "dataset7_EIR168_age_05.rds",
  "dataset8_EIR168_age_all.rds"
)

# Build the path to the dataset
dataset_path <- paste0("/scicore/home/pothin/GROUP/OM_emulator/datasets_22_07_vary_fut_cov/", datasets[dataset_index])

# Load the data
load(dataset_path)

# Extract the dataset name (without extension)
dataset_name <- sub(".rds$", "", datasets[dataset_index])

# Assign the dataset to a variable
data <- get(dataset_name)


data <- subset(data, year > 2022 & year < 2028)
data_wide <- data %>%
  select(nTreatments1, year, prevalenceRate, ID, seed, futITNcov2023, futITNcov2024, futITNcov2025, futITNcov2026, futITNcov2027, futIRScov2023, futIRScov2024, futIRScov2025, futIRScov2026, futIRScov2027, access2023, access2024, access2025, access2026, access2027) %>%
  pivot_wider(names_from = year, values_from = c(nTreatments1, prevalenceRate), names_glue = "{.value}_{year}")

# Manually select the X and Y columns
X_cols <- c("futITNcov2023", "futITNcov2024", "futITNcov2025", "futITNcov2026", "futITNcov2027", 
            "futIRScov2023", "futIRScov2024", "futIRScov2025", "futIRScov2026", "futIRScov2027", 
            "access2023", "access2024", "access2025", "access2026", "access2027") 

Y_cols <- c("nTreatments1_2023", "nTreatments1_2024", "nTreatments1_2025", "nTreatments1_2026", 
            "nTreatments1_2027", "prevalenceRate_2023", "prevalenceRate_2024", "prevalenceRate_2025",
            "prevalenceRate_2026", "prevalenceRate_2027")
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
#fonction to predict
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
#fonction to evaluate
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

# Function to sample by group
sampleByGroup <- function(df, sampleFraction = 0.1, split_group = "seed", s = 1) {
  df %>%
    mutate(idx = row_number()) %>%
    group_by_at(split_group) %>%
    group_split() %>%
    map(function(x) {
      set.seed(s)
      x[sample(1:nrow(x), floor(nrow(x) * sampleFraction)), ]
    }) %>%
    bind_rows()
}
# Function to create training and test datasets
create_folds <- function(data_wide, num_folds = 10, fraction = 0.1, split_group = "seed") {
  train_data_folds <- list()
  test_data_folds <- list()
  
  for (i in 1:num_folds) {
    set.seed(123 + i)
    train_data <- sampleByGroup(data_wide, sampleFraction = fraction, split_group = split_group, s = 123 + i)
    
    test_data <- data_wide %>%
      filter(!ID %in% train_data$ID)
    
    # Store the datasets for each fold
    train_data_folds[[paste0("train_data_fold_", i)]] <- train_data
    test_data_folds[[paste0("test_data_fold_", i)]] <- test_data
  }
  
  list(train_data_folds = train_data_folds, test_data_folds = test_data_folds)
}

# Function to train the models
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

# Function to make predictions with the trained models
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

# Function to evaluate the predictions
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

# Call the functions to create the datasets, train the models, make predictions, and evaluate
folds <- create_folds(data_wide, num_folds = 5, fraction = 0.2)
train_data_folds <- folds$train_data_folds
test_data_folds <- folds$test_data_folds

# Train the models
all_models <- train_models(train_data_folds, X_cols, Y_cols, num_folds = 5)

# Make predictions
all_predictions <- predict_models(all_models, test_data_folds, X_cols, Y_cols, num_folds = 5)

# Evaluate the predictions
all_evaluations <- evaluate_predictions(all_predictions, num_folds = 5)
# Créer un dossier de sortie pour chaque dataset
output_dir <- paste0("JobOut/dataset_", dataset_index)
dir.create(output_dir, recursive = TRUE)

# Sauvegarder les modèles entraînés
saveRDS(all_models, file = paste0(output_dir, "/trained_models.rds"))

# Sauvegarder les prédictions
saveRDS(all_predictions, file = paste0(output_dir, "/predictions.rds"))

# Sauvegarder les évaluations
saveRDS(all_evaluations, file = paste0(output_dir, "/evaluations.rds"))