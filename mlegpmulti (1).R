library(mlegp)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Load the dataset
data2 <- readRDS("/scicore/home/pothin/GROUP/OM_emulator/emulator_training01_20240618 (1).rds")
meta <- data2$metadata
data <- data2$trainingData
data <- subset(data, Region == 'GreaterAccra')

# Verify the structure of the data to ensure it has been loaded correctly
str(data)

# Manually select X and Y columns
X_cols <- c("EIR", "histITNcov2018", "histITNcov2019", "histITNcov2020")
Y_cols <- c("prevalenceRate2020", "prevalenceRate2021", "prevalenceRate2022")

# Split the data into training and test sets
set.seed(123)
n <- nrow(data)
train_idx <- sample(seq_len(n), size = 0.01 * n)
train_data <- data[train_idx, ]
test_data <- data[-train_idx, ]

#' Train mlegp models with principal component weights
#'
#' This function trains Gaussian process models using principal component weights.
#'
#' @title Train mlegp models
#' @description This function trains Gaussian process models using principal component weights.
#' @param X_cols A vector of column names for the predictors.
#' @param Y_cols A vector of column names for the response variables.
#' @param train_data The training data set.
#' @param numPCs The number of principal components to use. Default is 2.
#' @return An mlegp model fitted with principal component weights.
#' @export
train_mlegp_pc <- function(X_cols, Y_cols, train_data, numPCs = 2) {
  if (!all(c(X_cols, Y_cols) %in% colnames(train_data))) {
    stop("Some of the specified columns do not exist in the data.")
  }
  
  X <- as.matrix(train_data[, X_cols])
  Y <- as.matrix(train_data[, Y_cols])
  
  # Fit GP models to principal component weights
  start_time <- Sys.time()
  mlegp_model <- mlegp(X, t(Y), PC.num = numPCs)
  end_time <- Sys.time()
  training_time <- end_time - start_time
  cat("Training time for mlegp model with principal component weights:", training_time, "seconds\n")
  
  return(mlegp_model)
}

# Train mlegp models
mlegp_model_pc <- train_mlegp_pc(X_cols, Y_cols, train_data)

#' Predict using the fitted mlegp model
#'
#' This function makes predictions using the fitted mlegp model.
#'
#' @title Predict with mlegp model
#' @description This function makes predictions using the fitted mlegp model.
#' @param model The fitted mlegp model.
#' @param new_data The new data set for prediction.
#' @param X_cols A vector of column names for the predictors.
#' @param numPCs The number of principal components used in the model. Default is 2.
#' @return Predicted values.
#' @export
predict_mlegp_pc <- function(model, new_data, X_cols, numPCs = 2) {
  if (!all(X_cols %in% colnames(new_data))) {
    stop("Some of the specified columns do not exist in the data.")
  }
  
  X_new <- as.matrix(new_data[, X_cols])
  
  # Make predictions for each principal component
  Vprime_pred <- matrix(0, numPCs, nrow(X_new))
  for (i in 1:numPCs) {
    Vprime_pred[i, ] <- predict(model[[i]], X_new)
  }
  
  # Reconstruct the original output predictions
  UD <- model$UD
  
  if (ncol(UD) != nrow(Vprime_pred)) {
    stop("Dimensions of UD and Vprime_pred do not conform for multiplication.")
  }
  
  predY <- UD %*% Vprime_pred
  
  return(predY)
}

# Make predictions with the test data
preds <- predict_mlegp_pc(mlegp_model_pc, test_data, X_cols)

#' Calculate Mean Squared Error (MSE)
#'
#' This function calculates the Mean Squared Error (MSE) between true and predicted values.
#'
#' @title Calculate MSE
#' @description This function calculates the Mean Squared Error (MSE) between true and predicted values.
#' @param true The true values.
#' @param predicted The predicted values.
#' @return The MSE value.
#' @export
calc_mse <- function(true, predicted) {
  return(mean((true - predicted)^2))
}

# Calculate MSE for each prevalence year
mse_2020 <- calc_mse(test_data$prevalenceRate2020, preds[1, ])
mse_2021 <- calc_mse(test_data$prevalenceRate2021, preds[2, ])
mse_2022 <- calc_mse(test_data$prevalenceRate2022, preds[3, ])

cat("MSE for 2020:", mse_2020, "\n")
cat("MSE for 2021:", mse_2021, "\n")
cat("MSE for 2022:", mse_2022, "\n")

#' Calculate R-squared (R²)
#'
#' This function calculates the R-squared (R²) between true and predicted values.
#'
#' @title Calculate R-squared (R²)
#' @description This function calculates the R-squared (R²) between true and predicted values.
#' @param true The true values.
#' @param predicted The predicted values.
#' @return The R² value.
#' @export
calc_r2 <- function(true, predicted) {
  ss_res <- sum((true - predicted)^2)
  ss_tot <- sum((true - mean(true))^2)
  return(1 - ss_res / ss_tot)
}

# Calculate R² for each prevalence year
r2_2020 <- calc_r2(test_data$prevalenceRate2020, preds[1, ])
r2_2021 <- calc_r2(test_data$prevalenceRate2021, preds[2, ])
r2_2022 <- calc_r2(test_data$prevalenceRate2022, preds[3, ])

cat("R² for 2020:", r2_2020, "\n")
cat("R² for 2021:", r2_2021, "\n")
cat("R² for 2022:", r2_2022, "\n")

# Visualize predictions vs actual values
df_plot <- data.frame(
  Year = rep(c(2020, 2021, 2022), each = nrow(test_data)),
  Actual = c(test_data$prevalenceRate2020, test_data$prevalenceRate2021, test_data$prevalenceRate2022),
  Predicted = c(preds[1, ], preds[2, ], preds[3, ])
)

# Open a PDF graphics device
pdf("ActVSPred.pdf", width = 8, height = 10)

ggplot(df_plot, aes(x = Actual, y = Predicted, color = as.factor(Year))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ Year) +
  labs(title = "Predictions vs Actual Values", x = "Actual Values", y = "Predictions") +
  theme_minimal()

dev.off()
