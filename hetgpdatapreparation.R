# Load necessary libraries
library(hetGP)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(sensitivity)
library(purrr)

# Accept the dataset index as an argument (for array usage on the cluster)
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

# Prepare the data
data <- dataset1_EIR8_age_05
data <- subset(data, year > 2022 & year < 2028)
#select the columns needed

data_wide <- data %>%
  select(nTreatments1, year, prevalenceRate, ID, seed, futITNcov2023, futITNcov2024, futITNcov2025, futITNcov2026, futITNcov2027, futIRScov2023, futIRScov2024, futIRScov2025, futIRScov2026, futIRScov2027, access2023, access2024, access2025, access2026, access2027) %>%
  pivot_wider(names_from = year, values_from = c(nTreatments1, prevalenceRate), names_glue = "{.value}_{year}")
data <- data_wide
# Manually select the X and Y columns
X_cols <- c("futITNcov2023", "futITNcov2024", "futITNcov2025", "futITNcov2026", "futITNcov2027", 
            "futIRScov2023", "futIRScov2024", "futIRScov2025", "futIRScov2026", "futIRScov2027", 
            "access2023", "access2024", "access2025", "access2026", "access2027") 

Y_cols <- c("nTreatments1_2023", "nTreatments1_2024", "nTreatments1_2025", "nTreatments1_2026", 
            "nTreatments1_2027", "prevalenceRate_2023", "prevalenceRate_2024", "prevalenceRate_2025",
            "prevalenceRate_2026", "prevalenceRate_2027")
#Check that those columns have a numeric or integer type and that they have some variation
# Check data types of X columns
x_data_types <- sapply(data[, X_cols], class)

# Check data types of Y columns
y_data_types <- sapply(data[, Y_cols], class)

# Print the results
print("X Columns Data Types:")
print(x_data_types)
print("Y Columns Data Types:")
print(y_data_types)
# Check variation in X columns
x_variation <- sapply(train_data[, X_cols], function(col) length(unique(col)))

# Check variation in Y columns
y_variation <- sapply(train_data[, Y_cols], function(col) length(unique(col)))

# Print the results
print("X Columns Variation (Number of Unique Values):")
print(x_variation)
print("Y Columns Variation (Number of Unique Values):")
print(y_variation)

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
# Split the data into training and test sets
set.seed(123)
n <- nrow(data)
train_idx <- sample(seq_len(n), size = 0.1 * n)
train_data <- data[train_idx, ]
test_data <- data[-train_idx, ]

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