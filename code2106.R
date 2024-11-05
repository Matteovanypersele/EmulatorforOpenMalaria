library(ggplot2)
library(gridExtra)
library(reshape2)
library(hetGP)  
library(sensitivity)

# Charger la base de données
data2 <- readRDS("C:/Users/vanyma/Documents/emulator_training01_20240618 (1).rds") 
#data2 <- readRDS("/scicore/home/pothin/GROUP/OM_emulator/emulator_training01_20240618 (1).rds")
meta <- data2$metadata
data <- data2$trainingData
data <- subset(data, Region == 'GreaterAccra')

# Vérifier la structure des données pour s'assurer qu'elles ont été correctement chargées
str(data)

# Sélectionner manuellement les colonnes X et Y
X_cols <- c("EIR", "histITNcov2018", "histITNcov2019", "histITNcov2020")
Y_cols <- c("prevalenceRate2020", "prevalenceRate2021", "prevalenceRate2022")

# Diviser les données en ensemble d'entraînement et de test
set.seed(123)
n <- nrow(data)
train_idx <- sample(seq_len(n), size = 0.01 * n)
train_data <- data[train_idx, ]
test_data <- data[-train_idx, ]

#############################ENTRAINEMENT HETGP####################################
train_hetGP <- function(X_cols, Y_cols, train_data) {
  if (!all(c(X_cols, Y_cols) %in% colnames(train_data))) {
    stop("Certaines des colonnes spécifiées n'existent pas dans les données.")
  }
  nvar <- length(X_cols)
  
  models <- lapply(Y_cols, function(Y_col) {
    start_time <- Sys.time()
    hetgp_model <- mleHetGP(X = as.matrix(train_data[, X_cols]), 
                            Z = train_data[[Y_col]], 
                            lower = rep(0, nvar), 
                            upper = rep(1, nvar), 
                            covtype = "Matern5_2")
    end_time <- Sys.time()
    training_time <- end_time - start_time
    cat("Temps d'entraînement du modèle HetGP pour", Y_col, ":", training_time, "secondes\n")
    return(hetgp_model)
  })
  names(models) <- Y_cols
  return(models)
}

# Entraîner le modèle HetGP
hetgp_models <- train_hetGP(X_cols, Y_cols, train_data)
inversematricecov <- hetgp_models$prevalenceRate2020$Ki

summary(hetgp_models$prevalenceRate2020)

#############################PREDICTIONS HETGP####################################
predict_hetGP <- function(hetgp_models, X_cols, Y_cols, test_data) {
  if (!all(c(X_cols, Y_cols) %in% colnames(test_data))) {
    stop("Certaines des colonnes spécifiées n'existent pas dans les données.")
  }
  
  results_list <- lapply(Y_cols, function(Y_col) {
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
      SD2 = hetgp_predictions$sd2  # Ajouter SD2 aux résultats
    )
    
    rmse_hetgp <- sqrt(mean((results$Actual - results$Predicted_HetGP)^2))
    mae_hetgp <- mean(abs(results$Actual - results$Predicted_HetGP))
    
    return(list(results = results, rmse = rmse_hetgp, mae = mae_hetgp))
  })
  names(results_list) <- Y_cols
  return(results_list)
}

# Effectuer les prédictions avec le modèle HetGP
hetgp_results_list <- predict_hetGP(hetgp_models, X_cols, Y_cols, test_data)

# Nouvelle fonction pour calculer les métriques de performance
calculate_loss <- function(hetgp_results_list) {
  loss_metrics <- do.call(rbind, lapply(names(hetgp_results_list), function(Y_col) {
    metrics <- hetgp_results_list[[Y_col]]
    rmse_hetgp <- metrics$rmse
    mae_hetgp <- metrics$mae
    
    return(data.frame(Y_col = Y_col, RMSE = rmse_hetgp, MAE = mae_hetgp))
  }))
  return(loss_metrics)
}

# Calculer les métriques de performance
loss_metrics <- calculate_loss(hetgp_results_list)

########################################RESULTATS########################################
for (Y_col in Y_cols) {
  cat("Résultats pour", Y_col, ":\n")
  print(hetgp_results_list[[Y_col]]$results)
}

cat("Metrics de performance:\n")
print(loss_metrics)

########################################COVARIANCE X,YPredicted########################################
calculate_covariance <- function(hetgp_results_list, test_data, X_cols) {
  covariance_metrics <- do.call(rbind, lapply(names(hetgp_results_list), function(Y_col) {
    results <- hetgp_results_list[[Y_col]]$results
    cov_vals <- sapply(X_cols, function(X_col) {
      cov(test_data[[X_col]], results$Predicted_HetGP)
    })
    cov_data <- data.frame(t(cov_vals))
    colnames(cov_data) <- X_cols
    cov_data$Y_col <- Y_col
    return(cov_data)
  }))
  return(covariance_metrics)
}

# Calculer les covariances
covariance_metrics <- calculate_covariance(hetgp_results_list, test_data, X_cols)

# Convertir les résultats des covariances en format long pour le graphique de heatmap
covariance_long <- melt(covariance_metrics, id.vars = 'Y_col')
pdf("Covariance.pdf", width = 8, height = 10)
# Tracer la heatmap des covariances
heatmap_plot <- ggplot(covariance_long, aes(x = variable, y = Y_col, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(min(covariance_long$value), max(covariance_long$value)), 
                       space = "Lab", name="Covariance") +
  theme_minimal() +
  labs(title = "Heatmap of Covariances between X Variables and Predicted Y Values",
       x = "X Variables",
       y = "Y Variables")


# Fermer l'appareil graphique PDF
dev.off()

# Exemple de fonction pour calculer l'objectif (remplacez par votre propre fonction)
objective_function <- function(X_cols) {
  # Simuler votre modèle avec les paramètres spécifiés
  # Ici, nous supposons que 'params' est un vecteur contenant les valeurs des paramètres
  # Remplacez cette partie par votre propre calcul de l'objectif
  predicted <- Predicted_HetGP
  # Calculer l'objectif (par exemple, la somme des carrés des résidus)
  return(sum((test_data$Observed - predicted)^2))
}

# Définir les limites des paramètres (remplacez par vos propres limites)
param_bounds <- data.frame(
  param1 = c(min = 0, max = 1),
  param2 = c(min = 0, max = 1),
  param3 = c(min = 0, max = 1)
)

# Générer des échantillons pour l'analyse de Sobol
# Utiliser la fonction 'sobol2002' pour générer les échantillons et calculer les indices de Sobol
n <- 1000  # Nombre d'échantillons (ajustez selon vos besoins)
X1 <- data.frame(lapply(param_bounds, function(x) runif(n, x[1], x[2])))
X2 <- data.frame(lapply(param_bounds, function(x) runif(n, x[1], x[2])))

sobol_results <- sobol2002(model = objective_function, X1 = X1, X2 = X2, nboot = 100)

# Résultats des indices de Sobol
print(sobol_results)

# Convertir les résultats pour la visualisation
sobol_long <- melt(sobol_results$S, id.vars = 'variable')

# Créer un PDF avec la heatmap des indices de Sobol
pdf("Sobol_Analysis.pdf", width = 8, height = 10)

# Tracer la heatmap des indices de Sobol
sobol_plot <- ggplot(sobol_long, aes(x = variable, y = index, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(min(sobol_long$value), max(sobol_long$value)), 
                       space = "Lab", name="Sobol Index") +
  theme_minimal() +
  labs(title = "Heatmap of Sobol Indices",
       x = "Parameters",
       y = "Objective")

print(sobol_plot)  # Imprimer le graphique dans le PDF

# Fermer l'appareil graphique PDF
dev.off()




pdf("ActualvsPredicted.pdf", width = 8, height = 10)
########################################ACTUAL VS PREDICTED########################################
plots <- lapply(Y_cols, function(Y_col) {
  results <- hetgp_results_list[[Y_col]]$results
  p <- ggplot(results, aes(x = Actual, y = Predicted_HetGP)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    labs(title = paste("Actual vs Predicted with Heteroscedastic Gaussian Process Model for", Y_col),
         x = paste("Actual", Y_col),
         y = paste("Predicted", Y_col)) +
    theme_minimal()
  return(p)
})

# Afficher les graphiques ensemble
do.call(grid.arrange, c(plots, ncol = 1))
# Fermer l'appareil graphique PDF
dev.off()
#############################ANALYSE DE SENSIBILITÉ####################################
sensitivity_analysis <- function(hetgp_model, X_cols, test_data, var_to_vary) {
  X_mean <- colMeans(test_data[, X_cols])
  sensitivity_data <- data.frame()
  
  for (value in seq(min(test_data[[var_to_vary]]), max(test_data[[var_to_vary]]), length.out = 100)) {
    X_mean[var_to_vary] <- value
    hetgp_predictions <- predict(hetgp_model, x = as.matrix(t(X_mean)))
    sensitivity_data <- rbind(sensitivity_data, data.frame(Value = value, Predicted = hetgp_predictions$mean))
  }
  
  return(sensitivity_data)
}

# Tracer la variation de Y en fonction de la variation d'une variable X
plot_sensitivity <- function(hetgp_models, X_cols, test_data, var_to_vary) {
  sensitivity_plots <- lapply(names(hetgp_models), function(Y_col) {
    hetgp_model <- hetgp_models[[Y_col]]
    sensitivity_data <- sensitivity_analysis(hetgp_model, X_cols, test_data, var_to_vary)
    p <- ggplot(sensitivity_data, aes(x = Value, y = Predicted)) +
      geom_line() +
      labs(title = paste("Sensitivity Analysis for", Y_col, "with varying", var_to_vary),
           x = var_to_vary,
           y = paste("Predicted", Y_col)) +
      theme_minimal()
    return(p)
  })
  return(sensitivity_plots)
}
pdf("Sensibility.pdf", width = 8, height = 10)
# Choisissez une variable à faire varier, par exemple "EIR"
var_to_vary <- "histITNcov2020"
sensitivity_plots1 <- plot_sensitivity(hetgp_models, X_cols, test_data, var_to_vary)
# Afficher les graphiques de sensibilité ensemble
do.call(grid.arrange, c(sensitivity_plots1, ncol = 4))
# Fermer l'appareil graphique PDF
dev.off()
# Fonction pour tracer l'analyse de sensibilité
plot_sensitivity <- function(sensitivity_data, var_to_vary, y_col) {
  ggplot(sensitivity_data, aes(x = Value, y = Predicted)) +
    geom_line(color = 'blue', size = 1) +
    labs(title = paste("Sensitivity Analysis for", var_to_vary, "on", y_col),
         x = var_to_vary,
         y = "Predicted Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
}

# Variables à analyser
vars_to_vary <- c("EIR", "histITNcov2018", "histITNcov2019", "histITNcov2020")

# Ouvrir un appareil graphique PDF
pdf("Sensitivity.pdf", width = 8, height = 10)

# Générer et tracer les graphiques de sensibilité pour chaque variable à analyser
for (var_to_vary in vars_to_vary) {
  sensitivity_plots <- lapply(Y_cols, function(Y_col) {
    hetgp_model <- hetgp_models[[Y_col]]
    sensitivity_data <- sensitivity_analysis(hetgp_model, X_cols, test_data, var_to_vary)
    plot_sensitivity(sensitivity_data, var_to_vary, Y_col)
  })
  
  # Afficher les graphiques de sensibilité ensemble sur une nouvelle page du PDF
  do.call(grid.arrange, c(sensitivity_plots, ncol = 2))
}

# Fermer l'appareil graphique PDF
dev.off()

