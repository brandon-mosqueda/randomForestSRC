tune.zap.rfsrc <- function(X, y,
                           ntree_theta = c(500),
                           mtry_theta = NULL,
                           nodesize_theta = c(5, 15),
                           ntree_lambda = c(500),
                           mtry_lambda = NULL,
                           nodesize_lambda = c(5, 15),
                           type = "original",
                           importance = FALSE,
                           loss_function = mse,

                           cross_validation = "k_fold",
                           number_of_folds = 5,
                           proportion_of_testing = 0.2,
                           sample_proportion = 1,
                           verbose = TRUE) {

  validate_tune_zap_params(X, y, ntree_theta, mtry_theta,
                           nodesize_theta, ntree_lambda,
                           mtry_lambda, nodesize_lambda,
                           loss_function, cross_validation,
                           number_of_folds, proportion_of_testing,
                           sample_proportion, type, importance)

  mtry <- floor(ncol(X) / 3)
  default_mtry <- c(ceiling(mtry / 2), mtry, mtry + ceiling(mtry / 2))
  if (is.null(mtry_theta)) {
    mtry_theta <- default_mtry
  }
  if (is.null(mtry_lambda)) {
    mtry_lambda <- default_mtry
  }

  zeros_percentaje <- round(sum(y == 0) / length(y) * 100, 1)

  if (zeros_percentaje < 20) {
    warning("It looks that your data has only ", zeros_percentaje, " percentaje ",
            "of zeros, in this case maybe ZAP random forest could not be your ",
            "best option.")
  }

  folds <- get_folds(cross_validation, number_of_folds, proportion_of_testing,
                     n_records=length(y))
  n_folds <- length(folds)
  loss_name <- ifelse(identical(loss_function, mse), "mse", "loss_function")

  all_combinations <- get_combinations(ntree_theta, mtry_theta,
                                       nodesize_theta, ntree_lambda,
                                       mtry_lambda, nodesize_lambda,
                                       sample_proportion, type)
  Combinations <- data.frame()
  n_combinations <- length(all_combinations)

  if (verbose && n_combinations > 1) {
    cat("\tTotal combinations: ", n_combinations, "\n")
  }

  for (i in 1:n_combinations) {
    if (verbose && n_combinations > 1) {
      cat("\t\tCombination:", i, "/", n_combinations, "\n")
    }
    flags <- all_combinations[[i]]

    inner_losses <- c()
    for (j in 1:n_folds) {
      if (verbose && n_combinations > 1) {
        cat("\t\t\tFold:", j, "/", n_folds, "\n")
      }
      fold <- folds[[j]]

      X_training <- X[fold$training, ]
      y_training <- y[fold$training]
      X_testing <- X[fold$testing, ]
      y_testing <- y[fold$testing]

      model <- suppressWarnings(zap.rfsrc(
        X_training, y_training, importance = FALSE,
        ntree_theta = flags$ntree_theta,
        mtry_theta = flags$mtry_theta,
        nodesize_theta = flags$nodesize_theta,
        ntree_lambda = flags$ntree_lambda,
        mtry_lambda = flags$mtry_lambda,
        nodesize_lambda = flags$nodesize_lambda))
      predictions <- predict(model, X_testing, type = flags$type)
      inner_losses <- c(inner_losses,
                        loss_function(actual = y_testing,
                                      predicted = predictions$predicted))
    }

    flags[[loss_name]] <- mean(inner_losses)
    Combinations <- rbind(Combinations, as.data.frame(flags))
  }

  model <- list(combinations = Combinations)
  index_best_params <- which.min(Combinations[[loss_name]])
  model$best_params <- as.list(Combinations[index_best_params, ])

  model$zap_fitted_model <- suppressWarnings(zap.rfsrc(
    X, y,
    importance = importance,
    ntree_theta = model$best_params$ntree_theta,
    mtry_theta = model$best_params$mtry_theta,
    nodesize_theta = model$best_params$nodesize_theta,
    ntree_lambda = model$best_params$ntree_lambda,
    mtry_lambda = model$best_params$mtry_lambda,
    nodesize_lambda = model$best_params$nodesize_lambda
  ))

  model <- structure(model, class = "tune.zap.rfsrc")

  return(model)
}

predict.tune.zap.rfsrc <- function(model, X) {
  return(predict(model$zap_fitted_model, X, type = model$best_params$type))
}