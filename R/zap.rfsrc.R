g <- function(x, y) x / (1 - exp(-x)) - y

zap.rfsrc <- function(X, y,
                      ntree_theta = 500,
                      mtry_theta = ncol(X) / 3,
                      nodesize_theta = 5,
                      ntree_lambda = 500,
                      mtry_lambda = ncol(X) / 3,
                      nodesize_lambda = 5,
                      importance = FALSE) {
  validate_zap_params(X, y, ntree_theta, mtry_theta,
                      nodesize_theta, ntree_lambda, mtry_lambda,
                      nodesize_lambda, importance)

  zeros_percentaje <- round(sum(y == 0) / length(y) * 100, 1)

  if (zeros_percentaje < 20) {
    warning("It looks that your data has only ", zeros_percentaje, " percentaje ",
            "of zeros, in this case maybe ZAP random forest could not be your ",
            "best option.")
  }

  Data <- as.data.frame(X)
  predictors <- colnames(Data)
  Data$y <- y
  Data$y0 <- Data$y == 0

  # Binary random forest
  theta_forest <- rfsrc(
    get_formula("y0", predictors),
    Data,
    ntree = ntree_theta,
    mtry = mtry_theta,
    nodesize = nodesize_theta,
    importance = importance
  )

  # For lambda include only the observations with values greater than zero.
  Data <- Data[Data$y > 0, ]

  lambda_forest <- rfsrc(
    get_formula("y", predictors),
    Data,
    ntree = ntree_lambda,
    mtry = mtry_lambda,
    nodesize = nodesize_lambda,
    splitrule = "custom1",
    importance = importance
  )

  model <- structure(list(theta_forest = theta_forest,
                          lambda_forest = lambda_forest),
                     class = "zap.rfsrc")

  return(model)
}

predict_theta <- function(theta_forest, X) {
  return(predict(theta_forest, newdata=X)$predicted)
}

predict_lambda <- function(lambda_forest, X) {
  yhat0 <- predict(lambda_forest, newdata=X)$predicted
  N <- nrow(X)

  lambda_predictions <- rep(0, N)

  for(i in 1:N) {
    lambda_predictions[i] <- uniroot(g, interval=c(-1, yhat0[i] + 2),
                                     y=yhat0[i])$root
    lambda_predictions[i] <- lambda_predictions[i] *
                             (lambda_predictions[i] > 0)
  }

  return(lambda_predictions)
}

predict.zap.rfsrc <- function(model, X, type="original") {
  thetas <- predict_theta(model$theta_forest, X)
  lambdas <- predict_lambda(model$lambda_forest, X)
  N <- nrow(X)

  if (type == "original") {
    predictions <- ((1 - thetas) * lambdas) / (1 - exp(-lambdas))
  } else if (type == "custom") {
    predictions <- rep(0, N)

    no_zero_predictions_indices <- thetas < 0.5
    predictions[no_zero_predictions_indices] <-
      lambdas[no_zero_predictions_indices]
  } else {
    stop("type parameter only can be original or custom, not", type)
  }

  return(list(predicted = round(predictions),
              lambda = round(lambdas),
              theta = round(thetas)))
}