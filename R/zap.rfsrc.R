g <- function(x, y) x / (1 - exp(-x)) - y

#' @title Zero Altered Poisson (ZAP) random forest
#'
#' @description
#' Train a Zero Altered Poisson (ZAP) ranfom forest model. It only works for
#' count data and it is recommended to use when such data has an excess of
#' zeros.
#'
#' @param X (\code{matrix} or \code{data.frame}) The predictors.
#' @param y (\code{numeric}) The count response variable.
#' @param ntree_theta (\code{numeric}) The number of trees used to train the
#'        binary random forest model. 500 by default.
#' @param mtry_theta (\code{numeric}) The number of independent variables that
#'        are going to be sampled in each iteration used to train the binary
#'        random forest model. ncol(X) / 3 by default.
#' @param nodesize_theta (\code{numeric}) The number of samples in the final
#'        nodes used to train the binary random forest model. 5 by default.
#' @param ntree_lambda (\code{numeric}) The number of trees used to train the
#'        numeric random forest model. 500 by default.
#' @param mtry_lambda (\code{numeric}) The number of independent variables that
#'        are going to be sampled in each iteration used to train the numeric
#'        random forest model. ncol(X) / 3 by default.
#' @param nodesize_lambda (\code{numeric}) The number of samples in the final
#'        nodes used to train the numeric random forest model. 5 by default.
#' @param importance (\code{logical}) A flag indicating if the importance of
#'        predictors have to be assessed.
#'
#' @return An object of S3 class "zap.rfsrc" which is a list with the both
#'         fitted random forest models theta_forest and lambda_forest.
#' @examples
#' \dontrun{
#' y <- rpois(nrow(new_iris), 10)
#' y[sample(nrow(new_iris), nrow(new_iris) * 0.2)] <- 0
#' X <- iris[, 1:4]
#'
#' model <- zap.rfsrc(X, y,
#'                    ntree_lambda=600, mtry_lambda=4, nodesize_lambda=8,
#'                    ntree_theta=600, mtry_theta=4, nodesize_theta=8)
#' predictions <- predict(model, X)$
#' predictions$predicted
#' }
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

  Data <- data.frame(X)
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