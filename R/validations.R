validate_zap_params <- function(X, y, ntree_theta, mtry_theta,
                                nodesize_theta, ntree_lambda, mtry_lambda,
                                nodesize_lambda, importance) {
  expect_data_frame(X)
  expect_numeric(y)

  if (!all(is_discrete(y))) {
    stop("ZAP random forest only can be used with discrete data")
  } else if (nrow(X) != length(y)) {
    stop("X must have the same number of rows that elements in y")
  }

  expect_int(ntree_theta, lower = 1)
  expect_int(mtry_theta, lower = 1)
  expect_int(nodesize_theta, lower = 1)

  expect_int(ntree_lambda, lower = 1)
  expect_int(mtry_lambda, lower = 1)
  expect_int(nodesize_lambda, lower = 1)

  expect_logical(importance)
}

validate_tune_zap_params <- function(X, y, ntree_theta, mtry_theta,
                                     nodesize_theta, ntree_lambda,
                                     mtry_lambda, nodesize_lambda,
                                     loss_function, cross_validation,
                                     number_of_folds, proportion_of_testing) {
  expect_data_frame(X)
  expect_numeric(y)

  if (!all(is_discrete(y))) {
    stop("ZAP random forest only can be used with discrete data")
  } else if (nrow(X) != length(y)) {
    stop("X must have the same number of rows that elements in y")
  }

  expect_numeric(ntree_theta, lower = 1)
  expect_numeric(mtry_theta, lower = 1, null.ok = TRUE)
  expect_numeric(nodesize_theta, lower = 1)

  expect_numeric(ntree_lambda, lower = 1)
  expect_numeric(mtry_lambda, lower = 1, null.ok = TRUE)
  expect_numeric(nodesize_lambda, lower = 1)

  expect_logical(importance)

  expect_function(loss_function)

  expect_string(cross_validation)
  cross_validation <- tolower(cross_validation)
  expect_subset(cross_validation, c("k_fold", "random_partition"))
  if (cross_validation == "k_fold") {
    expect_int(number_of_folds, lower = 2)
  } else if (cross_validation == "random_partition") {
    expect_int(number_of_folds, lower = 1)
    expect_numeric(proportion_of_testing, lower = 0.01, upper = 0.99)
  }
}