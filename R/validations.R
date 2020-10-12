expect_cross_validation <- function(cross_validation, number_of_folds,
                                    proportion_of_testing) {
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

expect_pheno <- function(Pheno) {
  expect_data_frame(Pheno, min.cols = 2)
  expect_subset(c("Line", "Env"), colnames(Pheno),
                empty.ok = FALSE,
                label = "Pheno",
                info = "Pheno does not contain Line, Hybrid_Name and/or Env column(s)")
}

expect_geno <- function(Geno, Markers) {
  if (is.null(Geno) && is.null(Markers)) {
    stop("You must provide the Geno or Markers matrix")
  } else if (!is.null(Geno) && !is.null(Markers)) {
    stop("You must provide only the Geno or Markers matrix, not both")
  }

  with_markers <- !is.null(Markers)
  if (is.null(Geno)) {
    Geno <- Markers
    Markers <- NULL
  }

  expect_data_frame(Geno)
  if (is.null(Geno$Line)) {
    name <- ifelse(with_markers, "Markers", "Geno")
    stop(name, " does not contains Line column")
  }
}

validate_zap_params <- function(X, y, ntree_theta, mtry_theta,
                                nodesize_theta, ntree_lambda, mtry_lambda,
                                nodesize_lambda, importance) {
  expect_data_frame(X)
  expect_numeric(y)

  if (anyNA(y)) {
    stop("Some elements in the response variable has NA, please remove them")
  } else if (!all(is_discrete(y))) {
    stop("ZAP random forest only can be used with discrete data")
  } else if (nrow(X) != length(y)) {
    stop("X must have the same number of rows that elements in y")
  }

  expect_number(ntree_theta, lower = 1)
  expect_number(mtry_theta, lower = 1, null.ok = TRUE)
  expect_number(nodesize_theta, lower = 1)

  expect_number(ntree_lambda, lower = 1)
  expect_number(mtry_lambda, lower = 1, null.ok = TRUE)
  expect_number(nodesize_lambda, lower = 1)

  expect_logical(importance)
}

validate_tune_zap_params <- function(X, y, ntree_theta, mtry_theta,
                                     nodesize_theta, ntree_lambda,
                                     mtry_lambda, nodesize_lambda,
                                     loss_function, cross_validation,
                                     number_of_folds, proportion_of_testing,
                                     sample_proportion, type, importance) {
  expect_data_frame(X)
  expect_numeric(y)

  if (anyNA(y)) {
    stop("Some elements in the response variable has NA, please removed them")
  } else if (!all(is_discrete(y))) {
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

  expect_function(loss_function)

  expect_cross_validation(cross_validation, number_of_folds,
                          proportion_of_testing)

  expect_numeric(sample_proportion, lower = 0.01, upper = 1)

  expect_character(type)
  type <- tolower(type)
  expect_subset(type, c("original", "custom"))

  expect_logical(importance, len = 1)
}

validate_gene_zap_params <- function(
  Pheno, y, Geno, Markers, with_interaction, mult_env_anal,
  ntree_theta, mtry_theta, nodesize_theta,
  ntree_lambda, mtry_lambda, nodesize_lambda,
  importance, type, loss_function,
  cross_validation, number_of_folds, proportion_of_testing,
  type_of_tuning, tuning_cross_validation, tuning_number_of_folds,
  tuning_proportion_of_testing, sample_proportion,
  results_dir, verbose, digits, seed) {
  expect_int(digits, lower = 1)
  expect_number(seed, null.ok = TRUE)
  expect_logical(verbose, len = 1)
  expect_logical(with_interaction, len = 1)
  expect_logical(mult_env_anal, len = 1)
  expect_string(results_dir)
  expect_string(type_of_tuning)
  type_of_tuning <- tolower(type_of_tuning)
  expect_subset(type_of_tuning, c("global", "local"))

  expect_pheno(Pheno)
  expect_geno(Geno, Markers)

  expect_cross_validation(cross_validation, number_of_folds,
                          proportion_of_testing)

  validate_tune_zap_params(Pheno, y, ntree_theta, mtry_theta=1,
                           nodesize_theta, ntree_lambda,
                           mtry_lambda=1, nodesize_lambda,
                           loss_function, tuning_cross_validation,
                           tuning_number_of_folds, tuning_proportion_of_testing,
                           sample_proportion, type, importance)

  expect_numeric(mtry_lambda, lower = 0.01, upper = 1, null.ok = TRUE)
  expect_numeric(mtry_theta, lower = 0.01, upper = 1, null.ok = TRUE)
}